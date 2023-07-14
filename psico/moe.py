"""
(c) 2022 Thomas Holder

License: BSD-2-Clause
"""

from pymol import cmd


def dict_from_tag_value(items: list) -> dict:
    return {item["tag"]: item["value"] for item in items}


def parse_str(s: str) -> str:
    return s.replace("!", " ")


TYPES = {
    'i': int,
    'ix': lambda s: int(s, 16),
    'r': float,
    't': parse_str,  # token or text?
    'tt': parse_str,  # ?
    'c': parse_str,  # character?
    's': parse_str,  # string?
}


def typedvalue(typ: str, value):
    try:
        return TYPES[typ](value)
    except ValueError as ex:
        raise ValueError(f"typ={typ} {ex}")


class MoeParserError(Exception):
    pass


class MoeParser:

    def parse(self, contents: bytes, filename: str = "<buffer>") -> list:
        """
        Args:
          contents: File content string
          filename: File path for error reporting

        Return:
          List of data blocks
        """
        self.filename = filename
        self.linenumber = 0

        if isinstance(contents, bytes):
            contents = contents.decode("utf-8", errors="replace")

        if isinstance(contents, str):
            return self.parse_lines(contents.splitlines())

        assert contents is None

        with open(filename) as handle:
            return self.parse_lines(handle)

    def make_parser_exception(self, msg: str) -> MoeParserError:
        return MoeParserError(msg, f"{self.filename}:{self.linenumber}")

    def gen_lines(self, line_it: iter) -> iter:
        for line in line_it:
            self.linenumber += 1
            yield line

    def gen_words(self) -> iter:
        while True:
            yield from next(self.linehandle).split()

    def parse_attr(self, headerline: str) -> dict:
        a = headerline.split()
        assert a[0].startswith("#")
        return {
            "key": a[0][1:],
            "count": int(a[1]),
            "cols": list(zip(a[2::2], a[3::2])),
        }

    def parse_concat_value(self) -> str:
        count = int(next(self.wordhandle))
        return "".join(next(self.wordhandle) for _ in range(count))

    def parse_array(self, typ="*") -> list:
        count = int(next(self.wordhandle))
        return self.parse_fixed_array(count, typ)

    def parse_fixed_array(self, count: int, typ="*") -> list:
        return [self.parse_value(typ) for _ in range(count)]

    def parse_value(self, typ: str):
        if typ == "*":
            typ = next(self.wordhandle)

        if typ in TYPES:
            value = next(self.wordhandle)
            if value == "~":
                value = self.parse_concat_value()
            return typedvalue(typ, value)

        if typ[1:] == "*":
            return self.parse_array(typ[0])

        if typ[1:].isdigit():
            return self.parse_fixed_array(int(typ[1:]), typ[0])

        raise self.make_parser_exception(typ)

    def parse_values(self, attr: dict) -> iter:
        for key, typ in attr["cols"]:
            typ, delim, default = typ.partition("=")
            if delim:
                value = typedvalue(typ, default)
            else:
                value = self.parse_value(typ)
            yield key, value

    def parse_values_block(self, attr: dict) -> iter:
        for _ in range(attr["count"]):
            yield dict(self.parse_values(attr))

    def parse_lines_raw(self, line_it: iter) -> list:
        self.linehandle = self.gen_lines(line_it)
        self.wordhandle = self.gen_words()
        blocks = []

        line = next(self.linehandle)
        assert line.startswith("#moe")

        for line in self.linehandle:
            assert line.startswith("#")
            if line.startswith("#end"):
                continue

            attr = self.parse_attr(line)
            values = list(self.parse_values_block(attr))

            blocks.append((attr["key"], values))

        return blocks

    def parse_lines(self, line_it: iter) -> list:
        try:
            return self.parse_lines_raw(line_it)
        except MoeParserError:
            raise
        except Exception as ex:
            raise self.make_parser_exception(str(ex))


ATOM_1to1 = {
    "aName": "name",
    "aElement": "symbol",
    "aTempFactor": "b",
    "aOccupancy": "q",
    "aIon": "formal_charge",
    "aCharge": "partial_charge",
}

RESIDUE_1to1 = {
    "rName": "resn",
    "rUID": "resi_number",
    "rINS": "ins_code",
}


@cmd.extend
def load_moe(filename: str,
             oname: str,
             state: int = 0,
             *,
             contents: bytes = None,
             partial: int = 0,
             discrete: int = -1,
             quiet: int = 1,
             zoom: int = -1,
             _self=cmd):
    """
DESCRIPTION

    Load a molecular object from an MOE file.

ARGUMENTS

    partial = 0/1: If 1, then don't create "set" selections {default: 0}
    """
    import chempy
    from chempy import models

    molecule = None
    residues = {}
    chains = {}
    rEnd = 0
    systemdata = {}
    atom_ranges = {}
    collections = []

    blocks = MoeParser().parse(contents, filename)

    for blockname, blockvalues in blocks:
        if blockname == "system":
            systemdata = dict_from_tag_value(blockvalues)

        elif blockname == "molecule":
            moleculedata = dict_from_tag_value(blockvalues)

            molecule = models.Indexed()
            molecule.atom.extend(chempy.Atom()
                                 for _ in range(moleculedata["atoms"]))

        elif blockname == "bond":
            assert molecule is not None

            for bonddata in blockvalues:
                bond = chempy.Bond()
                bond.index = [bonddata["a"] - 1, bonddata["b"] - 1]
                bond.order = bonddata.get("o", 1)
                molecule.bond.append(bond)

        elif blockname == "attr":
            assert molecule is not None

            for row in blockvalues:
                if "rAtomCount" in row:
                    rStart = rEnd
                    rEnd += row["rAtomCount"]
                    residues[row["ID"]] = molecule.atom[rStart:rEnd]
                    atom_ranges[row["ID"]] = (rStart, rEnd)

                if "cResidueCount" in row:
                    chains[row["ID"]] = row

                for key, value in row.items():
                    if key in ("rINS", ):
                        value = value.strip()

                    if key in ATOM_1to1:
                        setattr(molecule.atom[row["ID"] - 1], ATOM_1to1[key],
                                value)
                    elif key == "aPosX":
                        molecule.atom[row["ID"] - 1].coord = [
                            value,
                            row["aPosY"],
                            row["aPosZ"],
                        ]
                    elif key in RESIDUE_1to1:
                        for atom in residues[row["ID"]]:
                            setattr(atom, RESIDUE_1to1[key], value)
                    elif key == "rType":
                        hetatm = value not in ("amino", "rna")
                        for atom in residues[row["ID"]]:
                            atom.hetatm = hetatm

        elif blockname == "collection":
            collections = blockvalues

    assert molecule is not None
    assert len(residues) > 0
    resiter = iter(residues.values())
    cAtomEnd = 0

    for chainindex, row in enumerate(chains.values()):
        chain_name = _self.get_legal_name(row["cName"]) + f"_{chainindex}"
        segi_name = _self.get_legal_name(row["cTag"])
        cAtomBegin = cAtomEnd

        for _ in range(row["cResidueCount"]):
            for atom in next(resiter):
                atom.chain = chain_name
                atom.segi = segi_name

                if row.get("cRGB", 0):
                    atom.trgb = row["cRGB"]

                cAtomEnd += 1

            atom_ranges[row['ID']] = (cAtomBegin, cAtomEnd)

    if "CellParameters" in systemdata:
        molecule.spacegroup, dims, angles = systemdata["CellParameters"]
        molecule.cell = dims + angles

    _self.load_model(molecule,
                     oname,
                     state,
                     discrete=int(discrete) > 0,
                     quiet=quiet,
                     zoom=zoom)

    _self.set("retain_order", 1, oname)
    _self.alter(oname, "elec_radius = vdw * 1.2")

    if int(partial):
        return

    selections = set(_self.get_names("selections"))
    tmpname = _self.get_unused_name("_tmp_sele_merge_or_rename")

    for row in collections:
        id_list = []
        for ID in row["objects"]:
            if ID < len(molecule.atom):
                id_list.append(ID - 1)
            elif ID in atom_ranges:
                id_list.extend(range(*atom_ranges[ID]))

        _self.select_list(tmpname, oname, id_list, mode="rank")

        # Merge existing selections, but don't overwrite existing objects
        if row["name"] in selections:
            _self.select(row["name"], tmpname, merge=True)
            _self.delete(tmpname)
        else:
            _self.set_name(
                tmpname, _self.get_unused_name(row["name"],
                                               alwaysnumber=False))


# vi:expandtab:smarttab
