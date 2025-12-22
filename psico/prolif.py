"""
(c) 2025 Thomas Holder

License: BSD-2-Clause
"""

import warnings

from pymol import cmd

from psico.querying import iterate_state_to_list

warnings.filterwarnings("ignore", message="pkg_resources is deprecated", category=UserWarning, module="prolif")
warnings.filterwarnings("ignore", message="tables has been moved", category=DeprecationWarning, module="MDAnalysis")

CURRENT_STATE = -1

CONTACT_COLORS = {
    "VdWContact": "magenta",
    "Hydrophobic": "forest",
}


@cmd.extendaa(
    cmd.auto_arg[1]["distance"],
    cmd.auto_arg[2]["distance"],
)
def prolif_interactions(
    sele_lig: str,
    sele_pro: str = "polymer",
    *,
    state: int = CURRENT_STATE,
    prefix: str = "prolif.",
    prefix_pseudo: str = "pseudo.",
    rep_pro: str = "licorice",
    quiet: int = 0,
    _self=cmd,
):
    """
DESCRIPTION

    Find interactions with ProLIF

    https://prolif.readthedocs.io/

USAGE

    prolif_interactions sele_lig [, sele_pro]

ARGUMENTS

    sele_lig = str: Ligand selection

    sele_pro = str: Protein selection (default: polymer)

    rep_pro = str: Atom representation for contacting protein residues (default: licorice)

EXAMPLE

    fetch 1eve
    h_add
    prolif_interactions resn E20
    """
    import prolif as plf
    from rdkit import Chem

    idx_pro = iterate_state_to_list(state, sele_pro, 'f"{model}`{index}"', _self=_self)
    idx_lig = iterate_state_to_list(state, sele_lig, 'f"{model}`{index}"', _self=_self)

    rdkit_mol_pro = Chem.MolFromPDBBlock(_self.get_str("pdb", sele_pro, state), removeHs=False, sanitize=False)
    rdkit_mol_lig = Chem.MolFromMolBlock(_self.get_str("mol", sele_lig, state), removeHs=False, sanitize=False)

    for rdkit_mol in [rdkit_mol_pro, rdkit_mol_lig]:
        failed_flag = Chem.SanitizeMol(rdkit_mol, catchErrors=True)
        if failed_flag:
            print(f"Sanitize failed for {rdkit_mol} on flag {failed_flag}")

    fp = plf.Fingerprint()
    fp.run_from_iterable([plf.Molecule(rdkit_mol_lig)], plf.Molecule(rdkit_mol_pro), n_jobs=1)

    if prefix.endswith("."):
        _self.group(prefix.removesuffix("."))

    if prefix_pseudo.endswith("."):
        _self.group(prefix_pseudo.removesuffix("."))

    assert len(fp.ifp) == 1

    i_pro_all = set()

    for internum, interactions in enumerate(fp.ifp[0].values(), 1):
        for intertype, interaction_tuple in interactions.items():
            for interobj in interaction_tuple:
                indices = interobj["parent_indices"]

                if not int(quiet):
                    print(f"{intertype:11s} {interobj['distance']:.2f} {indices}")

                i_lig = indices["ligand"]
                i_pro = indices["protein"]
                s_pro = " ".join(idx_pro[i] for i in i_pro)
                s_lig = " ".join(idx_lig[i] for i in i_lig)

                i_pro_all.update(i_pro)

                if len(i_lig) > 1:
                    s_lig_pseudo = f"{prefix_pseudo}{internum}.lig"
                    _self.pseudoatom(s_lig_pseudo, s_lig, state=state)
                    s_lig = s_lig_pseudo

                if len(i_pro) > 1:
                    s_pro_pseudo = f"{prefix_pseudo}{internum}.pro"
                    _self.pseudoatom(s_pro_pseudo, s_pro, state=state)
                    s_pro = s_pro_pseudo

                objname = f"{prefix}{intertype}"
                _self.distance(objname, s_pro, s_lig)

                color = CONTACT_COLORS.get(intertype)
                if color:
                    _self.color(color, objname)

    if rep_pro:
        s_pro = " ".join(idx_pro[i] for i in i_pro_all)
        _self.show(rep_pro, f"byres ({s_pro})")
