"""
Implementation of AGGRESCAN according to Conchillo-Solé et al. 2007 and the
corresponding web tool at http://bioinf.uab.es/aap/

Implementation details are tweaked to match the web tool results and may
devitate from the paper.

(c) 2022 Thomas Holder
(c) 2022 Leukocare AG

License: BSD-2-Clause
"""

from pymol import cmd, CmdException

# Amino-acid aggregation-propensity value (a3v)
# N. Sánchez de Groot et al.
# Prediction of "hot spots" of aggregation in disease-linked polypeptides.
# BMC Structural Biology 2005, 5:18
AGGRESCAN_A3V = {
    "I": 1.822,
    "F": 1.754,
    "V": 1.594,
    "L": 1.380,
    "Y": 1.159,
    "W": 1.037,
    "M": 0.910,
    "C": 0.604,
    "A": -0.036,
    "T": -0.159,
    "S": -0.294,
    "P": -0.334,
    "G": -0.535,
    "K": -0.931,
    "H": -1.033,
    "Q": -1.231,
    "R": -1.240,
    "N": -1.302,
    "E": -1.412,
    "D": -1.836,
}

# termini get average scores of charged residues, according to publication
AGGRESCAN_A3V["^"] = (AGGRESCAN_A3V["K"] + AGGRESCAN_A3V["R"]) / 2
AGGRESCAN_A3V["$"] = (AGGRESCAN_A3V["D"] + AGGRESCAN_A3V["E"]) / 2

# FIXME values in web tool appear to be swapped!
AGGRESCAN_A3V["^"], AGGRESCAN_A3V["$"] = AGGRESCAN_A3V["$"], AGGRESCAN_A3V["^"]

# Hot Spot Threshold (HST)
AGGRESCAN_HST = -0.02

# Minimum number of residues per hotspot
HS_MIN_RES_COUNT = 5


def mean(values: list) -> float:
    """Compute the arithmetic mean"""
    return sum(values) / len(values)


def trapezoidal_integration(values: list) -> float:
    """Trapezoidal integration (midpoint rule)"""
    return sum(values[1:-1]) + (values[0] + values[-1]) / 2


def aggrescan1d_scores(seq: str, winsize: int = None) -> dict:
    """Aggrescan reimplementation.

    Equivalent to http://bioinf.uab.es/aggrescan/

    Conchillo-Solé et al. AGGRESCAN: a server for the prediction and evaluation
    of "hot spots" of aggregation in polypeptides. BMC Bioinformatics 8, 65
    (2007). https://doi.org/10.1186/1471-2105-8-65

    Args:
      seq: Amino acid sequence in single letter code.
      winsize: Window size, must be odd and larger than 2.

    Returns:
      Dictionary with aggregation and hotspot scores.

    Dictionary keys:
      - a3vSA: a3v Sequence Average
      - nHS: Number of Hot Spots
      - NnHS: Normalized number of Hot Spots for 100 residues
      - AAT: Area of the Aggregation Profile above the hot-spot Threshold
      - THSA: Total Hot-Spot Area
      - TA: Total Area of the aggregation profile
      - AATr: AAT per residue
      - THSAr: THSA per residue
      - a4vSS: a4v Sequence Sum
      - Na4vSS: Normalized a4vSS for 100 residues
      - a3v: amino-acid aggregation-propensity value
      - a4v: a3v window average
      - HSA: Hot-Spot Area
      - NHSA: Normalized Hot-Spot Area per residue
      - a4vAHS: a4v average in the Hot Spot
    """
    assert isinstance(seq, str)
    non_aa = set(seq).difference("ACDEFGHIKLMNPQRSTVWY")
    if non_aa:
        print(f"Warning: unexpected codes in sequence: {non_aa}")

    N = len(seq)

    if winsize is None:
        if N <= 75:
            winsize = 5
        elif N <= 175:
            winsize = 7
        elif N <= 300:
            winsize = 9
        else:
            winsize = 11

    assert winsize > 2 and winsize % 2 == 1

    # amino-acid aggregation-propensity value (a3v)
    a3v = [AGGRESCAN_A3V.get(aa, AGGRESCAN_A3V["A"]) for aa in seq]
    a3v_with_termini = [AGGRESCAN_A3V["^"]] + a3v + [AGGRESCAN_A3V["$"]]
    a4v = [
        mean(a3v_with_termini[i:i + winsize]) for i in range(N - winsize + 3)
    ]

    # "The remaining off-centre N- and C-terminal residues are assigned the a4v
    # calculated for the first and last window centres, respectively."
    offcenter = winsize // 2 - 1
    a4v = [a4v[0]] * offcenter + a4v + [a4v[-1]] * offcenter

    assert len(a3v) == N
    assert len(a4v) == N

    a4v_minus_hst = [(v - AGGRESCAN_HST) for v in a4v]

    nHS = 0  # Number of Hot Spots (nHS)
    HSA = [0.0] * N  # Hot-Spot Area (HSA)
    NHSA = [0.0] * N  # Normalized Hot-Spot Area per residue
    a4vAHS = [0.0] * N  # a4v average in the Hot Spot
    run_hs = 0

    for (j, aa) in enumerate(seq + "P"):
        if aa != "P" and a4v_minus_hst[j] > 0:
            run_hs += 1
            continue

        if run_hs > 0:
            i = j - run_hs
            hsa = sum(a4v_minus_hst[i:j])  # trapezoidal_integration?

            # FIXME Needed to match values in web tool
            if i == 0:
                hsa += a4v_minus_hst[0] / 2

            HSA[i:j] = [hsa] * run_hs
            a4vAHS[i:j] = [mean(a4v[i:j])] * run_hs

            if run_hs >= HS_MIN_RES_COUNT:
                NHSA[i:j] = [hsa / run_hs] * run_hs
                nHS += 1

        run_hs = 0

    ta = trapezoidal_integration(a4v_minus_hst)
    aat = sum(max(0, v) for v in a4v_minus_hst)  # trapezoidal_integration?
    thsa = sum(NHSA)
    a4vSS = sum(a4v)

    return {
        "a3vSA": mean(a3v),
        "nHS": nHS,
        "NnHS": nHS / N * 100,
        "AAT": aat,
        "THSA": thsa,
        "TA": ta,
        "AATr": aat / N,
        "THSAr": thsa / N,
        "a4vSS": a4vSS,
        "Na4vSS": a4vSS / N * 100,
        "table": {
            "a3v": a3v,
            "a4v": a4v,
            "HSA": HSA,
            "NHSA": NHSA,
            "a4vAHS": a4vAHS,
        },
    }


@cmd.extend
def aggrescan1d(selection="polymer",
                key="a4v",
                palette="white forest",
                minimum=0,
                maximum=None,
                *,
                var='b',
                quiet=1,
                _self=cmd):
    """
DESCRIPTION

    Aggrescan (sequence based).

ARGUMENTS

    key = a3v | a4v | HSA | NHSA | a4vAHS: Score used for coloring. For the raw
    aggregation profile, use a4v. For hotspots, use NHSA. {default: a4v}

    var = b | q | ...: PyMOL property to assign the score to {default: b}
    """
    from psico.querying import iterate_to_list

    if maximum is None:
        maximum = {
            "a3v": 1.8,
            "a4v": 0.7,
            "HSA": 7,
            "NHSA": 0.5,
            "a4vAHS": 0.5,
        }.get(key, 5)

    ids_key = "(model, segi, chain), resi, oneletter"
    ids_all = iterate_to_list(f"({selection}) & guide", ids_key)
    data = {}

    def gen_chain_ranges():
        i, prev = -1, ()
        for j, idtuple in enumerate(ids_all + [(None, )]):
            current = idtuple[0]
            if prev != current:
                if i != -1:
                    yield ids_all[i:j]
                i, prev = j, current

    for ids in gen_chain_ranges():
        seq = "".join(idtuple[-1] for idtuple in ids)
        agg = aggrescan1d_scores(seq)
        table = agg.pop("table")
        assert len(ids) == len(table[key])
        data.update(zip(ids, table[key]))

        if not int(quiet):
            print('/{0}/{1}/{2}'.format(*ids[0][0]))

            for akey, values in table.items():
                print(f'Max {akey:7s}: {max(values):7.3f}')

            for akey, value in agg.items():
                if isinstance(value, float):
                    print(f"{akey:7s}: {value:7.3f}")
                else:
                    print(f"{akey:7s}: {value:7}")

    _self.alter(selection,
                f"{var} = data.get(({ids_key}), 0)",
                space={"data": data})

    _self.spectrum(var,
                   palette,
                   selection,
                   minimum=minimum,
                   maximum=maximum,
                   quiet=quiet)


@cmd.extend
def aggrescan3d(selection: str = "polymer",
                palette: str = "white forest",
                minimum: float = 0,
                maximum: float = 2,
                *,
                ph: float = 7.0,
                missing: float = 0.0,
                color_missing: str = "yellow",
                state: int = -1,
                var: str = 'b',
                quiet: int = 1,
                exe: str = "aggrescan",
                _self=cmd):
    """
DESCRIPTION

    Aggrescan3d (structure based).
    """
    quiet = int(quiet)

    import os, tempfile, shutil, subprocess, csv
    from psico.querying import iterate_to_sele

    ph_arg = ["--ph", str(ph)]

    # sanity check
    for chain in _self.get_chains(selection):
        if len(chain) > 1:
            raise CmdException(f"Chain ID longer than 1: {chain!r}")

    # temporary directory
    tempdir = tempfile.mkdtemp()
    if not quiet:
        print(' Tempdir:', tempdir)

    # filenames
    pdbfile = os.path.join(tempdir, 'my-input.pdb')

    try:
        _self.save(pdbfile, selection, state=state)

        subprocess.check_call(
            [exe, "--protein", pdbfile, "--work_dir", tempdir] + ph_arg)

        with open(f"{tempdir}/A3D.csv") as handle:
            scores = {(row['chain'], str(row['residue'])): float(row["score"])
                      for row in csv.DictReader(handle)}
    finally:
        shutil.rmtree(tempdir)

    if not quiet:
        sum_scores = sum(scores.values())
        sum_scores_pos = sum(s for s in scores.values() if s > 0)

        print(" Aggrescan3D"
              f" sum(all): {sum_scores:.2f}"
              f" sum(positive): {sum_scores_pos:.2f}")

    _self.alter(selection,
                f"{var} = scores.get((chain, resi), {missing!r})",
                space={"scores": scores})

    _self.spectrum(var, palette, selection, minimum=minimum, maximum=maximum)

    sele_missing = iterate_to_sele(selection,
                                   "(chain, resi) not in scores",
                                   space={"scores": scores})
    if sele_missing:
        _self.color(color_missing, sele_missing)


# tab-completion of arguments

cmd.auto_arg[0].update([
    ('aggrescan1d', cmd.auto_arg[1]['select']),
    ('aggrescan3d', cmd.auto_arg[1]['select']),
])

# vi:expandtab:smarttab
