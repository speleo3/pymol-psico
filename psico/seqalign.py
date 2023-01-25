'''
API only sequence alignment stuff.

Only load on demand to prevent unnecessary "Bio" import (biopython).

(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

import functools
from Bio.SeqIO import _FormatToIterator
from pymol import CmdException

def _assert_package_import():
    if not __name__.endswith('.seqalign'):
        raise CmdException("Must do 'import psico.seqalign' instead of 'run ...'")


DEFAULT_GAPOPEN_EMBOSS = 10
DEFAULT_GAPEXTEND_EMBOSS = 0.5

DEFAULT_GAPOPEN_BLASTP = 11
DEFAULT_GAPEXTEND_BLASTP = 1

DEFAULT_GAPOPEN = DEFAULT_GAPOPEN_EMBOSS
DEFAULT_GAPEXTEND = DEFAULT_GAPEXTEND_EMBOSS


@functools.lru_cache(maxsize=None)
def _get_aligner_BLOSUM62() -> "Bio.Align.PairwiseAligner":
    from Bio.Align import PairwiseAligner, substitution_matrices
    blosum62 = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner(internal_open_gap_score=-DEFAULT_GAPOPEN,
                              extend_gap_score=-DEFAULT_GAPEXTEND,
                              substitution_matrix=blosum62)
    assert aligner.mode == "global"
    return aligner


def _msa_from_pairwise(
    pairwise: "Bio.Align.PairwiseAlignment"
) -> "Bio.Align.MultipleSeqAlignment":
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    seqs = pairwise.format().splitlines()[::2]
    return MultipleSeqAlignment([
        SeqRecord(Seq(seqs[0]), "s1"),
        SeqRecord(Seq(seqs[1]), "s2"),
    ])


def needle_alignment(s1, s2):
    '''
DESCRIPTION

    Does a Needleman-Wunsch Alignment of sequence s1 and s2 and
    returns a Bio.Align.MultipleSeqAlignment object.
    '''
    try:
        # requires biopython >= 1.72
        aligner = _get_aligner_BLOSUM62()
    except ImportError as ex:
        print(f" needle_alignment-Warning: {ex}")
        return needle_alignment_Bio_pairwise2(s1, s2)
    alns = aligner.align(s1, s2)
    return _msa_from_pairwise(alns[0])


def needle_alignment_Bio_pairwise2(s1, s2):
    '''
    Same as `needle_alignment` but uses deprecated Bio.pairwise2 module.
    '''
    from Bio import pairwise2
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    try:
        from Bio.Align import substitution_matrices
    except ImportError:
        from Bio.SubsMat.MatrixInfo import blosum62
    else:
        blosum62 = substitution_matrices.load("BLOSUM62")

    def match_callback(c1, c2):
        return blosum62.get((c1, c2), 1 if c1 == c2 else -4)

    alns = pairwise2.align.globalcs(s1, s2,
            match_callback, -10., -.5,
            one_alignment_only=True)

    return MultipleSeqAlignment([
        SeqRecord(Seq(alns[0][0]), "s1"),
        SeqRecord(Seq(alns[0][1]), "s2"),
    ])

def needle_alignment_emboss(s1, s2):
    '''
    Same as `needle_alignment` but uses EMBOSS' "needle" executable.
    '''
    import subprocess
    from Bio.Emboss.Applications import NeedleCommandline
    from Bio import AlignIO
    cline = NeedleCommandline(auto=True, sprotein=True, stdout=True, gapopen=10, gapextend=1)
    cline.asequence = "asis:" + s1
    cline.bsequence = "asis:" + s2
    process = subprocess.Popen(str(cline), shell=True, stdout=subprocess.PIPE,
            universal_newlines=True)
    return AlignIO.read(process.stdout, "emboss")

def alignment_mapping(seq1, seq2):
    '''
DESCRIPTION

    Returns an iterator with seq1 indices mapped to seq2 indices

    >>> mapping = dict(alignment_mapping(s1, s2))
    '''
    i, j = -1, -1
    for a, b in zip(seq1, seq2):
        if a != '-': i += 1
        if b != '-': j += 1
        if a != '-' and b != '-': yield i, j

def aln_magic_format(infile):
    '''
DESCRIPTION

    Guess alignment file format.
    '''
    with open(infile) as handle:
        for line in handle:
            if len(line.rstrip()) > 0:
                break
    if line.startswith('CLUSTAL') or line.startswith('MUSCLE'):
        informat = 'clustal'
    elif line.startswith('>P1;'):
        informat = 'pir'
    elif line.startswith('>'):
        informat = 'fasta'
    elif line.startswith('# STOCKHOLM'):
        informat = 'stockholm'
    elif line.startswith('Align '):
        informat = 'fatcat'
    elif line.startswith('# ProSMART Alignment File'):
        informat = 'prosmart'
    else:
        informat = 'emboss'
    return informat

def aln_magic_read(infile, format=None):
    '''
DESCRIPTION

    Wrapper for Bio.AlignIO.read that guesses alignment file format.
    '''
    from Bio import AlignIO
    if not format:
        format = aln_magic_format(infile)
    with open(infile) as handle:
        return AlignIO.read(handle, format)

def FatCatIterator(handle):
    '''
DESCRIPTION

    Biopython FatCat alignment file support
    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    try:
        line = handle.readline()
        a = line.split()
        ids = a[1], a[4]
        ids = [(id[:-4] if id.endswith('.pdb') else id) for id in ids]
    except ImportError:
        ids = 'chain1', 'chain2'

    seqs = [[], []]
    for line in handle:
        if line.startswith('Chain 1:'): seqs[0].append(line[14:].rstrip())
        if line.startswith('Chain 2:'): seqs[1].append(line[14:].rstrip())

    for seq, id in zip(seqs, ids):
        yield SeqRecord(Seq(''.join(seq)), id=id, name=id)


def ProSMARTIterator(handle):
    '''
DESCRIPTION

    Biopython ProSMART alignment file support

    Warning: ProSMART alignment files do not contain gaps, so they only
    contain the residues which are actually aligned!
    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    _assert_package_import()
    from . import one_letter

    try:
        code1, chain1, code2, chain2 = handle.name.rsplit('.', 1)[0].split('_')
        ids = [code1 + '_' + chain1, code2 + '_' + chain2]
    except:
        ids = 'chain1', 'chain2'

    seqs = [[], []]
    for line in handle:
        if line.startswith('#') or line.startswith('Res1'):
            continue
        a = line.split('\t')
        try:
            aa1 = one_letter[a[2]]
            aa2 = one_letter[a[3]]
        except:
            print(' Warning: Exception while parsing ProSMART alignment')
            continue
        seqs[0].append(aa1)
        seqs[1].append(aa2)

    for seq, id in zip(seqs, ids):
        yield SeqRecord(Seq(''.join(seq)), id=id, name=id)


def POAIterator(handle):
    '''
DESCRIPTION

    Biopython POA (POSA) alignment file support

    Warning: POA alignment files do not contain gaps, so they only
    contain the residues which are actually aligned!
    '''
    import re
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    pattern = re.compile(r'\{ *\d+\}')

    for line in handle:
        if line.startswith('#'):
            continue
        id, seq = line.split(None, 1)
        seq = pattern.sub('-', seq.rstrip())
        yield SeqRecord(Seq(''.join(seq)), id=id, name=id)


if 'fatcat' not in _FormatToIterator:
    _FormatToIterator['fatcat'] = FatCatIterator

if 'prosmart' not in _FormatToIterator:
    _FormatToIterator['prosmart'] = ProSMARTIterator

if 'poa' not in _FormatToIterator:
    _FormatToIterator['poa'] = POAIterator

# vi:expandtab:smarttab
