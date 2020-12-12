import os
import psico.seqalign

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
FILENAME_FATCAT = os.path.join(DATA_DIR, '1ubqA-1hezE.fatcat')


def test_needle_alignment():
    A = psico.seqalign.needle_alignment("ACDEFGHIKLMN", "DEYGHVVVVIKLMN")
    assert str(A[0].seq) == "ACDEFGH----IKLMN"
    assert str(A[1].seq) == "--DEYGHVVVVIKLMN"


def test_alignment_mapping():
    seq1 = "ACDEFGH----IKLMN"
    seq2 = "--DEYGHVVVVIKLMN"
    mapping = psico.seqalign.alignment_mapping(seq1, seq2)
    assert dict(mapping) == {
        2: 0,
        3: 1,
        4: 2,
        5: 3,
        6: 4,
        7: 9,
        8: 10,
        9: 11,
        10: 12,
        11: 13,
    }


def test_aln_magic_format():
    assert "fatcat" == psico.seqalign.aln_magic_format(FILENAME_FATCAT)


def test_aln_magic_read():
    # indirectly also tests FatCatIterator
    A = psico.seqalign.aln_magic_read(FILENAME_FATCAT)
    S = [
        "IFVKTLTGKTITLEVEPSDT--IENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVL",
        "VNLIFADGKIQTAEFKGTFEEATAEAYRYADLLAKVNGEYTADLED-----------GGNH-----MNIKF",
    ]
    assert str(A[0].seq) == S[0]
    assert str(A[1].seq) == S[1]


# def needle_alignment_emboss(s1, s2):
# def FatCatIterator(handle):
# def ProSMARTIterator(handle):
# def POAIterator(handle):
