'''
Show UniProt-SNPs within PyMOL
Show dbSNP-SNPs within PyMOL

Commands: snp_uniprot, snp_ncbi

(c) 2010 Thomas Holder

License: BSD-2-Clause
'''

if not __name__.endswith('.snp'):
    raise Exception("Must do 'import psico.snp' instead of 'run ...'")

from pymol import cmd

# Email for Entrez connections
emailaddress = "pymol@psico.snp"


def snp_common(record, selection, label, name, quiet, *, _self=cmd):
    '''
    Common part of snp_uniprot and snp_ncbi.
    Argument `record' must be a Bio.SwissProt.Record object with `sequence',
    `entry_name' and `features' fields defined.
    '''
    from . import one_letter
    from .seqalign import needle_alignment, alignment_mapping

    label = int(label)
    quiet = int(quiet)
    pdbids = _self.get_object_list(selection)
    chains = _self.get_chains(selection)

    if len(pdbids) != 1:
        print('please select one object')
        return

    snpi = set()
    snpi_str = []
    labels = dict()

    for chain in chains:
        print('chain ' + chain)
        res_list = []
        _self.iterate('(%s) and chain %s and name CA' % (selection, chain),
                'res_list.append((resn,resv))', space=locals())
        seq = ''.join([one_letter.get(res[0], 'X') for res in res_list])
        align = needle_alignment(record.sequence, seq)
        if not quiet:
            align._records[0].id = record.entry_name
            align._records[1].id = pdbids[0] + '_' + chain
            print(align.format('clustal'))
        map1 = dict(alignment_mapping(*align))
        for feature in record.features:
            if feature[0] != 'VARIANT' or feature[1] != feature[2]:
                continue
            i = feature[1]
            if (i - 1) not in map1:
                if not quiet:
                    print('not mapped', feature)
                continue
            resi = res_list[map1[i - 1]][1]
            snpi.add(resi)
            if not quiet:
                print('%s`%d' % res_list[map1[i - 1]], feature[2:4])
            if label:
                labels.setdefault((chain, resi), []).append(feature[3].split(' (')[0])
        if len(snpi) > 0:
            snpi_str.append('(chain %s and resi %s)' % (chain, '+'.join(map(str, snpi))))

    for chain, resi in labels:
        lab = ', '.join(labels[(chain, resi)])
        _self.label('(%s) and chain %s and resi %d and name CA' % (selection, chain, resi), repr(lab))

    if len(snpi_str) == 0:
        print('no missense variants')
        return

    if name == '':
        name = _self.get_unused_name('nsSNPs')
    _self.select(name, '(%s) and (%s)' % (selection, ' or '.join(snpi_str)))


def snp_uniprot(uniprotname, selection='(all)', label=1, name='', quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    Selects all UniProt annotated nsSNPs (natural variants) in given
    structure. Does a sequence alignment of UniProt sequence and PDB
    sequence.

USAGE

    snp_uniprot uniprotname [, selection [, label [, name [, quiet]]]]

ARGUMENTS

    uniprotname = string: UniProt reference (like HBB_HUMAN or P68871)

    selection = string: atom selection

    label = 0 or 1: Label CA atoms of nsSNPs with mutation {default: 1}

    name = string: name of new selection {default: nsSNPs}

EXAMPLE

    fetch 3HBT
    snp_uniprot ACTG_HUMAN, chain A

SEE ALSO

    snp_ncbi
    '''
    from Bio import ExPASy
    from Bio import SwissProt
    handle = ExPASy.get_sprot_raw(uniprotname)
    record = SwissProt.read(handle)
    snp_common(record, selection, label, name, quiet, _self=_self)


def snp_ncbi(query, selection='(all)', label=1, name='', quiet=0, *, _self=cmd):
    '''
DESCRIPTION

    Selects all nsSNPs from NCBI (dbSNP) in given structure. Does a sequence
    alignment of reference protein sequence and PDB sequence.

USAGE

    snp_ncbi query [, selection [, label [, name [, quiet]]]]

ARGUMENTS

    query = string: Entrez query, for example a protein accession

    ... see snp_uniprot

EXAMPLE

    fetch 3HBT
    snp_ncbi NP_001092.1[accn], chain A

SEE ALSO

    snp_uniprot
    '''
    from Bio import Entrez, SeqIO
    from Bio.SwissProt import Record
    from lxml import etree

    Entrez.email = emailaddress
    ns = {'ds': 'http://www.ncbi.nlm.nih.gov/SNP/docsum'}
    features = set()

    # get first protein that matches query
    handle = Entrez.esearch(db="protein", term='(%s) AND srcdb_refseq[PROP]' % (query), retmax=1)
    record = Entrez.read(handle)
    if int(record['Count']) == 0:
        print('no such protein')
        return
    id = record['IdList'][0]
    handle = Entrez.efetch(db="protein", id=id, rettype="gb", retmode="text")
    seq = SeqIO.read(handle, 'gb')
    print('Protein match: %s (%s)' % (seq.id, seq.description))
    accn = seq.id
    try:
        protein_acc, protein_ver = accn.split('.')
    except ValueError:
        print('no refseq accession found')
        return

    # get snp-list for protein
    handle = Entrez.esearch(db="snp", term="%s[accn]" % (accn), retmax=100)
    record = Entrez.read(handle)
    print('Number of SNP records: ' + record['Count'])
    if int(record['Count']) > int(record['RetMax']):
        print('Warning: Maximum number of records exceeded (%s)' % (record['RetMax']))
    idlist = ','.join(record['IdList'])

    # xml path to fxnSet nodes
    addr = 'ds:Assembly[@reference="true"]//ds:FxnSet[@protAcc="%s" and @protVer="%s" and @fxnClass="%%s"]' % \
            (protein_acc, protein_ver)

    # download SNPs
    handle = Entrez.efetch(db="snp", id=idlist, rettype="xml", mode="xml")

    document = etree.parse(handle)
    root = document.getroot()
    for rs in root.xpath('/ds:ExchangeSet/ds:Rs', namespaces=ns):
        rsId = rs.get('rsId')
        # reference allele
        try:
            node = rs.xpath(addr % ('reference'), namespaces=ns)[0]
            refRes = node.get('residue')
            refPos = node.get('aaPosition')
            intPos = int(refPos) + 1
        except KeyError:
            print('no ref')
            continue
        # snp alleles
        for node in rs.xpath(addr % ('missense'), namespaces=ns):
            pos = node.get('aaPosition')
            assert pos == refPos
            features.add(('VARIANT', intPos, intPos, '%s -> %s (in dbSNP:rs%s)' %
                    (refRes, node.get('residue'), rsId)))

    # make SwissProt like record
    record = Record()
    record.entry_name = seq.id
    record.sequence = str(seq.seq)
    record.features = sorted(features)
    snp_common(record, selection, label, name, quiet, _self=_self)


cmd.extend('snp_uniprot', snp_uniprot)
cmd.extend('snp_ncbi', snp_ncbi)

# tab-completion of arguments
cmd.auto_arg[1].update({
    'snp_uniprot': cmd.auto_arg[0]['zoom'],
    'snp_ncbi': cmd.auto_arg[0]['zoom'],
})

# vi:expandtab
