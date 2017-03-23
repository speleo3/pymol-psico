'''
(c) 2010-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def select_pepseq(pattern, selection='all', name='sele', state=1, quiet=1,
        cutoff=4.0, one_letter=None):
    '''
DESCRIPTION

    Find a amino acid sequence pattern (regular expression) in given atom
    selection. Does not span gaps (unless matched by a wildcard).

USAGE

    select_pepseq pattern [, selection [, name [, state ]]]

ARGUMENTS

    pattern = string: amino acid sequence in one letter code, can be a
    regular expression pattern.

    selection = string: atom selection of protein (non protein atoms in
    selection are silently ignored) {default: all})

    name = a unique name for the selection {default: sele}

EXAMPLE

    fetch 1a00, async=0
    select_pepseq AL[EG]R, chain A+B, sele1
    select_pepseq ([LIVAM]{3,}), all, sele2

SEE ALSO

    There is the pepseq/ps. selection operator. Example:
    select actsite, protein and ps. ADFG

    Similar scripts:
    http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/seq_select.py
    http://pymolwiki.org/index.php/FindSeq
    '''
    import re
    from chempy import cpv

    if not one_letter:
        from . import one_letter

    state, quiet = int(state), int(quiet)
    cutoff = float(cutoff)

    seq_list = []
    idx_list = []
    prev = [1e300, 1e300, 1e300]

    def callback(model, index, resn, coord):
        if cpv.distance(coord, prev) > cutoff:
            seq_list.append('#')
            idx_list.append(None)
        seq_list.append(one_letter.get(resn, '#'))
        idx_list.append((model,index))
        prev[:] = coord

    cmd.iterate_state(state, '(%s) and guide' % (selection),
            'callback(model, index, resn, (x, y, z))', space=locals())

    matches = list(re.finditer(pattern.upper(), ''.join(seq_list)))
    if not quiet:
        if len(matches) == 0:
            print(' select_pepseq: Pattern not found in selection')
        else:
            print(' select_pepseq: Pattern found %d time(s)' % (len(matches)))

    sel_list = []
    for m in matches:
        start, stop = m.span()
        sel_list.extend('%s`%d' % idx for idx in idx_list[start:stop] if idx is not None)

    return cmd.select(name, '(' + selection + ') and byres (none ' + ' '.join(sel_list) + ')')


def select_nucseq(pattern, selection='all', name='sele', state=1, quiet=1):
    '''
DESCRIPTION

    Find a nucleic acid sequence pattern in given atom selection.
    '''
    one_letter = dict(["AA", "CC", "TT", "GG", "UU"])
    return select_pepseq(pattern, selection, name, state, quiet, 6.5, one_letter)


def select_sspick(selection, name=None, caonly=0, quiet=0):
    '''
DESCRIPTION

    Extend selection by connected secondary structure elements.

    Also available as wizard (type: "wizard sspick").

USAGE

    select_sspick selection [, name [, caonly [, quiet ]]]

ARGUMENTS

    selection = string: selection-expression

    name = string: create a named atom selection if not None {default: None}
    '''
    caonly, quiet = int(caonly), int(quiet)

    qkeys = set()
    cmd.iterate('bycalpha (%s)' % (selection),
            'qkeys.add(((model,segi,chain,ss), resv))', space={'qkeys': qkeys})

    def in_intervals(i, intervals):
        for interval in intervals:
            if interval[0] <= i <= interval[1]:
                return True
        return False

    elements = dict()
    for key, resv in qkeys:
        element = elements.setdefault(key, [])
        if in_intervals(resv, element):
            continue
        resv_set = set()
        cmd.iterate('/%s/%s/%s//CA and ss "%s"' % key, 'resv_set.add(resv)',
                space={'resv_set': resv_set})
        resv_min = resv
        resv_max = resv
        while (resv_min - 1) in resv_set:
            resv_min -= 1
        while (resv_max + 1) in resv_set:
            resv_max += 1
        element.append((resv_min, resv_max))

    sele_list = []
    ss_names = {'S': 'Strand', 'H': 'Helix', '': 'Loop', 'L': 'Loop'}
    for key in elements:
        model,segi,chain,ss = key
        for resv_min,resv_max in elements[key]:
            sele = '/%s/%s/%s/%d-%d' % (model, segi, chain, resv_min, resv_max)
            if caonly:
                sele += '/CA'
            sele_list.append(sele)
            if not quiet:
                print("%-6s %s" % (ss_names.get(ss, ss), sele))

    sele = ' '.join(sele_list)
    if name is not None:
        cmd.select(name, sele)
    elif not quiet:
        print(' Selection: ' + sele)

    return sele

def diff(sele1, sele2, byres=1, name=None, operator='in', quiet=0):
    '''
DESCRIPTION

    Difference between two molecules

ARGUMENTS

    sele1 = string: atom selection

    sele2 = string: atom selection

    byres = 0/1: report residues, not atoms (does not affect selection)
    {default: 1}

    operator = in/like/align: operator to match atoms {default: in}

SEE ALSO

    symdiff
    '''
    byres, quiet = int(byres), int(quiet)
    if name is None:
        name = cmd.get_unused_name('diff')
    if operator == 'align':
        alnobj = cmd.get_unused_name('__aln')
        cmd.align(sele1, sele2, cycles=0, transform=0, object=alnobj)
        sele = '(%s) and not %s' % (sele1, alnobj)
        cmd.select(name, sele)
        cmd.delete(alnobj)
    else:
        sele = '(%s) and not ((%s) %s (%s))' % (sele1, sele1, operator, sele2)
        cmd.select(name, sele)
    if not quiet:
        if byres:
            seleiter = 'byca ' + name
            expr = 'print "/%s/%s/%s/%s`%s" % (model,segi,chain,resn,resi)'
        else:
            seleiter = name
            expr = 'print "/%s/%s/%s/%s`%s/%s" % (model,segi,chain,resn,resi,name)'
        cmd.iterate(seleiter, expr)
    return name

def symdiff(sele1, sele2, byres=1, name=None, operator='in', quiet=0):
    '''
DESCRIPTION

    Symmetric difference between two molecules

SEE ALSO

    diff
    '''
    byres, quiet = int(byres), int(quiet)
    if name is None:
        name = cmd.get_unused_name('symdiff')
    tmpname = cmd.get_unused_name('__tmp')
    diff(sele1, sele2, byres, name, operator, quiet)
    diff(sele2, sele1, byres, tmpname, operator, quiet)
    cmd.select(name, tmpname, merge=1)
    cmd.delete(tmpname)
    return name

def collapse_resi(selection='(sele)', quiet=1):
    '''
DESCRIPTION

    Returns a compact selection macro for the given selection on
    residue number (resi) level.

    Rewrite of http://pymolwiki.org/index.php/CollapseSel

ARGUMENTS

    selection = string: atom selection {default: (sele)}
    '''
    from collections import defaultdict
    s_dict = defaultdict(set)
    cmd.iterate(selection, 's_dict[model,segi,chain].add(resv)', space=locals())
    r_all = []
    for key, s in s_dict.items():
        s = sorted(s)
        r = [[s[0], s[0]]]
        for i in s[1:]:
            if i <= r[-1][1] + 1:
                r[-1][1] = i
            else:
                r.append([i,i])
        resi = '+'.join(('%d-%d' % (f,t) if f != t else '%d' % (f)) for (f,t) in r)
        r_all.append('/%s/%s/%s/' % key + resi)
    if not int(quiet):
        for r in r_all:
            print(' collapse_resi: ' + str(r))
    return ' '.join(r_all)

def wait_for(name, state=0, quiet=1):
    '''
DESCRIPTION

    Wait for "name" to be available as selectable object.
    '''
    if cmd.count_atoms('?' + name, 1, state) == 0:
        s = cmd.get_setting_boolean('suspend_updates')
        if s: cmd.set('suspend_updates', 0)
        cmd.refresh()
        if s: cmd.set('suspend_updates')

def select_distances(names='', name='sele', state=1, selection='all', cutoff=-1, quiet=1):
    '''
DESCRIPTION

    Turns a distance object into a named atom selection.

ARGUMENTS

    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}

    name = a unique name for the selection {default: sele}

    state = int: object state (-1: current, 0: all states) {default: 1}

SEE ALSO

    get_raw_distances
    '''
    from collections import defaultdict
    from .querying import get_raw_distances

    state, cutoff, quiet = int(state), float(cutoff), int(quiet)
    states = [state] if state else list(range(1, cmd.count_states(selection)+1))

    sele_dict = defaultdict(set)
    for state in states:
        distances = get_raw_distances(names, state, selection)
        for idx1, idx2, dist in distances:
            if cutoff <= 0.0 or dist <= cutoff:
                sele_dict[idx1[0]].add(idx1[1])
                sele_dict[idx2[0]].add(idx2[1])

    cmd.select(name, 'none')
    tmp_name = cmd.get_unused_name('_')

    r = 0
    for model in sele_dict:
        cmd.select_list(tmp_name, model, list(sele_dict[model]), mode='index')
        r = cmd.select(name, tmp_name, merge=1)
        cmd.delete(tmp_name)

    if not quiet:
        print(' Selector: selection "%s" defined with %d atoms.' % (name, r))
    return r

# commands
cmd.extend('select_pepseq', select_pepseq)
cmd.extend('select_nucseq', select_nucseq)
cmd.extend('select_sspick', select_sspick)
cmd.extend('symdiff', symdiff)
cmd.extend('diff', diff)
cmd.extend('collapse_resi', collapse_resi)
cmd.extend('select_distances', select_distances)

# autocompletion
cmd.auto_arg[0].update([
    ('select_sspick', cmd.auto_arg[0]['align']),
    ('symdiff',     cmd.auto_arg[0]['align']),
    ('diff',        cmd.auto_arg[0]['align']),
    ('collapse_resi', cmd.auto_arg[0]['zoom']),
    ('select_distances', [
        lambda: cmd.Shortcut(cmd.get_names_of_type('object:measurement')),
        'distance object', '']),
])
cmd.auto_arg[1].update([
    ('select_pepseq', cmd.auto_arg[1]['select']),
    ('select_nucseq', cmd.auto_arg[1]['select']),
    ('symdiff',     cmd.auto_arg[1]['align']),
    ('diff',        cmd.auto_arg[1]['align']),
])

# vi:expandtab:smarttab
