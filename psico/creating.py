'''
(c) 2010 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def join_states(name, selection='all', discrete=-1, zoom=0, quiet=1):
    '''
DESCRIPTION

    The reverse of split_states

ARGUMENTS

    name = string: name of object to create or modify
 
    selection = string: atoms to include in the new object

    discrete = -2: match atoms by sequence alignment
    discrete = -1: Assume identical input objects (matching all atom
        identifiers) but also check for missing atoms and only include atoms
        that are present in all input objects {default}
    discrete = 0: Assume identical input objects
    discrete = 1: Input object may be (totally) different
    '''
    discrete, quiet = int(discrete), int(quiet)
    if discrete == -2:
        from .selecting import wait_for
        aln_obj = cmd.get_unused_name('_')
    models = cmd.get_object_list('(' + selection + ')')
    for i in range(len(models)):
        if discrete == -1 and i > 0:
            cmd.remove('(%s) and not (alt A+ and (%s) in (%s))' % (name, name, models[i]))
            cmd.create(name, '(%s) in (%s)' % (models[i], name), 1, i+1, 0, 0, quiet)
        elif discrete == -2 and i > 0:
            cmd.align(models[i], name, cycles=0, transform=0, object=aln_obj)
            wait_for(aln_obj)
            cmd.remove('(%s) and not (%s)' % (name, aln_obj))
            cmd.create(name, name, 1, i+1, 0, 0, quiet)
            cmd.update(name, '(%s) and (%s)' % (models[i], aln_obj), i+1, 1, 0, quiet)
            cmd.delete(aln_obj)
        else:
            cmd.create(name, models[i], 1, i+1, discrete == 1, 0, quiet)
    if int(zoom):
        cmd.zoom(name, state=0)

sidechaincenteratoms = {
    'GLY': ('CA',),
    'ALA': ('CB',),
    'VAL': ('CG1', 'CG2'),
    'ILE': ('CD1',),
    'LEU': ('CD1', 'CD2'),
    'SER': ('OG',),
    'THR': ('OG1', 'CG2'),
    'ASP': ('OD1', 'OD2'),
    'ASN': ('OD1', 'ND2'),
    'GLU': ('OE1', 'OE2'),
    'GLN': ('OE1', 'NE2'),
    'LYS': ('NZ',),
    'ARG': ('NE', 'NH1', 'NH2'),
    'CYS': ('SG',),
    'MET': ('SD',),
    'MSE': ('SE',),
    'PHE': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
    'TYR': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'),
    'TRP': ('CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3'),
    'HIS': ('CG', 'ND1', 'CD2', 'CE1', 'NE2'),
    'PRO': ('CB', 'CG', 'CD'),
}

def sidechaincenters(object='scc', selection='all', name='PS1', method='bahar1996'):
    '''
DESCRIPTION

    Creates an object with sidechain representing pseudoatoms for each residue
    in selection.

    Sidechain interaction centers as defined by Bahar and Jernigan 1996
    http://www.ncbi.nlm.nih.gov/pubmed/9080182

USAGE

    sidechaincenters object [, selection]

ARGUMENTS

    object = string: name of object to create

    selection = string: atoms to consider {default: (all)}

    name = string: atom name of pseudoatoms {default: PS1}

SEE ALSO

    sidechaincentroids, pseudoatom
    '''
    from chempy import Atom, cpv, models

    atmap = dict()
    if method == 'bahar1996':
        modelAll = cmd.get_model('(%s) and resn %s' % (selection, '+'.join(sidechaincenteratoms)))
        for at in modelAll.atom:
            if at.name in sidechaincenteratoms[at.resn]:
                atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)
    elif method == 'centroid':
        modelAll = cmd.get_model('(%s) and polymer and not (hydro or name C+N+O)' % selection)
        for at in modelAll.atom:
            atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)
    else:
        print 'Error: unknown method:', method
        raise CmdException

    model = models.Indexed()
    for centeratoms in atmap.itervalues():
        center = cpv.get_null()
        for at in centeratoms:
            center = cpv.add(center, at.coord)
        center = cpv.scale(center, 1./len(centeratoms))
        atom = Atom()
        atom.coord = center
        atom.index = model.nAtom + 1
        atom.name = name
        for key in ['resn','chain','resi','resi_number','hetatm','ss','segi']:
            atom.__dict__[key] = at.__dict__[key]
        model.add_atom(atom)
    model.update_index()
    if object in cmd.get_object_list():
        cmd.delete(object)
    cmd.load_model(model, object)
    return model

def sidechaincentroids(object='scc', selection='all', name='PS1'):
    '''
DESCRIPTION

    Sidechain centroids. Works like "sidechaincenters", but the
    pseudoatom is the centroid of all atoms except hydrogens and backbone atoms
    (N, C and O).

NOTE

    If you want to exclude C-alpha atoms from sidechains, modify the selection
    like in this example:

    sidechaincentroids newobject, all and (not name CA or resn GLY)

SEE ALSO

    sidechaincenters
    '''
    return sidechaincenters(object, selection, name, method='centroid')

cmd.extend('join_states', join_states)
cmd.extend('sidechaincenters', sidechaincenters)
cmd.extend('sidechaincentroids', sidechaincentroids)

cmd.auto_arg[1].update({
    'join_states'          : cmd.auto_arg[0]['zoom'],
    'sidechaincenters'     : cmd.auto_arg[0]['zoom'],
    'sidechaincentroids'   : cmd.auto_arg[0]['zoom'],
})

# vi: ts=4:sw=4:smarttab:expandtab
