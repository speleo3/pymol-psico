'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException
from pymol.movie import produce_mode_dict, produce_mode_sc

def _assert_package_import():
    if not __name__.endswith('.moving'):
        raise CmdException("Must do 'import psico.moving' instead of 'run ...'")

def frames2states(selection, specification, *, _self=cmd):
    '''
DESCRIPTION

    Map specific object states to frames.
    "specification" is a sequence of "frame:state" mappings.

EXAMPLE

    fetch 1d7q 1nmr, async=0
    mset 1-20
    # reverse state order just for 1nmr
    frames2states 1nmr, 1:20 20:1
    '''
    names = _self.get_object_list('(' + selection + ')')
    for x in specification.split():
        frame, state = x.split(':')
        _self.frame(int(frame))
        for name in names:
            _self.mview('store', object=name, state=int(state))

def save_movie_mpeg1(filename, mode='', first=0, last=0, preserve=0,
        fps=25, twopass=1, vbitrate=16000, quiet=1, exe='mencoder',
        *, _self=cmd):
    '''
DESCRIPTION

    Save movie as MPEG1

    This will not be the best quality possible for the file size, but
    we were successfull to play those movies in PowerPoint.

    Requires mencoder executable from http://mplayerhq.hu

ARGUMENTS

    filename = string: Filename, should end on '.mpeg'

    mode = normal | draw | ray

    first, last = integer: first and last frame number

    preserve = 0 or 1: delete temporary images if 0 {default: 0}

    fps = integer: frames per second {default: 25}

    twopass = 0 or 1: two pass mode encoding (improves quality) {default: 1}

    vbitrate = integer: average video bitrate {default: 16000}
    WARNING: 4-16000 (in kbit), 16001-24000000 (in bit)

SEE ALSO

    cmd.movie.produce
    '''
    import os, subprocess, tempfile

    first, last, quiet = int(first), int(last), int(quiet)
    fps, twopass, vbitrate = int(fps), int(twopass), int(vbitrate)

    if isinstance(mode, str):
        if mode == '':
            if _self.pymol.invocation.options.no_gui \
                    or _self.get_setting_boolean('ray_trace_frames'):
                mode = 'ray'
            else:
                mode = 'draw'
        mode = produce_mode_sc.auto_err(mode, 'mode')
        mode = produce_mode_dict[mode]
    mode = int(mode)

    try:
        subprocess.call([exe])
    except OSError:
        raise CmdException('Cannot execute "%s"' % (exe))

    if not quiet:
        print(' save_movie: Rendering frames...')

    tmp_path = tempfile.mkdtemp()
    prefix = os.path.join(tmp_path, 'frame')
    _self.mpng(prefix, first, last, preserve, mode=mode)

    mpeg1line = '-mf type=png:fps=%d -ovc lavc -forceidx -noskip -of rawvideo' \
            + ' -mpegopts format=mpeg1 -lavcopts vcodec=mpeg1video:vbitrate=%d' \
            + ':vhq:trell:keyint=25'
    mpeg1line = mpeg1line % (fps, vbitrate)
    cmdline = exe + ' -quiet mf://' + prefix + '* ' + mpeg1line

    if not quiet:
        print(' save_movie: Running mencoder...')

    if twopass:
        if not quiet:
            print(' save_movie: First pass...')
            cmdline1 = cmdline + ':vpass=1'
        subprocess.call(cmdline1.split() + ['-o', os.devnull])
        if not quiet:
            print(' save_movie: Second pass...')
        cmdline = cmdline + ':vpass=2'

    subprocess.call(cmdline.split() + ['-o', filename])

    if not preserve:
        import shutil
        shutil.rmtree(tmp_path)
    elif not quiet:
        print(' save_movie: Not deleting temporary directory: ' + tmp_path)

    if not quiet:
        print(' save_movie: Done')

def matrix_to_ttt(names, reverse=0, state=-1, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Objects can have state matrices and view (frames) matrices. This function
    takes the total matrix and stores it either as view matrix or as state
    matrix (reverse=1). For movie frames, movie_auto_store must be set.
    '''
    _assert_package_import()
    from . import querying
    reverse, state, quiet = int(reverse), int(state), int(quiet)
    ostate = state
    for object in _self.get_object_list('(' + names + ')'):
        if ostate < 1:
            state = querying.get_object_state(object, _self=_self)
        matrix = _self.get_object_matrix(object, state)
        for i in range(_self.count_states(object)):
            _self.matrix_reset(object, i+1)
        if reverse:
            _self.reset(object)
            _self.transform_object(object, matrix, homogenous=1)
        else:
            _self.set_object_ttt(object, matrix)

def get_keyframes(quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Get the list of camera keyframes for the current movie.
    '''
    s = _self.get_session("none")
    viewelem_list = s["movie"][6]
    if not viewelem_list:
        return []
    r = [i for (i,v) in enumerate(viewelem_list, 1) if v[12] == 2]
    if not int(quiet):
        print(r)
    return r

def get_closest_keyframe(*, _self=cmd):
    '''
    Get the closest movie keyframe, or None if no keyframes defined.
    '''
    keyframes = get_keyframes(_self=_self)
    if not keyframes:
        return None
    current = _self.get_frame()
    return min(keyframes, key=lambda i: abs(i - current))

def closest_keyframe(quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Jump to the closest movie keyframe.
    '''
    r = get_closest_keyframe(_self=_self)
    if r is not None:
        _self.frame(r)
    if not int(quiet):
        print(' Closest Keyframe: ' + str(r))
    return r

def next_keyframe(quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Jump to the next movie keyframe.
    '''
    keyframes = get_keyframes(_self=_self)
    current = _self.get_frame()
    keyframes = [i for i in keyframes if i > current]
    if keyframes:
        r = keyframes[0]
        _self.frame(r)
    else:
        r = None
    if not int(quiet):
        print(' Next Keyframe: ' + str(r))
    return r

def prev_keyframe(quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Jump to the previous movie keyframe.
    '''
    keyframes = get_keyframes(_self=_self)
    current = _self.get_frame()
    keyframes = [i for i in keyframes if i < current]
    if keyframes:
        r = keyframes[-1]
        _self.frame(r)
    else:
        r = None
    if not int(quiet):
        print(' Previous Keyframe: ' + str(r))
    return r

def get_mdo_commands(quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Get the list of mdo commands.
    '''
    s = _self.get_session("none")
    commands = s["movie"][5] or []
    if not int(quiet):
        for frame, command in enumerate(commands, 1):
            if command:
                print('mdo {}: {}'.format(frame, command.strip(';')))
    return commands

def dump_mviews(*, _self=cmd):
    '''
DESCRIPTION

    Dump the current movie as 'set_view' with 'mview store' commands.
    '''
    for frame in get_keyframes(_self=_self) or ():
        _self.frame(frame)
        print(_self.get_view(3).strip())
        print('mview store, {}'.format(frame))

cmd.extend('frames2states', frames2states)
cmd.extend('save_movie_mpeg1', save_movie_mpeg1)
cmd.extend('matrix_to_ttt', matrix_to_ttt)
cmd.extend("get_keyframes", get_keyframes)
cmd.extend("closest_keyframe", closest_keyframe)
cmd.extend("next_keyframe", next_keyframe)
cmd.extend("prev_keyframe", prev_keyframe)
cmd.extend("get_mdo_commands", get_mdo_commands)
cmd.extend("dump_mviews", dump_mviews)

cmd.auto_arg[0]['matrix_to_ttt'] = cmd.auto_arg[0]['disable']
cmd.auto_arg[0]['frames2states'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['save_movie_mpeg1'] = [produce_mode_sc, 'mode', '']

# vi: ts=4:sw=4:smarttab:expandtab
