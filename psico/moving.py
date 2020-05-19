'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException
from pymol.movie import produce_mode_dict, produce_mode_sc

def _assert_package_import():
    if not __name__.endswith('.moving'):
        raise CmdException("Must do 'import psico.moving' instead of 'run ...'")

def frames2states(selection, specification):
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
    names = cmd.get_object_list('(' + selection + ')')
    for x in specification.split():
        frame, state = x.split(':')
        cmd.frame(int(frame))
        for name in names:
            cmd.mview('store', object=name, state=int(state))

def save_movie_mpeg1(filename, mode='', first=0, last=0, preserve=0,
        fps=25, twopass=1, vbitrate=16000, quiet=1, exe='mencoder'):
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

    cmd.movie.produce, http://www.freemol.org
    '''
    import os, subprocess, tempfile

    first, last, quiet = int(first), int(last), int(quiet)
    fps, twopass, vbitrate = int(fps), int(twopass), int(vbitrate)

    if cmd.is_string(mode):
        if mode == '':
            if cmd.pymol.invocation.options.no_gui \
                    or cmd.get_setting_boolean('ray_trace_frames'):
                mode = 'ray'
            else:
                mode = 'draw'
        mode = produce_mode_sc.auto_err(mode, 'mode')
        mode = produce_mode_dict[mode]
    mode = int(mode)

    try:
        subprocess.call([exe])
    except OSError:
        print(' Error: Cannot execute "%s"' % (exe))
        raise CmdException

    if not quiet:
        print(' save_movie: Rendering frames...')

    tmp_path = tempfile.mkdtemp()
    prefix = os.path.join(tmp_path, 'frame')
    cmd.mpng(prefix, first, last, preserve, mode=mode)

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

def matrix_to_ttt(names, reverse=0, state=-1, quiet=1):
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
    for object in cmd.get_object_list('(' + names + ')'):
        if ostate < 1:
            state = querying.get_object_state(object)
        matrix = cmd.get_object_matrix(object, state)
        for i in range(cmd.count_states(object)):
            cmd.matrix_reset(object, i+1)
        if reverse:
            cmd.reset(object)
            cmd.transform_object(object, matrix, homogenous=1)
        else:
            cmd.set_object_ttt(object, matrix)

def get_keyframes(quiet=1):
    '''
DESCRIPTION

    Get the list of camera keyframes for the current movie.
    '''
    s = cmd.get_session("none")
    viewelem_list = s["movie"][6]
    if not viewelem_list:
        return []
    r = [i for (i,v) in enumerate(viewelem_list, 1) if v[12] == 2]
    if not int(quiet):
        print(r)
    return r

def get_closest_keyframe():
    '''
    Get the closest movie keyframe, or None if no keyframes defined.
    '''
    keyframes = get_keyframes()
    if not keyframes:
        return None
    current = cmd.get_frame()
    return min(keyframes, key=lambda i: abs(i - current))

def closest_keyframe(quiet=1):
    '''
DESCRIPTION

    Jump to the closest movie keyframe.
    '''
    r = get_closest_keyframe()
    if r is not None:
        cmd.frame(r)
    if not int(quiet):
        print(' Closest Keyframe: ' + str(r))
    return r

def next_keyframe(quiet=1):
    '''
DESCRIPTION

    Jump to the next movie keyframe.
    '''
    keyframes = get_keyframes()
    current = cmd.get_frame()
    keyframes = [i for i in keyframes if i > current]
    if keyframes:
        r = keyframes[0]
        cmd.frame(r)
    else:
        r = None
    if not int(quiet):
        print(' Next Keyframe: ' + str(r))
    return r

def prev_keyframe(quiet=1):
    '''
DESCRIPTION

    Jump to the previous movie keyframe.
    '''
    keyframes = get_keyframes()
    current = cmd.get_frame()
    keyframes = [i for i in keyframes if i < current]
    if keyframes:
        r = keyframes[-1]
        cmd.frame(r)
    else:
        r = None
    if not int(quiet):
        print(' Previous Keyframe: ' + str(r))
    return r

def get_mdo_commands(quiet=1):
    '''
DESCRIPTION

    Get the list of mdo commands.
    '''
    s = cmd.get_session("none")
    commands = s["movie"][5]
    if not int(quiet):
        for frame, command in enumerate(commands, 1):
            if command:
                print('mdo {}: {}'.format(frame, command.strip(';')))
    return commands

def dump_mviews():
    '''
DESCRIPTION

    Dump the current movie as 'set_view' with 'mview store' commands.
    '''
    for frame in get_keyframes() or ():
        cmd.frame(frame)
        print(cmd.get_view(3).strip())
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
