'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException
from pymol.movie import produce_mode_dict, produce_mode_sc

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
        fps=25, twopass=1, vbitrate=16000, quiet=1):
    '''
DESCRIPTION

    Save movie as MPEG1

    This will not be the best quality possible for the file size, but
    we were successfull to play those movies in PowerPoint.

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

    if not quiet:
        print ' save_movie: Rendering frames...'

    tmp_path = tempfile.mkdtemp()
    prefix = os.path.join(tmp_path, 'frame')
    cmd.mpng(prefix, first, last, preserve, mode=mode)

    mencoder = 'mencoder -quiet'
    mpeg1line = '-mf type=png:fps=%d -ovc lavc -forceidx -noskip -of rawvideo' \
            + ' -mpegopts format=mpeg1 -lavcopts vcodec=mpeg1video:vbitrate=%d' \
            + ':vhq:trell:keyint=25'
    mpeg1line = mpeg1line % (fps, vbitrate)
    cmdline = mencoder + ' mf://' + prefix + '* ' + mpeg1line

    if not quiet:
        print ' save_movie: Running mencoder...'

    if twopass:
        if not quiet:
            print ' save_movie: First pass...'
            cmdline1 = cmdline + ':vpass=1'
        subprocess.call(cmdline1.split() + ['-o', os.devnull])
        if not quiet:
            print ' save_movie: Second pass...'
        cmdline = cmdline + ':vpass=2'

    subprocess.call(cmdline.split() + ['-o', filename])

    if not preserve:
        import shutil
        shutil.rmtree(tmp_path)
    elif not quiet:
        print ' save_movie: Not deleting temporary directory:', tmp_path

    if not quiet:
        print ' save_movie: Done'

cmd.extend('frames2states', frames2states)
cmd.extend('save_movie_mpeg1', save_movie_mpeg1)

cmd.auto_arg[0]['frames2states'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['save_movie_mpeg1'] = [produce_mode_sc, 'mode', '']

# vi: ts=4:sw=4:smarttab:expandtab
