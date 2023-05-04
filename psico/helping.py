'''
(c) 2011 Thomas Holder, MPI for Developmental Biology
(c) 2005-2009 Ezequiel Panepucci

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def grepset(regexp='', *, _self=cmd):
    '''
DESCRIPTION

    "grepset" greps through the list of settings using a python
    regular expression as defined in the 're' module.
    It returns a list of settings/values matching the regexp.
    No regexp returns every setting.

    Author: Ezequiel Panepucci
    http://pymolwiki.org/index.php/Grepset

USAGE

    grepset [regexp]

EXAMPLE

    grepset line
    grepset ray
    grepset (^line|color$)

SEE ALSO

    apropos
    '''
    import re
    import pymol.setting
    count = 0
    regexp = re.compile(regexp)
    matches = []
    for a in pymol.setting.get_index_list():
        setting = pymol.setting._get_name(a)
        if regexp.search(setting):
            count += 1
            matches.append((setting, _self.get_setting_text(a, '', -1)))
    # max length of the setting names that matched
    maxlen = max([len(s[0]) for s in matches] + [0])
    fmt = "%%-%ds : %%s" % (maxlen,)
    for setting in matches:
        print(fmt % setting)
    print('%d settings matched' % (count,))


def apropos(regexp=''):
    '''
DESCRIPTION

    "apropos" searches through the documentation of all currently
    defined commands and lists those commands for which the keyword
    is either contained in the documentation or matches the command
    name itself.

    If an appropriate "DESCRIPTION" section is provided in the documentation
    of the command, the first 80 characters are listed as a summary.

    Author: Ezequiel Panepucci
    http://pymolwiki.org/index.php/Apropos

USAGE

    apropos [keyword or regexp]

SEE ALSO

    grepset
    '''
    import re

    count = 0
    matches_with_help = []
    matches_without_help = []

    maxcclen = 0
    for cc in cmd.keyword:
        if cc == regexp:
            print('\n###EXACT MATCH FOR: %s ==> try \'help %s\' at the prompt.' % (cc, cc))

        doc = cmd.keyword[cc][0].__doc__

        if doc is None:
            if re.search(regexp, cc, re.IGNORECASE):
                count += 1
                matches_without_help.append(cc)
            continue

        if re.search(regexp, doc, re.MULTILINE | re.IGNORECASE):
            count += 1
            if len(cc) > maxcclen:
                maxcclen = len(cc)

            docmatches = re.match(r"""^\s+DESCRIPTION\s+(.{0,80})\S*""", doc, re.IGNORECASE)
            if docmatches is None:
                desc = '>>>>>>>> Ooopsie, no DESCRIPTION found for this command!!! <<<<<<'
            else:
                desc = docmatches.group(1)
            matches_with_help.append((cc, desc))

    if len(matches_without_help) > 0:
        print('\n###The following commands are NOT documented.\n')
        for cc in matches_without_help:
            print('%*s' % (maxcclen, cc))

    if len(matches_with_help) > 0:
        print('\n###The following commands are documented.  \'help command\' \n')
        for cc, desc in matches_with_help:
            print('%*s : %s' % (maxcclen, cc, desc))


def api_info(name):
    '''
DESCRIPTION

    Get the full function name (incl. module) of given command.

ARGUMENTS

    name = string: name of a PyMOL command
    '''
    import sys
    name = cmd.kwhash.shortcut.get(name, name)
    try:
        func = cmd.keyword[name][0]
    except KeyError:
        raise CmdException('No such command')
    print(' CMD: ' + str(name))
    print(' API: %s.%s' % (func.__module__, func.__name__))
    if func == getattr(cmd, func.__name__, None):
        print(' API: cmd.' + func.__name__)
    print(' FILE: ' + str(sys.modules[func.__module__].__file__))
    return func


def write_html_ref(filename, prefix='psico', format='html'):
    '''
DESCRIPTION

    Write psico command reference to file.

SEE ALSO

    cmd.write_html_ref
    '''
    if prefix == 'psico' and not __name__.startswith(prefix):
        prefix = __name__.rsplit('.', 1)[0]
    ref = []
    for a in sorted(cmd.keyword):
        if a.startswith('_'):
            continue
        func = cmd.keyword[a][0]
        if (func.__module__, func.__name__) == ('pymol.helping', 'python_help'):
            continue
        if func.__module__.startswith(prefix) and hasattr(func, '__doc__'):
            if isinstance(func.__doc__, str):
                ref.append([a, func])
    f = open(filename, 'w')

    if format == 'txt_short':
        import re
        pattern_header = re.compile(r'^[A-Z]+', flags=re.MULTILINE)
        pattern_emptyline = re.compile(r'\s+$', flags=re.MULTILINE)
        for (a, func) in ref:
            doc = func.__doc__.partition('DESCRIPTION')[-1]
            doc = pattern_header.split(doc, 1)[0]
            doc = pattern_emptyline.sub('', doc)
            f.write('%s\t%s\n\n' % (a, doc))
        f.close
        return

    f.write("<html><head><style type='text/css'>")
    f.write("body {font-family:sans-serif}")
    f.write("p.api {font:small monospace;color:#999}")
    f.write("</style></head><body><h1>pymol psico reference</h1><ul>")
    for (a, _) in ref:
        f.write("<li><a href='#%s'>%s</a>" % (a, a))
    f.write("</ul>")
    for (a, func) in ref:
        doc = func.__doc__.strip().replace("<", "&lt;")
        f.write("<hr size=1><h2 id='%s'>%s</h2>" % (a, a))
        f.write("<pre>%s</pre>" % (doc))
        f.write("<p class='api'>api: %s.%s</p>" % (func.__module__, func.__name__))
    f.write("</body></html>")
    f.close()


def write_txt_ref(filename, prefix='psico'):
    '''
DESCRIPTION

    Write psico command and its DESCRIPTION to file as plain text.

SEE ALSO

    cmd.write_html_ref
    '''
    write_html_ref(filename, prefix, 'txt_short')


cmd.extend('grepset', grepset)
cmd.extend('apropos', apropos)
cmd.extend('api_info', api_info)

cmd.auto_arg[0]['api_info'] = [cmd.kwhash, 'command', '']

# vi:expandtab:smarttab
