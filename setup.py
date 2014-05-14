from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
import sys
from subprocess import call
from os.path import join, abspath
from os import chdir, getcwd

xtra_opts = {"--with-libcfu-inc":"libcfu headers at this location",
             "--with-libcfu-lib":"libcfu library at this location",
             "--with-libhts-inc":"htslib headers at this location",
             "--with-libhts-lib":"htslib library at this location"}


if '--help' not in sys.argv:
    # set the location of the compiled c library to the place where the headers will be
    for scheme in INSTALL_SCHEMES.values():
        scheme['data'] = join(scheme['platlib'], 'parsem', 'c', 'bam')

    # grab extra configuration arguments
    configure_args = []
    for opt in xtra_opts.keys():
        try:
            opt_idx = sys.argv.index(opt)
            configure_args.append(opt+"="+abspath(sys.argv[opt_idx+1]))
            del sys.argv[opt_idx+1]
            sys.argv.remove(opt)
        except ValueError:
            pass

    # configure and make the c portion of the program
    cur_dir = getcwd()
    chdir(join(cur_dir, 'c', 'bam'))
    call([join(getcwd(), "configure")] + configure_args)
    call(['make','clean'])
    call(['make'])
    chdir(cur_dir)
else:
    print
    print "Embedded C options (for building libPMBam.a) USE: --OPTION<space>PATH"
    for opt in xtra_opts.keys():
        print "  %s  %s"%(opt,xtra_opts[opt])
    print

# return to regular viewing
setup(
    name='ParseM',
    version='0.0.3',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['parsem', 'parsem.test'],
    url='http://pypi.python.org/pypi/ParseM/',
    license='GPLv3',
    description='ParseM',
    long_description=open('README.md').read(),
    install_requires=[],
    data_files=[('', ['c/bam/libPMBam.a'])]
)

