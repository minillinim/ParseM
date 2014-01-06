from distutils.core import setup

setup(
    name='ParseM',
    version='0.0.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['parsem', 'parsem.test'],
    scripts=['bin/ParseM'],
    url='http://pypi.python.org/pypi/ParseM/',
    license='GPLv3',
    description='ParseM',
    long_description=open('README.md').read(),
    install_requires=[],
)

