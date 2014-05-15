#ParseM

## Overview

ParseM is an attempt to put all the little bits of Python-based biological file parsing snippets into one place. I am aware that other groups are doing this (in a much nicer fashion than I) so I'm not going to try make the best Fasta parser ever, although there is one here. The main object of ParseM is to store a BAM-file parser I'm writing in C and wrapping in Python. I'm using htslib to do this. This part of the code is intended to provide a faster, more stable interface to parsing BAM files than PySam.

## Installation

Dependencies:

The BAM parsing is done using c and a few external libraries. This slightly complicates the makes ParseM installation process but fortunately not too much.

The notes here are for installing on Ubuntu, but should be transferrable to any other Linux system. I'm not sure if these notes are transferabble to fashionable overpriced systems with rounded rectangles and retina displays and I'm almost cetain you'll need some sysadmin-fu to get things going on a Windows system. If you do ge it all set up then please let me know how and I'll buy you a 5 shot venti, 2/5th decaf, ristretto shot, 1pump Vanilla, 1pump Hazelnut, breve,1 sugar in the raw, with whip, carmel drizzle on top, free poured, 4 pump mocha.

First, you need pip, git, zlib and a C-compiler. On Ubuntu this looks like:

    sudo apt-get -y install git build-essential python-pip zlib1g-dev

Next you'll need htslib (Samtools guts) and libcfu (hash objects).

###Notes on installing htslib:
Get the latest htslib from github:

    git clone https://github.com/samtools/htslib.git

For various resons we need to install a statically linked version of htslib. When making use this command instead of just 'make':

    make CFLAGS='-g -Wall -O2 -fPIC -static-libgcc -shared'

###Notes on installing libcfu:
I have built this librray around libcfu 0.03. It is available here:

    http://downloads.sourceforge.net/project/libcfu/libcfu/libcfu-0.03/libcfu-0.03.tar.gz

On my system I have trouble installing it becuase of inconsistencies with print statements. I fix this with sed like so and then install locally:

    tar -xvf libcfu-0.03.tar.gz
    cd libcfu-0.03
    sed -i -e "s/%d/%zd/g" examples/*.c
    sed -i -e "s/%u/%zu/g" examples/*.c
    # Remove the '--prefix=' part to install system wide
    ./configure --prefix=`pwd`/build
    # Code in the examples folder breaks compilation so nuke this from the Makefile
    sed -i -e "s/src examples doc/src doc/" Makefile
    # make and install
    make CFLAGS='-g -Wno-unused-variable -O2 -fPIC -static-libgcc -shared'
    make install

If you install these libraries to local folders (e.g. in your home folder) then you need to take note of where you installed them. If you installed them system-wide then it *should* be no hassle.

###Install ParseM
Get the latest version from github:

    git clone https://github.com/minillinim/ParseM.git

If you installed htslib and libcfu system-wide then installation is very straight forward:

    python setup.py install

If you installed one or more of these libraries locally then you need to tell setup.py where they are:

    python setup.py install --with-libhts-lib /path/to/htslib --with-libhts-inc /path/to/htslib --with-libcfu-inc /path/to/libcfu/include/ --with-libcfu-lib path/to/libcfu/lib/

Relative paths are OK. You can add the --prefix flag to setup.py to install ParseM locally. Once done, don't forget to add ParseM to your PYTHONPATH.

## Example usage

## Help

If you experience any problems using ParseM, open an [issue](https://github.com/minillinim/ParseM/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/ParseM

This software is currently unpublished

## Copyright

Copyright (c) 2014 Michael Imelfort, Donovan Parks. See LICENSE.txt for further details.
