# ParseM

## Overview

## Installation

Dependencies:

The BAM parsing is done using c and a few external libraries. At minimum you will need:

    gcc
    libcfu
    htslib

Notes on installing htslib:
For various resons we need to install a statically linked version of htslib. When making use this command instead of just 'make':

make CFLAGS='-g -Wall -O2 -fPIC -static-libgcc -shared'

Then it should be as simple as

    pip install ParseM

## Example usage

## Help

If you experience any problems using ParseM, open an [issue](https://github.com/minillinim/ParseM/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/ParseM

This software is currently unpublished

## Copyright

Copyright (c) 2014 Michael Imelfort, Donovan Parks. See LICENSE.txt for further details.
