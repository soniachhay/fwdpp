#!/usr/bin/env sh

## This script may help a user install libsequence, 
## but the main purpose is to make sure this dependency
## is met for travis-ci builds

##Install libseq 1.8.5
wget http://github.com/molpopgen/libsequence/archive/1.8.5.tar.gz
tar xzf 1.8.5.tar.gz
cd libsequence-1.8.5
##We don't compile the unit tests here. 
## We'll rely on libseq's own travis-ci setup to handle that in the future
./configure && make && sudo make install