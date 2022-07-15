#!/usr/bin/env fish
# 


set BIN_SETUP_SCRIPT (dirname (status --current-filename))/bin/bin_setup.bash

#echo 
#echo "PATH = $PATH"
#echo 

set -x PATH (string split ":" (bash -c "source $BIN_SETUP_SCRIPT -debug -check go-jwst-download-by-proposal-id.py -print" | tail -n 1))

type go-jwst-download-by-proposal-id.py

#echo 
#echo "PATH = $PATH"
#echo 


