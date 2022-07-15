#!/usr/bin/env bash
# 


Crab_BIN_SETUP_SCRIPT=$(dirname "${BASH_SOURCE[0]}")/bin/bin_setup.bash

if [[ $- == *i* ]]; then
    source "$Crab_BIN_SETUP_SCRIPT" -check go-jwst-download-by-proposal-id.py
else
    source "$Crab_BIN_SETUP_SCRIPT"
fi


