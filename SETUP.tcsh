#!/usr/bin/env tcsh
# 


set BIN_SETUP_SCRIPT = `dirname $0`/bin/bin_setup.bash

set PATH = `bash -c "source $BIN_SETUP_SCRIPT -print" | tail -n 1`

type go-jwst-download-by-proposal-id.py


