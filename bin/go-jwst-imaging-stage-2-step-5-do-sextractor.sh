#/bin/bash
#
# Run SExtractor to make catfile for {dataset}_cal.fits. 
# The output files are in a new "{dataset}_cal_run_sextractor_classic" directory in the same directory as the cal file.
# 
# This requires https://github.com/1054/Crab.Toolkit.SExtractorPlus.
# 

# Set necessary system variables
if [[ -z "$CRDS_PATH" ]]; then
    export CRDS_PATH=$HOME/jwst_crds_cache
fi
if [[ -z "$CRDS_SERVER_URL" ]]; then
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
fi

# Get script_dir
script_dir=$(dirname "${BASH_SOURCE[0]}")

# Check Crab.Toolkit.SExtractorPlus
if [[ -f "$HOME/Cloud/Github/Crab.Toolkit.SExtractorPlus/bin/sextractor_classic_go_find_sources.py" ]]; then
    export PATH="$PATH:$HOME/Cloud/Github/Crab.Toolkit.SExtractorPlus/bin"
fi
if [[ $(type sex 2>/dev/null | wc -l) -eq 0 ]]; then
    echo "Error! \"sex\" is not found in PATH! Please install SExtractor and make sure the command \"sex\" can be found in PATH!"
    exit 255
fi
if [[ $(type sextractor_classic_go_find_sources.py 2>/dev/null | wc -l) -eq 0 ]]; then
    echo "Error! \"sextractor_classic_go_find_sources.py\" is not found in PATH! Please download https://github.com/1054/Crab.Toolkit.SExtractorPlus and add the \"bin\" directory into PATH!"
    exit 255
fi

# Run 
if [[ $# -eq 0 ]]; then
    echo "Please input a JWST imaging cal file, e.g.,"
    echo "    ./this_script jw01837022001_02201_00002_nrca1_cal.fits"
    echo "Options:"
    echo "    see sextractor_classic_go_find_sources.py --help"
    exit 255
fi
if [[ "$1" == *"nrca"* ]] || [[ "$1" == *"nrcb"* ]] || [[ "$1" == *"nrclong"* ]] || [[ "$1" == *"miri"* ]]; then
    echo sextractor_classic_go_find_sources.py --detect-thresh 4.0 $@
    sextractor_classic_go_find_sources.py --detect-thresh 4.0 $@
    if [[ $? -ne 0 ]]; then
        echo "Error?! Please check previous messages."
        exit 255
    fi
else
    echo "Error! The input data is not a NIRCam or MIRI imaging data?"
    exit 255
fi




