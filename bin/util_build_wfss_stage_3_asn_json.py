#!/usr/bin/env python
#
# Build asn.json for input WFSS stage 3 cal files
# 
# see 'util_extract_wfss_source_from_cal_files.py' instead
# 
import os, sys, re, json, copy, glob, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from pprint import pprint
from tqdm import tqdm
from collections import OrderedDict

template_dict = OrderedDict([
    ("asn_type", "None"),
    ("asn_rule", "DMS_Level3_Base"),
    ("version_id", None),
    ("code_version", "1.20.2"),
    ("degraded_status", "No known degraded exposures in association."),
    ("program", "00000"),
    ("constraints", "No constraints"),
    ("asn_id", "001"),
    ("asn_pool", "none"),
    ("products", [{"name": "", "members": []},])
])


@click.command()
@click.argument('input_cal_files', type=str, nargs=-1, required=True)
@click.option('--output-name', '--output', 'output_name', type=str, default=None)
def main(input_cal_files, output_name):
    asn_dict = copy.deepcopy(template_dict)
    program = ""
    for cal_file in input_cal_files:
        asn_dict["products"][0]["members"].append(
            {"expname": cal_file, "exptype": "science"}
        )
        cal_filename = os.path.basename(cal_file)
        if program == "":
            regex_match = re.match(r'^(.*_|)jw([0-9]{5}).*', cal_filename)
            if regex_match:
                program = regex_match.group(2)
    
    if program != "":
        asn_dict["program"] = program

    if output_name is None:
        output_name = "output"
    asn_dict["products"][0]["name"] = output_name

    if os.path.exists("asn.json"):
        shutil.move("asn.json", "asn.json.backup")
    with open("asn.json", "w") as fp:
        json.dump(asn_dict, fp, indent=4)
    print("Output to {!r}".format("asn.json"))

    if os.path.exists("run_calwebb_spec3.py"):
        shutil.move("run_calwebb_spec3.py", "run_calwebb_spec3.py.backup")
    with open("run_calwebb_spec3.py", "w") as fp:
        fp.write("#!/usr/bin/env python\n")
        fp.write("#\n")
        fp.write("import os, sys, re, json, copy, datetime, time, glob, shutil\n")
        fp.write("from jwst import datamodels\n")
        fp.write("from jwst.pipeline import calwebb_spec3\n")
        fp.write("from jwst.pipeline import calwebb_image3\n")
        fp.write("from jwst.associations import asn_from_list\n")
        fp.write("from jwst.associations.lib.rules_level3_base import DMS_Level3_Base\n")
        fp.write("\n")
        fp.write("pipeline_object = calwebb_spec3.Spec3Pipeline()\n")
        fp.write("pipeline_object.output_file = \"{}\"\n".format(output_name))
        fp.write("pipeline_object.save_results = True\n")
        fp.write("pipeline_object.run(\"{}\")\n".format("asn.json"))
        fp.write("\n")
    print("Output to {!r}".format("run_calwebb_spec3.py"))



if __name__ == '__main__':
    main()



