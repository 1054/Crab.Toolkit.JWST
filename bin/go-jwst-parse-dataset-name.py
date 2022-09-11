#!/usr/bin/env python
#
"""
Parse JWST data set name, e.g., jwpppppooovvv_ggsaa_eeeee_detector_prodType.

"""

import re
import click
from collections import namedtuple

JWST_Dataset_Name = namedtuple('JWST_Dataset_Name', 
                               ['proposal_id', 'obs_num', 'visit_num', 
                                'visit_group', 'parallel', 'activity', 
                                'exposure', 'detector', 'prod_type'])

def parse_jwst_dataset_name(input_str, raise_exception=True):
    regex_format = r'^jw([0-9]{5})([0-9]{3})([0-9]{3})_([0-9]{2})([0-9]{1})([0-9]{2})_([0-9]{5})_([a-z0-9]+)(.*)$'
    regex_match = re.match(regex_format, input_str)
    if regex_match is not None:
        return JWST_Dataset_Name(*regex_match.groups())
    else:
        if raise_exception:
            raise Exception('Error! The input prefix does not seem to have the right format: {}'.format(regex_format))
        return None

@click.command()
@click.argument('dataset_name', type=str)
def main(
        dataset_name, 
    ):
    
    jwst_dataset = parse_jwst_dataset_name(dataset_name)
    print(jwst_dataset)



if __name__ == '__main__':
    main()



