#!/usr/bin/env python
#
"""
Sorting datasets into visit-instrument-filter groups.

By Daizhong Liu.

"""

import os, sys, re, datetime, glob, shutil
import numpy as np
from astropy.table import Table
from collections import namedtuple, OrderedDict

# Import JWST pipeline
from jwst import datamodels

# Import click
import click

# Setup logging
import logging

# Import multiprocessing
import multiprocessing

# Define utility functions
def get_script_dir():
    """Get current script file's directory path."""
    return os.path.abspath(os.path.dirname(__file__))

def get_script_name():
    """Get current script file name without the suffix and replaced some characters to underscores."""
    return re.sub(r'[^a-zA-Z0-9_]', r'_', os.path.splitext(os.path.basename(__file__))[0])

def setup_logger():
    logger_streamhandler = logging.StreamHandler()
    logger_streamhandler_formatter = logging.Formatter("[%(asctime)-8s] %(message)s", "%H:%M:%S")
    logger_streamhandler.setFormatter(logger_streamhandler_formatter)
    logger_streamhandler.setLevel(logging.DEBUG)

    log_file = get_script_name()
    log_time = datetime.datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
    logger_filehandler = logging.FileHandler(f"log_{log_file}_{log_time}.txt", mode='a')
    logger_filehandler_formatter = logging.Formatter("[%(asctime)-15s] %(message)s", "%Y-%m-%d %H:%M:%S")
    logger_filehandler.setFormatter(logger_filehandler_formatter)

    logger = logging.getLogger()
    while len(logger.handlers) > 0:
        del logger.handlers[0]
    logger.addHandler(logger_streamhandler)
    logger.addHandler(logger_filehandler)
    logger.setLevel(logging.DEBUG)
    
    return logger

# Parse dataset name (from "go-jwst-parse-dataset-name.py")
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
            raise Exception('Error! The input str {!r} does not seem to have the right format: {!r}'.format(input_str, regex_format))
        return None


def store_into_jwst_dataset_dict(i, d, all_jwst_uncal_files, all_jwst_dataset_name, all_jwst_dataset_info):
    jwst_uncal_file = all_jwst_uncal_files[i]
    with datamodels.open(jwst_uncal_file) as model:
        instrument_name = model.meta.instrument.name
        filter_name = model.meta.instrument.filter
        detector_name = model.meta.instrument.detector
        ra_v1, dec_v1, pa_v3 = model.meta.pointing.ra_v1, model.meta.pointing.dec_v1, model.meta.pointing.pa_v3
        s_region = model.meta.wcsinfo.s_region
    dataset_name = all_jwst_dataset_name[i]
    dataset_info = all_jwst_dataset_info[i]
    d['dataset_name'][i] = dataset_name
    d['proposal_id'][i] = dataset_info.proposal_id
    d['obs_num'][i] = dataset_info.obs_num
    d['visit_num'][i] = dataset_info.visit_num
    d['instrument'][i] = instrument_name
    d['filter'][i] = filter_name
    d['detector'][i] = detector_name
    d['parallel'][i] = dataset_info.parallel
    d['ra_v1'][i] = ra_v1
    d['dec_v1'][i] = dec_v1
    d['pa_v3'][i] = pa_v3
    d['s_region'][i] = s_region



# Main 
@click.command()
@click.argument('jwst_dataset_dir', required=False, default='.', type=click.Path(exists=True))
@click.option('--ncpu', type=int, default=None)
def main(
        jwst_dataset_dir,
        ncpu, 
    ):

    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Prepare output
    output_name = 'table_of_jwst_dataset_visit_instrument_filter_groups.csv'
    
    # check input ncpu
    if ncpu is None:
        ncpu = max(1, os.cpu_count()//2)
    
    # 
    all_jwst_uncal_files = []
    all_jwst_dataset_name = []
    all_jwst_dataset_info = []
    has_error = False
    for dataset_name in os.listdir(jwst_dataset_dir):
        dataset_info = parse_jwst_dataset_name(dataset_name, raise_exception=False)
        if dataset_info is not None:
            uncal_file = os.path.join(jwst_dataset_dir, dataset_name, 'uncals', dataset_name + '_uncal.fits')
            if not os.path.isfile(uncal_file):
                logger.error('Error! File not found: {}'.format(uncal_file))
            all_jwst_uncal_files.append(uncal_file)
            all_jwst_dataset_name.append(dataset_name)
            all_jwst_dataset_info.append(dataset_info)
    if has_error:
        raise Exception('Error! Some files are not found! See previous message.')
    
    # 
    ndataset = len(all_jwst_uncal_files)
    all_jwst_dataset_dict = OrderedDict()
    all_jwst_dataset_dict['dataset_name'] = []
    all_jwst_dataset_dict['proposal_id'] = []
    all_jwst_dataset_dict['obs_num'] = []
    all_jwst_dataset_dict['visit_num'] = []
    all_jwst_dataset_dict['instrument'] = []
    all_jwst_dataset_dict['filter'] = []
    all_jwst_dataset_dict['detector'] = []
    all_jwst_dataset_dict['parallel'] = []
    all_jwst_dataset_dict['ra_v1'] = []
    all_jwst_dataset_dict['dec_v1'] = []
    all_jwst_dataset_dict['pa_v3'] = []
    all_jwst_dataset_dict['s_region'] = []
    
    # for i, jwst_uncal_file in enumerate(all_jwst_uncal_files):
    # 
    #     with datamodels.open(jwst_uncal_file) as model:
    #         instrument_name = model.meta.instrument.name
    #         filter_name = model.meta.instrument.filter
    #         detector_name = model.meta.instrument.detector
    #         ra_v1, dec_v1, pa_v3 = model.meta.pointing.ra_v1, model.meta.pointing.dec_v1, model.meta.pointing.pa_v3
    #         s_region = model.meta.wcsinfo.s_region
    # 
    #     dataset_name = all_jwst_dataset_name[i]
    #     dataset_info = all_jwst_dataset_info[i]
    # 
    #     all_jwst_dataset_dict['dataset_name'].append(dataset_name)
    #     all_jwst_dataset_dict['proposal_id'].append(dataset_info.proposal_id)
    #     all_jwst_dataset_dict['obs_num'].append(dataset_info.obs_num)
    #     all_jwst_dataset_dict['visit_num'].append(dataset_info.visit_num)
    #     all_jwst_dataset_dict['instrument'].append(instrument_name)
    #     all_jwst_dataset_dict['filter'].append(filter_name)
    #     all_jwst_dataset_dict['detector'].append(detector_name)
    #     all_jwst_dataset_dict['parallel'].append(dataset_info.parallel)
    #     all_jwst_dataset_dict['ra_v1'].append(ra_v1)
    #     all_jwst_dataset_dict['dec_v1'].append(dec_v1)
    #     all_jwst_dataset_dict['pa_v3'].append(pa_v3)
    #     all_jwst_dataset_dict['s_region'].append(s_region)
    
    # 
    with multiprocessing.Manager() as manager:
        with multiprocessing.Pool(ncpu) as pool:
            d = manager.dict()
            for key in all_jwst_dataset_dict:
                d[key] = manager.list([None]*ndataset)
            for i in range(ndataset):
                pool.apply_async(store_into_jwst_dataset_dict, args=(i, d, all_jwst_uncal_files, all_jwst_dataset_name, all_jwst_dataset_info))
            pool.close()
            pool.join()
            for key in all_jwst_dataset_dict.keys():
                all_jwst_dataset_dict[key] = list(d[key])
    
    # 
    all_jwst_dataset_table = Table(all_jwst_dataset_dict)
    
    # 
    all_jwst_dataset_table.sort(['proposal_id', 'obs_num', 'visit_num', 'parallel', 'instrument', 'filter', 'detector'])
    
    # 
    # save "list_of_datasets_in_the_same_group.txt"
    output_list_name = 'list_of_jwst_datasets_in_the_same_group.txt'
    grouped = all_jwst_dataset_table.group_by(['proposal_id', 'obs_num', 'visit_num', 'instrument', 'filter'])
    all_jwst_group_id = []
    all_jwst_group_sizes = []
    group_id = 0
    total_group_number = len(grouped.groups)
    for key, group in zip(grouped.groups.keys, grouped.groups):
        group_id += 1
        group_text = '# group_id = {}'.format(group_id) + '\n'
        for irow in range(len(group)):
            group_text += group[irow]['dataset_name'] + '\n'
        for irow in range(len(group)):
            all_jwst_group_id.append(group_id)
            all_jwst_group_sizes.append(len(group))
            dataset_name = group[irow]['dataset_name']
            output_list_file = os.path.join(jwst_dataset_dir, dataset_name, output_list_name)
            if os.path.isfile(output_list_file):
                shutil.move(output_list_file, output_list_file+'.backup')
            with open(output_list_file, 'w') as fp:
                fp.write(group_text)
            logger.info('Output to {!r} (progress {}/{})'.format(output_list_file, group_id, total_group_number))
    
    # 
    all_jwst_dataset_table.add_column(all_jwst_group_id, index=0, name='group_id')
    all_jwst_dataset_table.add_column(all_jwst_group_sizes, index=1, name='group_size')
    
    # 
    output_file = os.path.join(jwst_dataset_dir, output_name)
    if os.path.isfile(output_file):
        shutil.move(output_file, output_file+'.backup')
    all_jwst_dataset_table.write(output_file)
    logger.info('Output to {!r}'.format(output_file))
    
        
        




# Main
if __name__ == '__main__':
    
    main()



