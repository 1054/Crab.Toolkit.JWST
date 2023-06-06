#!/usr/bin/env python
# 

import os, sys, re, json, glob
import numpy as np
import click
from collections import OrderedDict
from astropy.table import Table


@click.command()
@click.argument('input_stats_json', nargs=-1, required=True, type=click.Path(exists=True))
@click.argument('output_csv', nargs=1, required=True, type=click.Path(exists=False))
def main(
        input_stats_json, 
        output_csv, 
    ):


    result_dict = OrderedDict()
    result_dict['dataset'] = []
    result_dict['filter'] = []
    has_all_keys = False
    for json_file in input_stats_json:
        #dataset = re.sub(r'^output_aperture_statistics_(.*)_stats.json$', r'\1', os.path.basename(json_file))
        #print('/n23data2/mfranco/jwst/COSMOS-Web_april_2023/products/pipeline_level2/*/{}.fits'.format(dataset))
        #datafiles = glob.glob('/n23data2/mfranco/jwst/COSMOS-Web_april_2023/products/pipeline_level2/*/{}.fits'.format(dataset))
        #filter_name = os.path.basename(os.path.dirname(datafiles[0]))
        #dataset_filter_name = re.sub(r'^output_aperture_statistics_(.*)_stats.json$', r'\1', os.path.basename(json_file))
        #dataset, filter_name = dataset_filter_name.split('_')
        
        filter_name = os.path.basename(os.path.dirname(json_file))
        
        with open(json_file, 'r') as fp:
            statdict = json.load(fp)
            if isinstance(statdict, list):
                statdicts = statdict
            else:
                statdicts = statdict
            
            for statdict in statdicts:
                if not has_all_keys:
                    for key in statdict:
                        if key not in ['label']:
                            result_dict[key] = []
                    has_all_keys = True
                
                dataset = statdict['image']
                
                result_dict['dataset'].append(dataset)
                result_dict['filter'].append(filter_name)
                for key in result_dict:
                    if key not in ['dataset', 'filter']:
                        if key in result_dict:
                            result_dict[key].append(statdict[key])
                        else:
                            result_dict[key].append(np.nan)

    result_table = Table(result_dict)
    result_table.write(output_csv, format='csv', overwrite=True)
    print('Output to {!r}'.format(output_csv))



if __name__ == '__main__':
    main()


