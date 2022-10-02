#!/usr/bin/env python
# 

import os, sys, re, click
import requests
from astropy.table import Table
from bs4 import BeautifulSoup
from collections import OrderedDict


@click.command()
@click.argument('proposal_id', type=int)
@click.argument('output_table', type=click.Path(exists=False))
def main(
        proposal_id, 
        output_table, 
    ):

    #proposal_id = 1345
    
    url_str = 'https://www.stsci.edu/cgi-bin/get-visit-status?id={:d}&markupFormat=html&observatory=JWST'.format(proposal_id)
    
    #ret = requests.get(url_str)
    
    html_content = requests.get(url_str).text

    soup =  BeautifulSoup(html_content, "lxml" )
    
    #print(soup.prettify())
    
    rows = soup.table.find_all('tr')
    len(rows)

    colhead = [td.text.strip() for td in rows[0].find_all('td')]
    coldict = OrderedDict()
    for i in range(1, len(rows)):
        coldict2 = dict(zip(colhead, [td.text.strip() for td in rows[i].find_all('td')]))
        for colname in colhead:
            if colname not in coldict:
                coldict[colname] = []
            coldict[colname].append(coldict2[colname])

    coltable = Table(coldict)


if __name__ == '__main__':
    main()


