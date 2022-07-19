#!/usr/bin/env python
# 

import os, sys, re, click
from urllib.parse import quote, unquote


@click.command()
@click.argument('proposal_id', type=int)
def main(proposal_id):

    #proposal_id = 1345
    
    url_base = 'https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery='

    url_str = """
    {{
        "service": "CAOMFILTERED",
        "inputText":
        [
            {{
                "paramName": "proposal_id",
                "niceName": "proposal_id",
                "values":
                [],
                "valString": "{0}",
                "isDate": false,
                "freeText": "{0}",
                "displayString": "{0}"
            }}
        ],
        "position": "undefined, undefined, undefined",
        "paramsService": "Mast.Caom.Filtered",
        "title": "MAST:  Advanced Search 1",
        "tooltip": "{0}; ",
        "columns": "*"
    }}
    """.format(
        proposal_id
        )

    url_str_clean = re.sub(r'[\n ]', r'', url_str)
    print(url_base + url_str_clean)
    print(url_base + quote(url_str_clean))



if __name__ == '__main__':
    main()


