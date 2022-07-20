#!/usr/bin/env python
#
"""
Print CRDS default context.

By Daizhong Liu.

"""

import os
#assert os.environ["CRDS_PATH"] != ''
#assert os.environ["CRDS_SERVER_URL"] != ''
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# Import CRDS
import crds

print(crds.get_default_context())


