#!/usr/bin/env python
#
#  pull_facets.py
#  
#
#  Created by Danny Jacobs on 9/3/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#

import optparse,re,sys

o = optparse.OptionParser()
o.set_usage('pull_facets.py <files>')
o.set_description(__doc__)
o.add_option('--match',type='str',default='z.*_(\d+).*fits',
     help="Provide regular expression for extracting the facet number. [z.*_(\d+).*fits]")
opts,args = o.parse_args(sys.argv[1:])


prog = re.compile(opts.match)
facets = []
for filename in args:
    result = prog.search(filename)
    if result is None:continue
    facets.append(result.groups()[0])
print ' '.join(list(set(facets)))