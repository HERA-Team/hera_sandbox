#! /usr/bin/env python
"""
sets exit_status to 1 for all currently open tasks
WARNING: do not run if still execution is in progress unless you really really think its a good idea.
Even if you think its a good idea, it probably isn't.
At the very least, shutdown qdaemon first. -DCJ
"""

from PDB import *
import optparse
import sys
from time import time,sleep
o = optparse.OptionParser()
opts, args = o.parse_args()

pdb.db.query("update history set exit_status=1 where exit_status is NULL;")


