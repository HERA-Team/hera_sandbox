#! /usr/bin/env python

from PDB import *
import sys

pdb.delrow('orders', sys.argv[1])
