#!  /usr/bin/env python
import optparse,sys,os
from ddr_compress.dbi import DataBaseInterface

#open a connection to the db
dbi = DataBaseInterface()
#check that the tables 
