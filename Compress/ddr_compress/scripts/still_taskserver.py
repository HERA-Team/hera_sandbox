#! /usr/bin/env python
import ddr_compress as ddr
import sys,optparse
import logging; logging.basicConfig(level=logging.DEBUG)

#DATA_DIR = '/data' # where stills put the data they are working on
o = optparse.OptionParser()
o.set_usage('still_taskserver [options] *.uv')
o.set_description(__doc__)
o.add_option('--port',default=14204,type=int,
            help='set port number [default=14204]')
opts, args = o.parse_args(sys.argv[1:])
DATA_DIR = args[0]

dbi = ddr.dbi.DataBaseInterface()
task_server = ddr.task_server.TaskServer(dbi, data_dir=DATA_DIR,port=opts.port)
task_server.start()
