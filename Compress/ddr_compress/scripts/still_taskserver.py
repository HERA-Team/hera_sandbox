#! /usr/bin/env python
import ddr_compress as ddr
import sys
import logging; logging.basicConfig(level=logging.DEBUG)

#DATA_DIR = '/data' # where stills put the data they are working on
DATA_DIR = sys.argv[-1]

dbi = ddr.dbi.DataBaseInterface()
task_server = ddr.task_server.TaskServer(dbi, data_dir=DATA_DIR)
task_server.start()
