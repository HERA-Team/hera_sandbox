#! /usr/bin/env python
import ddr_compress as ddr
import logging; logging.basicConfic(level=logging.INFO)

DATA_DIR = '/data' # where stills put the data they are working on

dbi = ddr.dbi.DataBaseInterface()
task_server = ddr.task_server.TaskServer(dbi, data_dir=DATA_DIR)
task_server.start()
