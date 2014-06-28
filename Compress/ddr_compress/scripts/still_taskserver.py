#! /usr/bin/env python
import ddr_compress as ddr
import sys,optparse,os,configparser
import logging; logging.basicConfig(level=logging.DEBUG)
print ddr.__file__
#DATA_DIR = '/data' # where stills put the data they are working on
o = optparse.OptionParser()
o.set_usage('still_taskserver [options] *.uv')
o.set_description(__doc__)
o.add_option('--port',default=14204,type=int,
            help='set port number [default=14204]')
o.add_option('--logfile',
            help="optionally send logs to a file instead")
o.add_option('--configfile',
            help='Input a configuration file. see ddr_compress/configs/ for template')
opts, args = o.parse_args(sys.argv[1:])
configfile = os.path.expanduser('~/.ddr_compress/still.cfg')
logger = logging.getLogger('taskserver')
if len(args)==0 and os.path.exists(configfile):#todo config file an inputable thing above
    hostname = ddr.dbi.gethostname()
    config = configparser.ConfigParser()
    configfile = os.path.expanduser(configfile)
    if os.path.exists(configfile):
        logger.info('loading file '+configfile)
        config.read(configfile)
        DATA_DIR = config[hostname]['datadir']
    else:
        logging.info(configfile+" Not Found. Exiting")
        sys.exit()
else:
    DATA_DIR = args[0]
if not opts.logfile is None:
    fh = logging.FileHandler(opts.logfile)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    logger.debug('Starting log file')
dbi = ddr.dbi.DataBaseInterface()
logger.debug('testing db connection')
dbi.test_db()
task_server = ddr.task_server.TaskServer(dbi, data_dir=DATA_DIR,port=opts.port)
task_server.start()
