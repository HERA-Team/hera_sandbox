import ddr_compress as ddr
import sys,os
from ddr_compress.dbi import Observation
import sys,optparse,os,configparser
import logging; logging.basicConfig(level=logging.DEBUG)
"""
input: list of hosts
action: get all obsids for each host and kill any registered processes
"""

print ddr.__file__
#DATA_DIR = '/data' # where stills put the data they are working on
o = optparse.OptionParser()
o.set_usage('kill_all_tasks.py <host1> <host2> <etc>')
o.set_description(__doc__)
o.add_option('--port',type=int,
            help='set port number [no default]')
o.add_option('--logfile',
            help="optionally send logs to a file instead")
o.add_option('--configfile',
            help='Input a configuration file. see ddr_compress/configs/ for template')
opts, args = o.parse_args(sys.argv[1:])
if not opts.configfile is None:
    configfile = opts.configfile
else:
    configfile = os.path.expanduser('~/.ddr_compress/still.cfg')
#from ddr.task_server import logger
logger = logging.getLogger('taskserver')
logger.setLevel(logging.DEBUG)
if os.path.exists(configfile):
    config = configparser.ConfigParser()
    configfile = os.path.expanduser(configfile)
    if os.path.exists(configfile):
        logger.info('loading file '+configfile)
        config.read(configfile)
    else:
        logging.info(configfile+" Not Found. Exiting")
        sys.exit()
dbi = ddr.dbi.DataBaseInterface()
for host in sys.argv[1:]:
    #make a task client
    taskserverport = int(config[host]['port'])
    taskclient = ddr.task_server.TaskClient(dbi,host,port=taskserverport)
    #get obsnums executing on the indicated host
    s = dbi.Session()
    obsnums = map(int,[obs.obsnum for obs in s.query(Observation).filter(Observation.currentpid!=None,
                                                                Observation.stillhost==host)])
    print "killing processes for obsnums = {obsnums} on {host}".format(obsnums=','.join(map(str,obsnums)),host=host)
    for obsnum in obsnums:
        print "kill!",
        print dbi.get_obs_pid(obsnum)
        taskclient.tx_kill(obsnum)

    s.close()
