"""Tests for mf_tools.py"""

import nose.tools as nt
import os
import numpy as np
from hera_op.data import DATA_PATH
import hera_op.mf_tools as mt
import ConfigParser as configparser
from configparser import ConfigParser, ExtendedInterpolation

class TestMethods(object):
    def setUp(self):
        self.config_file = os.path.join(DATA_PATH, 'sample_config.cfg')
        self.obsids_pol = ['zen.2458000.12345.xx.uv', 'zen.2458000.12345.xy.uv',
                           'zen.2458000.12345.yx.uv', 'zen.2458000.12345.yy.uv']
        self.obsids_nopol = ['zen.2458000.12345.uv']
        self.pols = ['xx', 'xy', 'yx', 'yy']
        return

    def test_get_config_entry(self):
        # retreive config
        config = ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)

        # retrieve specific entry
        header = 'OMNICAL'
        item = 'prereqs'
        nt.assert_equal(["FIRSTCAL_METRICS"], mt.get_config_entry(config, header, item))

        # get nonexistent, but not required, entry
        header = 'OMNICAL'
        item = 'blah'
        nt.assert_equal([], mt.get_config_entry(config, header, item, required=False))

        # raise an error for a nonexistent, required entry
        nt.assert_raises(AttributeError, mt.get_config_entry, config, header, item)
        return

    def test_make_outfile_name(self):
        # define args
        obsid = self.obsids_pol[0]
        action = 'OMNICAL'
        pols = self.pols
        outfiles = set(['zen.2458000.12345.xx.uv.OMNICAL.xx.out', 'zen.2458000.12345.xx.uv.OMNICAL.xy.out',
                        'zen.2458000.12345.xx.uv.OMNICAL.yx.out', 'zen.2458000.12345.xx.uv.OMNICAL.yy.out'])
        nt.assert_equal(outfiles, set(mt.make_outfile_name(obsid, action, pols)))

        # run for no polarizations
        pols = []
        outfiles = ['zen.2458000.12345.xx.uv.OMNICAL.out']
        nt.assert_equal(outfiles, mt.make_outfile_name(obsid, action, pols))
        return

