#! /usr/bin/env python
import pstats, sys

s = pstats.Stats(sys.argv[-1])
s.sort_stats('cumulative').print_stats()

