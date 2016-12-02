#!/bin/bash
PW=`python -c "from ddr_compress.dbi import dbinfo; print dbinfo['password']"`
echo ${PW}
DB=paperdistiller
# Where is dbinfo is getting its info?  Should be getting from still.cfg, but it's not
#`python -c "from ddr_compress.dbi import dbinfo; print dbinfo['dbname']"`
echo ${DB}
mysql ${DB} -h shredder --user=obs --password=${PW} -e "set foreign_key_checks=0; truncate  neighbors; truncate observation; truncate file; truncate log; set foreign_key_checks=1;"
python -c "from ddr_compress.dbi import DataBaseInterface; dbi = DataBaseInterface(); dbi.test_db()"