#! /bin/bash

mysql < killDB.sql
mysql < createDB.sql
./initDB.py
./add_host.py qmaster 10.0.1.30 obs
