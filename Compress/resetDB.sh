#! /bin/bash

mysql < killDB.sql
mysql < createDB.sql
./initDB.py
./add_host.py qmaster 10.0.1.30 obs
./new_observation.py qmaster zen.245666.456456.xx.uv zen.245666.23452.xx.uv zen.245666.234626.xx.uv
