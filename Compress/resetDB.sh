#! /bin/bash

function header () {
    echo "#####################################"
    echo "${1}"
    echo "#####################################"
}

header "killing db"
mysql < killDB.sql
header "creating new db"
mysql < createDB.sql
./initDB.py
header "add_host.py"
./add_host.py GrumpyCat 127.0.0.1 damo 
header "new_observation.py"
./new_observation.py GrumpyCat `ls -d test_files/*`
header "md5update.py"
./md5update.py test_files/*
header "record_launch.py"
./record_launch.py test_files/zen.245666.234626.xx.uvR -i test_files/zen.245666.234626.xx.uv -d XRFI
header "add_file.py"
./add_file.py test_files/zen.2455666.234636.xx.uvR -i test_files/zen.2455666.234636.xx.uv
header "record_completion.py"
./record_completion.py test_files/zen.245666.234626.xx.uvR
