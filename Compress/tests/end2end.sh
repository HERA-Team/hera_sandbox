#! /bin/bash

TEST=minimumtest

hosts="pot0 still4 still5"
#Clear off the scratch space on the execution hosts
echo "deleting data from /data on hosts: "${hosts} 
for h in $hosts; do
    ssh ${h} "rm -rf /data/zen* 2> /dev/null" 
done

#clear out all the data from the test db
mysql test --host=10.0.1.20 --user=obs --password=P9ls4R*@ -e "truncate observations; truncate files; truncate orders; truncate history;"
echo Compress/tests/PopulateTestDB.sh
/home/obs/Compress/tests/PopulateTestDB.sh /home/obs/Compress/tests/${TEST}/zen*uv
/home/obs/Compress/get_files_to_distill.py 300
echo starting Compress/qdaemon.sh
/home/obs/Compress/qdaemon.sh
