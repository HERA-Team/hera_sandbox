#! /bin/bash

hosts="still4 still5"
for h in $hosts; do
    ssh ${h} "rm /data/zen*" 
done

mysql test --host=10.0.1.20 --user=obs --password=P9ls4R*@ -e "truncate observations; truncate files; truncate orders; truncate history;"
/home/obs/Compress/tests/PopulateTestDB.sh
/home/obs/Compress/tests/qdaemon.test.sh
