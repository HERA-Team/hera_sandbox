#! /bin/bash

TEST=micro

hosts="still4 still5"
#Clear off the scratch space on the execution hosts
for h in $hosts; do
    ssh ${h} "rm -rf /data/zen* 2> /dev/null" 
done
rm -r /home/obs/test_pot/test${TEST}/*[DEFz] 2> /dev/null

#clear out all the data from the test db
mysql test --host=10.0.1.20 --user=obs --password=P9ls4R*@ -e "truncate observations; truncate files; truncate orders; truncate history;"
echo Compress/tests/PopulateTestDB.sh
/home/obs/Compress/tests/PopulateTestDB.sh /home/obs/test_pot/testmicro/zen*uv
echo starting Compress/qdaemon.sh
/home/obs/Compress/qdaemon.sh
