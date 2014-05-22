#! /bin/bash

pot=pot0

initDB.py #Create the db

#add hosts
add_host.py qmaster 10.0.1.20 obs
add_host.py still4 10.0.1.24 obs
add_host.py pot0 10.0.1.30 obs
#add_host.py still3 10.0.1.20 obs 
add_host.py still5 10.0.1.30 obs #still5 is doing double duty as pot0


###
# The below is the prototype script for what would happen on paper1
####
#get the list of files to be added
#filelist=`ls -d /home/obs/test_pot/testmicro/zen.* | xargs`
filelist=$*
#Add them to the db.  
echo "adding "`echo ${filelist} | wc -w`" files to pdb from host" `hostname`
new_observation.py ${filelist}

#scp these files over to pot0 and begin compression sequence.
echo moving files to $pot
for f in $filelist; do
    echo ${f}
    infile=qmaster:${f}
    outfile=${pot}:/data/${f##*/}
    md5update.py ${infile}
    record_launch.py -i ${infile} -d '1-RSYNC' ${outfile}
    echo scp -r -c arcfour256 ${infile} ${outfile}
    scp -r -c arcfour256 ${infile} ${outfile}
    if [[ $? ]]; then
        ssh ${pot} "add_file.py ${outfile} -i ${infile}"
        record_completion.py ${outfile}
    else
        record_failure.py ${outfile}
    fi
done
