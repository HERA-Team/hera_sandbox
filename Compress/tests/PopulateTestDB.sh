#! /bin/bash

pot=still5

initDB.py

#add hosts
add_host.py qmaster 10.0.1.20 obs
add_host.py still4 10.0.1.24 obs
add_host.py pot0 10.0.1.30 obs
#this is a hack and needs to be fixed. 
#hostname and socket.gethostname() don't agree. I can't figure out why...
add_host.py still3 10.0.1.20 obs 
add_host.py still5 10.0.1.30 obs 

#add the test files in pot0:/data
filelist=`ls -d /home/obs/test_pot/zen.*`
filelist=`echo ${filelist}` #this is a hack to have ls print all files to a single line --- no \n
#echo "new_observation.py ${filelist}"
new_observation.py ${filelist}

#rsync these files over to pot0 and begin compression sequence.
for f in $filelist; do
    infile=still3:${f}
    outfile=${pot}:/data/${f##*/}
    md5update.py ${infile}
    record_launch.py -i ${infile} -d '1-RSYNC' ${outfile}
    scp ${infile} ${outfile}
    if [[ $? ]]; then
        ssh ${pot} add_file.py ${outfile} -i ${infile}
        record_completion.py ${outfile}
    else
        echo "DO SOMETHING!"
    fi
done
