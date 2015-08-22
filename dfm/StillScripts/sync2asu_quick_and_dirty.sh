#! /bin/bash

asu='damo@enterprise.sese.asu.edu'
finaldir='/data3/PAPER/psa128'
latestdirs=`ssh pot0 "ls -td /data0/psa6* | head -n 4"`
#latestdir=`echo ${latestdir} | cut -f 1 -d " "`

logfile=/home/obs/TransferLogs/ASU-`date "+%F-%H-%M"`.log
for latestdir in $latestdirs
do
    echo "Working on pot0:${latestdir}" >> ${logfile}
    echo >> ${logfile}
    
    files=`ssh pot0 "ls -1d ${latestdir}/zen*E*"` \
        || echo "No data to compress" >> ${logfile}
    
    nfiles=`echo $files | wc -w`
    
    echo "Transferring ${nfiles} files in ${latestdir} to ${asu}:${finaldir}" >> ${logfile}
    
    for f in $files; do
        echo $f >> ${logfile}
        todir=${finaldir}/psa${latestdir##*psa}
        ssh ${asu} test -e ${todir} || ssh ${asu} "mkdir ${finaldir}"
        echo rsync -avP -e "ssh -c arcfour128" ${f} ${asu}:${todir}/ >> ${logfile}
        ssh pot0 "rsync -av -e 'ssh -c arcfour128' ${f} ${asu}:${todir}/" >> ${logfile}
    done
done
