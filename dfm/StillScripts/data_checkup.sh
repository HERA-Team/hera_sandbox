#! /bin/bash

paperblog=teampaper.stats@blogger.com
#mailto="dfmuchicago@gmail.com wheresmytab@gmail.com jeaupenn@gmail.com"
#mailto="paper@lists.berkeley.edu"
mailto="paper@mail.astro.virginia.edu"
today=`date "+%F"`
report=""

newline=$'\n'
function add_to() {
    report=${report}${newline}${1}
}

error="!!!BAD!!!"
nerr=0

#Is there new data?
add_to "--- Transfer summary from paper1 ---"

for oe in {0..0}
do
    list_of_data=`ssh pot${oe} "ls -d /data${oe}/zen.*uv"`
    if [ $? -eq 0 ]
    then
        list_of_data=$(( `echo $list_of_data | wc -w` ))
        add_to "${list_of_data} files have been successfully transfered to pot${oe} today."
        if [ $list_of_data == 288 ]
        then
            add_to " - this is correct"
        else
            add_to $error; (( nerr += 1 ))
        fi
    fi
    [[ $report =~ *pot* ]] && add_to $error; (( nerr += 1)) 
done
#How much data is on the pots?
add_to "${newline}--- Free space on various machines ---"
for oe in {0..0}
do
    dfline=`ssh pot${oe} "df -h | grep data"`
    add_to "pot${oe}:${newline}${dfline}"
    p=`echo $dfline | awk '{print $5}'`
    [ ${p%?} -gt 95 ] && add_to $error; (( nerr += 1 ))
done
dfline=`ssh paper1 "df -h | grep data"`
add_to "paper1:${newline}${dfline}"
p=`echo $dfline | awk '{print $5}'`
[ ${p%?} -gt 95 ] && add_to $error; (( nerr += 1 ))

#Find Lost Lambs?
add_to "${newline}--- Lost Lambs? ---"
for oe in {0..1}
do
    latest=`ssh pot${oe} "ls -td /data${oe}/psa66*"` 
    latest=`echo $latest | cut -f 1 -d " "`
    cnt=0
    add_to " Latest mod to pot${oe} is on ${latest}"
    for f in `ssh pot${oe} "ls -d ${latest}/zen*uv"`
    do
        ssh pot${oe} test -e ${f}cRRE \
            || (( cnt += 1 )); echo $f >> /home/obs/MissingFilesLogs/${today}.log
    done
    add_to "${cnt} missing files."
    [ $cnt != 0 ] && add_to "A list can be found in qmaster:/home/obs/MissingFilesLogs/${today}.log"
done

#send a report.
add_to "${newline}This daily report is located in qmaster:/home/obs/AutomaticallyGeneratedReports/${today}.log"

echo "$report" >> /home/obs/AutomaticallyGeneratedReports/${today}.log
[ $nerr -gt 0 ] && echo "$report" | mail -s "Daily Data Checkup" $mailto
echo "$report" | mail -s "Daily Data Report" $paperblog 
