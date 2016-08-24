#! /bin/bash
workon PAPER
DATE=`date +"%d-%b-%y"`
DATAROOT=/data3/PAPER/psa128/
#get the latest data dir
DATADIR=`ls -d ${DATAROOT}/2456* | tail -n 1`
EMAIL=teampaper.stats@blogger.com
echo found $DATADIR
#check to see if RFI has been summarized
ls ${DATADIR}/*RFI*png 2>/dev/null
if [ $? -ne 0 ]
then 
#if theres no png we need to make one!
cd ${DATADIR}
python ~/scripts/dailyRFIreport.py *npz
PLOT=`ls -t ${DATADIR}/*RFI*png | head -n 1`
echo ${PLOT}
DATE=`python -c "print '${PLOT}'.split('_')[3].split('.')[0]"`
echo mutt  -s "RFI report for "${DATE} -a ${PLOT} -- ${EMAIL}
echo  | mutt  -s "RFI report for "${DATE} -a ${PLOT} -- ${EMAIL}

fi
