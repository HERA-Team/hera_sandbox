#$ -S /bin/bash
ARGS=`pull_args.py $*`

for FILE in $ARGS; do
    #xrfi.py -m val ${FILE} 
    /data1/paper/arp/scripts/xrfi_simple.py --dt=3.5 --combine -t 80 ${FILE}
done