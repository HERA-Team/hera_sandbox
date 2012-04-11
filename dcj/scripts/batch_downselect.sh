#$ -S /bin/bash
#$ -V
#$ -cwd 
#$ -o grid_output
#$ -j y
#$ -N PSAds
#$ -l h_vmem=500M


ARGS=`pull_args.py $*`
RFI_CHAN=0_130,755_777,1540,1704,1827,1868,1901_2047
for FILE in $ARGS
do
    #pull out the shortest EW spacing and slanty EW
    /home/jacobsda/capo/dcj/scripts/pull_gridspacing.py --correct $FILE

    #pull grid spacing outputs to the local dir, so strip the path from the input filenames
    WORKFILE=`python -c "import os;print os.path.basename('$FILE')"`

    #initial pass at xrfi. we are doing a --combine and outputting the npz instead of a flagged dataset.
    /data3/paper/pober/capo/jcp/scripts/xrfi_simple.py -a cross --dt=4 --df=6 -c $RFI_CHAN --combine -t 15  --npz ${WORKFILE}'G'

done
