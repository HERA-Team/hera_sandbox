#$ -S /bin/bash
ARGS=`pull_args.py $*`

for ARG in $ARGS ; do
    rm_npz9.py -C pgb015_v005 $ARG
    filter_src.py -C pgb015_v005 -s Sun -r 5 -d 5 --clean=1e-3 ${ARG}d
done
