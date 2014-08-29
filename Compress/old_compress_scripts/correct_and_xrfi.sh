#$ -S /bin/bash
#$ -V
#$ -l h_vmem=8G
#$ -cwd

ARGS=`pull_args.py /home/obs/Scratch/zen.*.uv`
for FILE in $ARGS; do
    echo "===> correct and first pass XRFI"
    echo correct_and_XRFI.py -a 1 -t 80 -n 5 --df=6 -c 0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023 ${FILE}
    time correct_and_XRFI.py -a 1 -t 80 -n 5 --df=6 -c 0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023 ${FILE}
done
