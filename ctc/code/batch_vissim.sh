#$ -S /bin/bash
#$ -V
#$ -cwd 
#$ -l h_vmem=12G
#$ -l paper
#$ -o /data2/home/cacheng/capo/ctc/code/vissim.out
#$ -e /data2/home/cacheng/capo/ctc/code/vissim.err

echo CODE START

startjd=2454500
endjd=2454501
inttime=5000
numchunks=4 #number of processors
nchunks=`seq ${numchunks}`
echo nchunks: ${nchunks}

ARGS=`python -c "import numpy; import aipy; print ' '.join(map(str,numpy.arange(${startjd},${endjd},${inttime}/aipy.const.s_per_day)))"`
echo times: ${ARGS}

numtimes=`python -c "import numpy; import aipy; print len(numpy.arange(${startjd},${endjd},${inttime}/aipy.const.s_per_day))"`
echo numtimes: ${numtimes}

startindex=1

chunksize=$((${numtimes}/${numchunks}+1))
echo numchunks: ${numchunks}
echo chunksize: ${chunksize}

loopargs=`pull_args.py ${nchunks}`

for chunk in ${loopargs}; do

    #echo ${ARGS} | cut -d " " -f ${startindex}-$((${startindex}+${chunksize}-1)) #list of times broken into chunks
    list="$(echo ${ARGS} | cut -d " " -f ${startindex}-$((${startindex}+${chunksize}-1)))" #list of times broken into chunks
    echo ${list}
    name=`echo ${list} | cut -d " " -f 1` 
    echo ${name}
    #vis_simulation_v4.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime 20000 --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec203lmax200/ --filename test.uv -C psa898_v003 -a 0_16 ${list}

    
    startindex=$((${startindex}+${chunksize}));
 
done


