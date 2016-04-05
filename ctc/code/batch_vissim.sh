#$ -S /bin/bash
#$ -V
#$ -cwd 
#$ -l h_vmem=16G
#$ -l paper
#$ -o /data2/home/cacheng/capo/ctc/code/gridoutput
#$ -e /data2/home/cacheng/capo/ctc/code/gridoutput

myargs=`pull_args.py $*`

echo my times: ${myargs}

name=`echo ${myargs} | cut -d " " -f 1`
echo first arg: ${name}

echo vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.14679803 --nchan 20 --inttime 42.9 --map pspec --mappath /home/cacheng/capo/ctc/images/pspecs/pspec20lmax300/ --filename /home/cacheng/capo/ctc/tables/20files_samemaps/pspec_${name}.uv -C psa898_v003 -a 0_16 ${myargs}

#vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.14679803 --nchan 20 --inttime 42.9 --map pspec --mappath /home/cacheng/capo/ctc/images/pspecs/pspec20lmax300/ --filename /home/cacheng/capo/ctc/tables/20files_samemaps/pspec_${name}.uv -C psa898_v003 -a 0_16 ${myargs}

vis_simulation_v4.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime 20000 --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec100lmax100/ --filename test_${name}.uv -C psa898_v003 -a 0_16 ${myargs}




#startjd=2454500
#endjd=2454501
#inttime=20000
#numchunks=4 #number of processors... MUST match qsub command

#ARGS=`python -c "import numpy; import aipy; print ' '.join(map(str,numpy.arange(${startjd},${endjd},${inttime}/aipy.const.s_per_day)))"`
#echo times: ${ARGS}

#numtimes=`python -c "import numpy; import aipy; print len(numpy.arange(${startjd},${endjd},${inttime}/aipy.const.s_per_day))"`
#echo numtimes: ${numtimes}

#chunksize=$((${numtimes}/${numchunks}+1))
#echo numchunks: ${numchunks}
#echo chunksize: ${chunksize}

#timeindex=`python -c "import numpy; print ' '.join(map(str,numpy.arange(1,${numtimes}+1,${chunksize})))"`

#loopargs=`pull_args.py ${timeindex}`

#for chunk in ${loopargs}; do
    
    #echo ${chunk}
    #list="$(echo ${ARGS} | cut -d " " -f ${chunk}-$((${chunk}+${chunksize}-1)))" #list of times broken into chunks
    #echo mylist: ${list}
    #name=`echo ${list} | cut -d " " -f 1` 
    #echo ${name}
    #echo vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.1 --nchan 203 --inttime ${inttime} --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec203lmax200/ --filename /data2/home/cacheng/capo/ctc/tables/203files/pspec_${name}.uv -C psa898_v003 -a 0_16 ${list}
    #vis_simulation_v4.py --sdf 0.00049261 --sfreq 0.1 --nchan 203 --inttime ${inttime} --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec203lmax200/ --filename /data2/home/cacheng/capo/ctc/tables/203files/pspec_${name}.uv -C psa898_v003 -a 0_16 ${list}
    #vis_simulation_v4.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime ${inttime} --map pspec --mappath /data2/home/cacheng/capo/ctc/images/pspecs/pspec100lmax100/ --filename test_${name}.uv -C psa898_v003 -a 0_16 ${list}
#done


