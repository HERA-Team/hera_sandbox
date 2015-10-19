#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=8G
#$ -l paper
#$ -N omni_v1

ARGS=`pull_args.py $*`

echo Processing: ${ARGS}

for f in ${ARGS}; do

    echo working on ${f}...
    
    TAG=${f##*/}
    TAG=${TAG%.*}
    TAG=${TAG%.*} #gets "zen.2456942.38963" from "zen.2456942.38963.xx.uvcRRE", for example
    
    #WITHOUT XTALK REMOVAL
    #echo omnical_PSA128.py -p xx -C psa6622_v003 -d ${TAG} -i /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/redundantinfo_first_cal_epoch3xx.bin -r /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3xx.p -o /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v1_xtalk/ -u -s ${f}
    #omnical_PSA128.py -p xx -C psa6622_v003 -d ${TAG} -i /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/redundantinfo_first_cal_epoch3xx.bin -r /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3xx.p -o /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v1_xtalk/ -u -s ${f}

    #WITH XTALK REMOVAL
    echo omnical_PSA128.py -p xx -C psa6622_v003 --add -d ${TAG} -i /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/redundantinfo_first_cal_epoch3xx.bin -r /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3xx.p -o /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v1_xtalk/ -u -s ${f}
    omnical_PSA128.py -p xx -C psa6622_v003 --add -d ${TAG} -i /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/redundantinfo_first_cal_epoch3xx.bin -r /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3xx.p -o /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v1_xtalk/ -u -s ${f}

done

