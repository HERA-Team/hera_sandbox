#$ -S /bin/bash
#$ -j y
#$ -N lst_bin
#$ -o /data1/paper/jacobsda/stokes/pgb050/grid_output
DLST=.3
LSTS=`python -c "import numpy as n; print ' '.join(['%f_%f' % (a, a+${DLST}) for a in n.arange(0,2*n.pi,${DLST})])"`
MYLSTS=`pull_args.py $LSTS`

if ! ls /scratch/paper &> /dev/null; then
    mkdir /scratch/paper
fi
LOCAL=/scratch/paper/pgb050/
if ! ls $LOCAL &> /dev/null; then
    echo Creating local data dir...
    mkdir $LOCAL
#    mv /tmp/lst*uv $LOCAL
fi
TARGS=
INP=$*
#DATA_ROOT=`python -c "print '${INP}'.split(' ')[0]"`
#DATA_ROOT=`python -c "import string as s; print s.join('${DATA_ROOT}'.split('/')[:-1],'/')+'/'"`
#echo Moving data from: $DATA_ROOT
for lst in $MYLSTS; do
    lst_h=`python -c "import string as s; print str(float('${lst}'.split('_')[0])*12/3.14159)+'_'+ str(float('${lst}'.split('_')[1])*12/3.1415)"`
    echo finding files in lst range $lst_h "(hours)"
 #   FS=`lst_select.py --ra=${lst_h} -C pgb015_v004 $*`
    FS=$*
    echo Preloading: ${FS}
    for F in $FS; do
        rsync -avz shredder:${F} $LOCAL
        TARGS=${LOCAL}`python -c "print '${F}'.split('/')[-1]"`" "${TARGS}
    done
 #   echo -e processing ${TARGS}
    echo lstbin.py -p yy,xx,xy,yx -a "all,-1" -l $lst $TARGS
    /data1/paper/jacobsda/scripts/lstbin.py -p yy,xx,yx,xy -a "all,-1" -l $lst $TARGS
done
#remove duplicate files from the list
#TARGS = `python -c "print ' '.join(set('${TARGS}'.split()))"`
#echo -e processing ${TARGS}
#for lst in $MYLSTS ; do
#    echo lstbin.py -p yy,xx,xy,yx -a "all,-1" -l $lst $TARGS 
#    lstbin.py -p yy,xx,yx,xy -a "all,-1" -l $lst $TARGS
#done
