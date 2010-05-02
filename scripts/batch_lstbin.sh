#$ -S /bin/bash
DLST=.3
LSTS=`python -c "import numpy as n; print ' '.join(['%f_%f' % (a, a+${DLST}) for a in n.arange(0,6.2832,${DLST})])"`
MYLSTS=`pull_args.py $LSTS`
for lst in $MYLSTS ; do
    echo lstbin.py -p yy -a "-1,1_1" -l $lst 
    lstbin.py -p yy -a "-1,1_1" -l $lst $*
done
