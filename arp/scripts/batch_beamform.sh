#$ -S /bin/bash
SRCS=`cat srclist.txt`
SRCS=`pull_args.py ${SRCS}`

for SRC in $SRCS ; do
    echo beamform.py -C psa746_v008 -a cross,-40,-55 --minuv=20 -p xx --cat=culgoora,parkes,misc -s $SRC $*
    beamform.py -C psa746_v008 -a cross,-40,-55 --minuv=20 -p xx --cat=culgoora,parkes,misc -s $SRC $*
    echo lstbin.py --nogaps -C psa746_v008 --lst_res=85.9 *bm_${SRC}
    lstbin.py --nogaps -C psa746_v008 --lst_res=85.9 *bm_${SRC}
    echo mv lst.2455*.uv lst747_${SRC}.uv
    mv lst.2455*.uv lst747_${SRC}.uv
    echo fringe_rate_filter.py --maxfr=.0001 --minfr=-.0001 --clean=1e-6 lst747_${SRC}.uv
    fringe_rate_filter.py --maxfr=.0001 --minfr=-.0001 --clean=1e-6 lst747_${SRC}.uv
done
