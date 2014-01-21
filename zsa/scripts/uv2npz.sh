CALFILE=psa746_v008
ARGS=`pull_args.py $*`

for FILE in $ARGS; do
    SRC=`echo ${FILE} | cut -d _ -f2 | sed 's/RmI//'`
    echo $SRC
    beamform_uv2npz.py ${FILE} -C ${CALFILE} -s ${SRC} -p xx --cat=helm,misc,culgoora,mrt
#    beamform_uv2npz.py ${FILE} -C ${CALFILE} -s ${SRC} -p xx
#    beamform_uv2npz.py ${FILE} -C ${CALFILE} -s ${SRC} -p xx --cat=mrt
done
