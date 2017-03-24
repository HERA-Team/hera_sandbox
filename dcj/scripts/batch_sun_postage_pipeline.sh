#$ -S /bin/bash
#. ~/.bashrc
#files to work on. Can be newest or all or something.
#$ -j y
#$ -N snap_Sun
#$ -o /data1/paper/gb/snapshots/grid_log/
CHN=70_75 #noisy?
#CHN=150_160
ANTS=cross,-0_5,-0_6,-1_8,-3_11,-3_12,-4_15,-5_6,-12_13
#ANTS=cross
SZ=400
RES=2.8
#RES=2.8
ALTMIN=-20
#ALTMIN=30
UNIF=.01
AIPY=/usr/global/paper/bin/
PYTHONPATH=/usr/global/paper/lib/python2.6/site-packages/:$PYTHONPATH
LSTS=`lst_select.py -C pgb966_v003 --sun -d - $*|pull_days.py`
if [ $? -ne 0 ]; then exit;fi


echo "rough_pipelining: "$LSTS
rough_sun_pipeline.sh $LSTS
if [ $? -ne 0 ]; then exit;fi

CALD=`echo $LSTS | sed 's/uvcb/uvcbdrx.e_Sunm/g'`
echo "Imaging: "$CALD
DAY=`python -c "import sys; print \"$CALD\".split('.')[1]"`
echo "DAY="$DAY
${AIPY}mk_img.py -c $CHN -C pgb966_v003 -p yy\
 -s Sun --fmt=snap_Sun_${DAY}_%01d \
 --size=$SZ --res=$RES --altmin=$ALTMIN -u $UNIF -a $ANTS  $CALD
if [ $? -ne 0 ]; then exit;fi

${AIPY}cl_img.py -d cln snap_Sun_${DAY}_0.d[i,b]m.fits
if [ $? -ne 0 ]; then exit;fi

mv snap*.fits snapshots/
