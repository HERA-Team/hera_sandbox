#$ -S /bin/bash
#$ -j y
#$ -o grid_output
#$ -e grid_output
#$ -cwd
#$ -V
ARGS=`pull_args.py $*`
CORRECT=correct_psa746_v002.py
CAL=psa746_v010
ANTS="all"
#WINDOW='blackman-harris'
WINDOW='hamming'

for FILE in $ARGS; do
    #$CORRECT $FILE -t /data3/paper/psa/Jul2011/psa746
    #/data3/paper/pober/capo/jcp/scripts/tempgain.py ${FILE}c
    #/usr/global/paper/bin/xrfi_simple.py -a 1 --combine -t 20 -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 ${FILE}
    #/usr/global/paper/bin/apply_cal.py -C $CAL ${FILE}R
    #pspec_prep_fg.py -C $CAL --clean=1e-9 --horizon=50 --window=$WINDOW --model --nophs --nogain ${FILE}
    #pspec_prep_fg.py -C $CAL --clean=1e-9 --horizon=50 --window=$WINDOW --nophs --nogain ${FILE}
    delay_null.py -C $CAL --clean=1e-9 --horizon=50 --window=$WINDOW --nophs --nogain ${FILE}
    #/usr/global/paper/bin/xrfi_simple.py -n 3 -c 0_350,1760_2047 ${FILE} #needs >2G of memory for uncompressed data
    #/data3/paper/pober/capo/dfm/scripts/xtalk4.py ${FILE}R
done
