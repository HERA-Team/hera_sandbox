#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=2G
#$ -j y
#$ -N lstbin
#$ -o grid_output

#usage 
# batch_lstbin.sh mylstbinrun.cfg
#to be defined
#export START_LST= <hours>
#export END_LST=  <hours>
#export DLST= <hours>
#export CAPO=/path/to/capo/
#export CAL=psaXXX_v00X
#export LST_RES=<seconds>
#export TFILE=<length of file in seconds>
#export FILES=your data ex `ls -d /data4/raw_data/psa903-930/psa*/*RREXCBR`

( #the parens keep the vars defined in the config local
echo HELLO.
echo RUNNING YOUR CFG SCRIPT
. $* #execute the cfg file as a script
shopt -s extglob
. /usr/global/paper/paperenv.sh
workon dcj-PAPER
~/scripts/pywhich psa898_v003
myLSTS=`${CAPO}/dcj/scripts/pull_args.py ${LSTS}`
for LST in $myLSTS
do
echo Working on $LST
LSTFILES=`${CAPO}/pspec_pipeline/lst_select.py -C ${CAL} --ra=${LST} --suntime=n ${FILES}`
echo working on  `echo ${LSTFILES} | wc -w`/ `echo ${FILES} | wc -w ` files
echo ${CAPO}/jcp/scripts/lstbin_v02.py -C ${CAL} -a 0_16 -p I ${LSTFILES}
${CAPO}/jcp/scripts/lstbin_v02.py -C ${CAL} -a 0_16 -p I ${LSTFILES}
done
)
