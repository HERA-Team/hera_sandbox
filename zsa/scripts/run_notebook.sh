#! /bin/bash 

export HERA_DATA_DIR=$1
jd=`python -c "print '${HERA_DATA_DIR}'.split('/')[-1]"`
WORKDIR=${PWD}
FILEROOT=firstcal_notebook_${jd}
#FILEROOT=omnical_notebook_${jd}
BASENOTEBOOK=firstcal_notebook.ipynb
#BASENOTEBOOK=omnical_notebook.ipynb
export FIRSTCAL_SCRIPT=~/src/mycapo/omni/firstcal.py
export OMNI_APPLY_SCRIPT=~/src/mycapo/omni/omni_apply.py
export OMNI_RUN_SCRIPT=~/src/mycapo/omni/omni_run.py
export EX_ANTS_X=81
export EX_ANTS_Y=81
export HERA_CAL_FILE=hsa7458_v000
export REDO_FIRSTCAL=0 #0 for dont redo, > 0 for redo
export REDO_REWRITE=0 #
export REDO_FIRSTCAL_APPLY=0 #
export REDO_OMNIRUN=0
export REDO_OMNIAPPLY=0
export SLACK_KEY=''  # slack key for messages

echo HERA_DATA_DIR=$HERA_DATA_DIR
echo REDO_FIRSTAL=$REDO_FIRSTCAL
echo REDO_REWRITE=$REDO_REWRITE
echo REDO_FIRSTCAL_APPLY=$REDO_FIRSTCAL_APPLY
echo REDO_OMNIRUN=$REDO_OMNIRUN
echo REDO_OMNIAPPLY=$REDO_OMNIAPPLY
echo WORKDIR=$WORKDIR
echo FILEROOT=$FILEROOT
echo BASENOTEBOOK=$BASENOTEBOOK
echo FIRSTCAL_SCRIPT=$FIRSTCAL_SCRIPT
echo OMNI_RUN_SCRIPT=$OMNI_RUN_SCRIPT
echo OMNI_APPLY_SCRIPT=$OMNI_APPLY_SCRIPT
echo EX_ANTS_X=$EX_ANTS_X
echo EX_ANTS_Y=$EX_ANTS_Y
echo HERA_CAL_FILE=$HERA_CAL_FILE

echo converting notebook: 
echo python -c "import capo; capo.zsa.run_nb('${WORKDIR}', '${FILEROOT}', '${BASENOTEBOOK}')"
python -c "import capo; capo.zsa.run_nb('${WORKDIR}', '${FILEROOT}', '${BASENOTEBOOK}')"
echo Finished converting notebook

echo Moving notebook to git repo
echo mv firstcal_notebook_${jd}.ipynb ~/HERA_plots/.
mv firstcal_notebook_${jd}.ipynb ~/HERA_plots/.
echo cd ~/HERA_plots
cd ~/HERA_plots
echo git add firstcal_notebook_${jd}.ipynb
git add firstcal_notebook_${jd}.ipynb
echo git commit -m "firstcal notebook for ${jd}"
git commit -m "firstcal notebook for ${jd}"

python -c "import slacker; slack=slacker.Slacker(${SLACK_KEY}); slack.chat.post_message('@zaki', 'Firstcal notebook for %s finished running'%'${jd}', '@zaki')"
