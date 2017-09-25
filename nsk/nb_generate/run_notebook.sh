#!/bin/bash 
# run_notebook.sh <hera_data_directory>

# assign data directory
export DATA_PATH=$1

# get JD
jd=`python -c "print '${DATA_PATH}'.split('/')[-1]"`

# get more env vars
WORKDIR=${PWD}
BASENBDIR=$HOME/HERA_plots
export CALFILE=hsa7458_v001
OUTPUT=data_inspect_"$jd".ipynb
OUTPUTDIR=$HOME/HERA_plots

# copy and run notebook
jupyter nbconvert --output=$OUTPUTDIR/$OUTPUT --to notebook --execute $BASENBDIR/Data_Inspect.ipynb

# cd to git repo
cd $OUTPUTDIR

# add to git repo
git add $OUTPUT
git commit -m "data inspect notebook for $jd"
git push