
DIR=$*
for dir in $DIR; do
    echo working in $dir
    for f in $dir/*.uvcRRE; do
        #echo pull_antpols.py -p xx -a "($ANT_LIST)_($ANT_LIST)" $file 
        correct_psa6240.py $f || echo $f cannot be corrected
    done;
done;
