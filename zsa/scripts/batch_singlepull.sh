DIR=$*
ANT_LIST=49,41,47,19,29,28,34,51,10,3,25,48,24,55,27,57,9,58,1,4,17,13,56,59,22,61,35,18,5,32,30,23
for dir in $DIR; do
    echo working in $dir
    for f in $dir/*.uvcRREc; do
        #echo pull_antpols.py -p xx -a "($ANT_LIST)_($ANT_LIST)" $file 
        pull_antpols.py -p xx -a "($ANT_LIST)_($ANT_LIST)" $f || rm -rf $f && correct_psa6240.py ${f%?} && pull_antpols.py -p xx -a "($ANT_LIST)_($ANT_LIST)" $f
    done;
done;
