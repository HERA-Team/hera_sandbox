DIR=`pull_args.py $*`
for dir in $DIR; do 
    echo working on $dir 
    for f in $dir/*.uvcRREcAz; do
        fname=`du -h $f | awk '{if ($1 ~ 1.1G) print $2}'`
        if [ $fname ]; then 
            echo xtalk3.py  $fname
            xtalk3.py $fname
        fi
    done
done
