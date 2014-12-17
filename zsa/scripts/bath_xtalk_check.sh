DIRS=$*

#echo start > bad_pullant_files.txt
for dir in $DIRS; do
    echo $dir
    one=`ls -d $dir/*.uvcRREcAz | wc -l`
    two=`ls -d $dir/*.uvcRREcAzx | wc -l`
    if [ $one -eq $two ] 
        then echo TRUE! $one equals $two
        else 
            echo $dir is bad $one not equal to $two
            for file in $dir/*.uvcRREcAz; do
                if [ ! -d $file'x' ] 
#                    then echo $file >> bad_pullant_files.txt
                    then echo \t $file 
                fi
            done
    fi
done
     

