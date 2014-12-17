DIRS=$*

#echo start > bad_pullant_files.txt
for dir in $DIRS; do
    echo $dir
    one=`ls -d $dir/*.uvcRREc | wc -l`
    two=`ls -d $dir/*.uvcRREcA | wc -l`
    if [ $one -eq $two ] 
        then echo TRUE! $one equals $two
        else 
            echo $dir is bad $one not equal to $two
            for file in $dir/*.uvcRREc; do
                if [ ! -d $file'A' ] 
#                    then echo $file >> bad_pullant_files.txt
                    then echo \t $file 
                fi
            done
    fi
done
     

