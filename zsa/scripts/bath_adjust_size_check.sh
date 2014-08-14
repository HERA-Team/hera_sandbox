DIRS=$*

#echo start > bad_pullant_files.txt
for dir in $DIRS; do
    echo $dir
    one=`ls -d $dir/*.uvcRREcAz | wc -l`
    twelve=12
    if [ $one -eq $twelve ]
        then echo $dir has $one z files. Good.
        else
            echo $dir has $one z files. Bad.
    fi
done


#    two=`ls -d $dir/*.uvcRREcA | wc -l`
#    if [ $one -eq $two ] 
#        then echo TRUE! $one equals $two
#        else 
#            echo $dir is bad $one not equal to $two
#            for file in $dir/*.uvcRREc; do
#                if [ ! -d $file'A' ] 
#                    then echo $file >> bad_pullant_files.txt
#                    then echo \t $file 
#                fi
#            done
#    fi
#done

