DIRS=$*

echo start > gooddir.txt
echo start > baddir.txt
for dir in $DIRS; do
    echo $dir
    one=`ls -d $dir/*.uvcRRE | wc -l`
    two=`ls -d $dir/*.uvcRREc | wc -l`
    if [ $one -eq $two ] 
        then echo TRUE! $one equals $two
        echo $dir >> gooddir.txt
        else 
            echo $dir is bad $one not equal to $two
            echo $dir >> baddir.txt
    fi
done
     

