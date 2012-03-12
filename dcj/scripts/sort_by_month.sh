#! /bin/bash
# Move listed files to monthly directories
for F in $*; do
    echo Moving $F
    MONTH=`month.py $F`
    if ! ls $MONTH &> /dev/null; then
        echo Creating $MONTH
        mkdir $MONTH
    fi
    mv $F $MONTH/

done
