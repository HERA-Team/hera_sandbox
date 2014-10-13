for file in /home/zakiali/psa_live/*; do
    if [ -d "$file" ]; then
 #       echo $file
        count=`ls -d $file/*.uvcRRE | wc -l`
#        echo $count
        if [ $count -eq 70 ]; then
            echo $file
        fi
    fi 
done    
