#given a list of input files and a final postfix (FINAL)
# find and remove any intermediate steps.
#  but only if the final product is done.
INTERMEDIATES=(B)
FINAL=BR
INPUTS=$*
unset DELETED
unset INCOMPLETE
unset COMPLETE
for FILE in $INPUTS
do
#check if the end file exists
    if [[ -e ${FILE}${FINAL} ]]
        then 
        #echo Found final product ${FILE}${FINAL}
	COMPLETE=("${COMPLETE[@]}" ${FILE}${FINAL})
	#cleanup!
	    for PF in $INTERMEDIATES
	    do
	        if [[ -e ${FILE}${PF} ]]
	        then
                    #echo Cleaning up ${FILE}${PF}
                    rm ${FILE}${PF}/*
                    rmdir ${FILE}${PF}
	            DELETED=("${DELETED[@]}" ${FILE}${PF})
                fi
	    done
    else
            INCOMPLETE=("${INCOMPLETE[@]}" $FILE)
    fi
done
echo "Found " ${#COMPLETE[@]} "files"
echo "Deleted" ${#DELETED[@]} "files"
echo "Found" ${#INCOMPLETE[@]} "incomplete files"
for FILE in ${INCOMPLETE[@]}
do
    echo $FILE
done
