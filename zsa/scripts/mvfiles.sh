DEST=/data4/paper/2012EoR/psa_live
DIRS=$*

for dir in $DIRS; do 
    for uv in $dir/*.uvcRRE; do
        if [ ! -d "${DEST}/${uv}" ] ; then
            echo cp -r ${uv} ${DEST}/${dir}/ 
            cp -r  ${uv} ${DEST}/${dir}/ 
            #echo copied ${uv} to ${DEST}/${dir}/ >> ~/copied_uvfiles_from_lost_lambs_to_psa_live_10_11_2013
            echo copied ${uv} to ${DEST}/${dir}/ >> ~/temp
        fi
    done;
done;
