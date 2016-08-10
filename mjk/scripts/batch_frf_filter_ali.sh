#!/bin/bash

seps='sep0,1 sep1,1 sep-1,1'
cal=psa6240_v003
frpad='1.0'
alietal='--alietal'
#alietal=''
declare -A ants
ants[sep0,1]=0_44
ants[sep1,1]=1_3
ants[sep-1,1]=1_48

paths='even odd'

printf 'This is Batch Fringe Rate Filter using:\n'
printf 'calfile: %s \n' $cal
printf 'frpad: %s\n' $frpad
printf 'alietal: %s\n' $alietal
sleep 1.5
for path in $paths; do
    for sep in $seps; do
        printf 'Checking for Old FRF Data\n'
        outpath='/home/mkolopanis/psa64/lstbin_psa64_ali_reconstruction/'$path
        if [[ $(ls -d $outpath/$sep/*.uvGAL) ]] ; then
            printf 'Found %s folders\n' $(ls -d $outpath/$sep/*.uvGAL| wc -l)
            printf 'First deleting old FRF Filtered Data \n'
            sleep 1
            printf 'I wil delete: \n'
            sleep 1
            printf '%s \n' $(ls -d $outpath/$sep/*.uvGAL)
            read -p "Are you sure you want to delete these? [y/N] " -n 1 -r
            printf '\n'
            if [[ $REPLY =~ ^[Yy]$ ]]
                then
                    rm -rf $(ls -d $outpath/$sep/*.uvGAL)
                    printf 'Files are deleted\nContinuting to Fringe Rate Filter\n'
                else
                    printf 'Nothing is deleted\ncontinuing\n'
                    continue
            fi
        fi
        printf 'Filtering %s by selecting ant %s \n' $sep ${ants[$sep]}
        if [[ $path == 'even' ]]
        then
            outpath='/home/mkolopanis/psa64/lstbin_psa64_ali_reconstruction'
            printf '/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py -C %s  %s -a %s --frpad %s --outpath=%s '  $cal $alietal ${ants[$sep]} $frpad $outpath
            printf '%s\n' $path/$sep
            #files=$(ls -d $path/$sep/lst.*242.[3456]*.uvGA)
            files=$(ls -d $path/$sep/*.uvGA)
            "/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal $alietal --frpad $frpad --outpath=$outpath $files 
        else
            outpath='/home/mkolopanis/psa64/lstbin_psa64_ali_reconstruction'
            printf '/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py -C %s  %s -a %s --frpad %s --outpath=%s '  $cal $alietal ${ants[$sep]} $frpad $outpath
            printf '%s\n' $path/$sep
            #files=$(ls -d $path/$sep/lst.*243.[3456]*.uvGA)
            files=$(ls -d $path/$sep/*.uvGA)
            "/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal $alietal --frpad $frpad --outpath=$outpath $files
        fi
    done
done
