#!/bin/bash

seps='sep0,1 sep1,1 sep-1,1'
cal=psa6240_v003
frpad='1.05'
#alietal='--alietal'
alietal=''
out='/home/mkolopanis/psa64/lstbin_psa64_noise_only_30Jy_frpad_1.05'
appelation='uvGAs'
declare -A ants
ants[sep0,1]=0_44
ants[sep1,1]=1_3
ants[sep-1,1]=1_48

paths='even odd'

printf 'This is Batch Fringe Rate Filter using:\n'
printf 'calfile: %s \n' $cal
printf 'frpad: %s\n' $frpad
sleep 1.5

if [ ! -d $out   ]; then
    printf 'This output director appears to be new\n'
    printf 'Creating Directoires now\n'
    for path in $paths; do
        for sep in $seps; do
            mkdir -p ${out}/${path}/${sep}
        done
   done
fi

for path in $paths; do
    for sep in $seps; do
        printf 'Checking for Old FRF Data\n'
        outpath=$out/$path
        if [[ $(ls -d $outpath/$sep/*.${appelation}L) ]] ; then
            printf 'Found %s folders\n' $(ls -d $outpath/$sep/*.${appelation}L| wc -l)
            printf 'First deleting old FRF Filtered Data \n'
            sleep 1
            printf 'I wil delete: \n'
            sleep 1
            printf '%s \n' $(ls -d $outpath/$sep/*.${appelation}L)
            read -p "Are you sure you want to delete these? [y/N] " -n 1 -r
            printf '\n'
            if [[ $REPLY =~ ^[Yy]$ ]]
                then
                    rm -rf $(ls -d $outpath/$sep/*.${appelation}L)
                    printf 'Files are deleted\nContinuting to Fringe Rate Filter\n'
                else
                    printf 'Nothing is deleted\ncontinuing\n'
                    continue
            fi
        fi
        if [[ -z ${alietal} ]]; then
        printf 'Filtering %s by selecting ant %s \n' $sep ${ants[$sep]}
        printf '/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py -C %s   -a %s --frpad %s --outpath=%s/'  $cal ${ants[$sep]} $frpad ${out}
        printf '%s\n' $path/$sep
        #files=$(ls -d $path/$sep/lst.*242.[3456]*.${appelation})
        files=$(ls -d $path/$sep/*.${appelation})
        "/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal $alietal --frpad $frpad $files --outpath=${out}

        else
        printf 'Filtering %s by selecting ant %s \n' $sep ${ants[$sep]}
        printf '/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py -C %s  %s -a %s --frpad %s --outpath=%s/'  $cal $alietal  ${ants[$sep]} $frpad ${out}
        printf '%s\n' $path/$sep
        #files=$(ls -d $path/$sep/lst.*242.[3456]*.uvGA)
        files=$(ls -d $path/$sep/*.${appelation})
        "/home/mkolopanis/src/capo/mjk/scripts/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal $alietal --frpad $frpad $files --outpath=${out}
        fi
    done
done
