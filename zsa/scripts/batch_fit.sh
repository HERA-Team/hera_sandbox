echo start > fittest.txt
bl=49_41

for i in -4 0 2.7 5 9; do
    echo fitmdl.py -C psa6240_v001 -p xx -a $bl -s pic -c 50_150_10 -x 20 -P "aa=tau_ew/$i" $*
    fitmdl.py -C psa6240_v001 -p xx -a $bl -s pic,for,crab -c 50_150_10 -x 20 -P "aa=tau_ew/$i" $* | tee -a fittest.txt
    
done
