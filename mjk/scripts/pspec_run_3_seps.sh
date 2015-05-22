VERSION=$*

for i in -1 0 1; do
    mkdir sep${i},1/v0${VERSION}
    cd sep${i},1/v0${VERSION}
    python /home/mkolopanis/src/capo/mjk/scripts/pspec_cov_v002.py -C psa6240_v003 --window=none -c 95_115 --sep=sep${i},1 -b 20
    #python /home/mkolopanis/src/capo/mjk/scripts/pspec_cov_v002.py -C psa6240_v003 -c 95_115 --sep=sep${i},1 -b 100 -a cross
#this rmbls was for the no omnical data. baseline was just bad.
#    python /Users/sherlock/src/capo/zsa/scripts/pspec_cov_v002.py -C psa6240_v003 --window=none -c 95_115 --sep=sep${i},1 -b 100 --rmbls=13111
#    python /Users/sherlock/src/capo/zsa/scripts/pspec_cov_v002.py -C psa6240_v003 --window=none -c 95_115 --sep=sep${i},1 -b 20
    cd ../../
done

mkdir v0${VERSION}
cd v0${VERSION}

#~/src/capo/zsa/scripts/pspec_cov_boot.py ../sep{-1,0,1},1/v0${VERSION}/pspec_boot00*.npz
#python ../pspec_cov_boot.py ../sep{0},1/v0${VERSION}/pspec_boot00*.npz --nocov
/home/mkolopanis/src/capo/arp/scripts/pspec_cov_boot.py ../sep{-1,0,1},1/v0${VERSION}/pspec_boot00*.npz

#pspec_plot_pk_k3pk.py pspec.npz --show --flux -f 1.6 --cov 
/home/mkolopanis/src/capo/mjk/scripts/plot_pk_k3pk_zsa_2.py pspec.npz --show --cov --afrf

#vi README
cd ../
