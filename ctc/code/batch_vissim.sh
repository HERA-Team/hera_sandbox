#$ -S /bin/bash
ARGS = `pull_args.py $*` #give a chunk of times

vis_simulation_v4.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime 20000 --map pspec --mappath /Users/carinacheng/capo/ctc/images/pspecs/pspec100lmax100/ --filename test.uv -C psa898_v003 -a 0_16 $ARGS



