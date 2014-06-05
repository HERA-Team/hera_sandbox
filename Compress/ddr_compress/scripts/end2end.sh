#! /bin/bash
#copy the test data to the pot
echo "initializing the pot for the test"
scp -R  ~/test_pot/sim4/z*uv pot0:
echo "starting the task servers on still machines"
echo "starting qdaemon"
echo "loading the data into the db"
ssh pot0 add_observations.py /data/z*uv  #NB this step won't work till ddr_compress is installed on the pot

