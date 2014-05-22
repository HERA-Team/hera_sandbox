#! /bin/bash
# Tests wait_for_still_process_start.py
echo test_wait_for_still_process_start.sh
echo WARNING: closing out all open tasks. If qdaemon is running, you probably just made a mistake 
/home/obs/Compress/close_out_open_tasks.py
testfile="qmaster:"`ls -d /home/obs/Compress/tests/minimumtest/z*uv | head -n 1`
echo "using test filename:"$testfile
echo starting wait_for_still_process_start.py
/home/obs/Compress/wait_for_still_process_start.py --operation="TEST" qmaster:/dev/null --timeout=0.5 && echo "process started"&
sleep 10
echo recording launch of TEST with output file qmaster:/dev/null
/home/obs/Compress/record_launch.py -d TEST -i ${testfile} qmaster:/dev/null
sleep 20
