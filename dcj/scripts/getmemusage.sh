#!/bin/sh
#GetmemUsageModified

USAGE="Usage&#58; $0 processName"

if &#91; $# -ne 1 &#93;; then
   echo $USAGE
   exit 1
fi

# In case the monitored process has not yet started
# keep searching until its PID is found
PROCESS_PID=""
while &#58;
do
   PROCESS_PID=`/sbin/pidof $1`

   if &#91; "$PROCESS_PID.X" != ".X" &#93;; then
      break
   fi
done

LOG_FILE="memusage.csv"

echo "ElapsedTime,VmSize,VmRSS" > $LOG_FILE

ELAPSED_TIME=`date`
PERIOD=2        # seconds

# Monitor memory usage forever until script is killed
while &#58;
do
   VM_SIZE=`awk '/VmSize/ &#123;print $2&#125;' < /proc/$PROCESS_PID/status`
   if &#91; "$VM_SIZE.X" = ".X" &#93;; then
      continue
   fi