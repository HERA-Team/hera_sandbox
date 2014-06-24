#!/bin/bash
#
RETVAL=0;

start() {
echo “Starting still_task_server”
}

stop() {
echo “Stopping still_task_server”
}

restart() {
stop
start
}

case “$1″ in
start)
  start-stop-daemon --start  -m --pidfile=/var/run/still_task_server.pid -n task_server -b task_server.py  
;;
stop)
  stop-stop-daemon --stop --pidfile=/var/run/still_task_server.pid -n task_server -b task_server.py
;;
restart)
  restart
;;
*)

echo $”Usage: $0 {start|stop|restart}”
