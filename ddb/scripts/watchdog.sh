#!/bin/sh
echo 'Watching' >>watch.txt
date -u >> watch.txt
ping -c 1 enterprise.sese.asu.edu >> watch.txt
scp -i .ssh/id_watch_rsa watch.txt ddeboer@enterprise.sese.asu.edu:.
time -p -o watchBand.dat scp -i .ssh/id_watch_rsa watch150M.dat ddeboer@enterprise.sese.asu.edu:.
scp -i .ssh/id_watch_rsa watchBand.dat ddeboer@enterprise.sese.asu.edu:.
