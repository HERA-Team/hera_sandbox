#!/bin/sh
echo 'Watching' >>watch.dat
date -u >> watch.dat
ping -c 1 morph.berkeley.edu >> watch.dat
scp -i .ssh/id_watch_rsa watch.dat ddeboer@morph.berkeley.edu:.
