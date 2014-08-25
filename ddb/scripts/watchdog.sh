#!/bin/sh
echo 'Watching' >>watch.txt
date -u >> watch.txt
ping -c 1 morph.berkeley.edu >> watch.txt
scp -i .ssh/id_watch_rsa watch.txt ddeboer@morph.berkeley.edu:.
