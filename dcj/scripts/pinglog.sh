#! /bin/bash
echo date >> pinglog.txt
ping -c 10 google.com >>pinglog.txt
echo "\n" >> pinglog.txt