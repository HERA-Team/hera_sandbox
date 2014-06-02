#! /bin/bash

rm -rf $1
echo scp -r -c arcfour256 $2 .
scp -r -c arcfour256 $2 .
