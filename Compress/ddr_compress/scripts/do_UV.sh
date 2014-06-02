#! /bin/bash

rm -rf $1
scp -r -c arcfour256 $2 .
