#! /bin/bash

for file in $* ; do
    rm -rf `python -c "import os; print os.path.basename('$file')"`
    scp -r -c arcfour256 $file .
done
