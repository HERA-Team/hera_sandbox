#! /bin/bash

for file in $* ; do
    rm -rf `import os; print os.path.basename('$file')`
    scp -r -c arcfour256 $file .
done
