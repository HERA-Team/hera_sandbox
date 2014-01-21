#!/bin/bash

echo This is a test script for catching things

awk 'BEGIN{
    FS=".";
    first=0;
}

{   
    if ( first == 0) {
        print "Working in directory " $1;
    }
    gap = $3 - first;
    if (((gap != 696) && (gap != 695))  && first!=0) {
        print $0 "is a bad file";
        print  "Gap= " $3-first;
        first=0;
    }
    else{
    first=$3;
    }
}'
echo finished awk. Quitting....

