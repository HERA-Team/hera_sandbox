#! /bin/bash

ARGS=`pull_args.py $*`
declare -a POLS=( "xx" "yy" "xy" "yx" )

for FILE in "${ARGS}"; do
	echo $FILE
	for POL in "${POLS[@]}"; do
		echo $POL
		echo pull_antpols.py $FILE -p $POL
		pull_antpols.py $FILE -p $POL
		echo mv ${FILE}A ${FILE:0:18}${POL}${FILE:17}
		mv ${FILE}A ${FILE:0:18}${POL}${FILE:17}
	done
done
