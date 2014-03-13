#!/bin/sh
list=`fgrep -l params.h *.F`
for file in $list
do
	ofile=${PVM_ARCH}/${file%.F}.o
	if [ ! -s $ofile ]; then
		echo "$0: error: cannot find $ofile"
	else
		rm $ofile
	fi
done
