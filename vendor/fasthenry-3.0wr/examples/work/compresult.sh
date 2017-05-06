#!/bin/bash
for f in *.mat
do
	echo $f
	cmp $f ./results/$f
done
