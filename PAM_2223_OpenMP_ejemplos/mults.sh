#! /bin/bash


for n in 512 1024 2048 4096;
do
	echo 'n=' $n
	for t in 1 2 4 8 16;
	do
		echo 't=' $t
		for f in 1 2 4 8;
		do
			echo 'f=' $f
            ./mulmatcua $n $f $t
		done
	done

done
