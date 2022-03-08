#! /bin/bash


for n in 128 256 512 1024 2048 4096;
do
	echo 'n=' $n
	for b in 4 8 16 32 64 128 256;
	do
		echo 'b=' $b
		for t in 1 2 3 6 9 11;
		do
			echo 'threads=' $t
			./mulmatcuablo $n $b $t
		done
	done

done


