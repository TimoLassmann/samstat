#!/bin/bash


declare -a files=("aux.fa" "c1.fa" "ce.fa" "xx.fa" "aux#aux.sam" "c1#pad3.sam" "ce#large_seq.sam" "ce#unmap2.sam" "xx#minimal.sam" "c1#bounds.sam" "ce#1.sam" "ce#tag_depadded.sam" "fieldarith.sam" "xx#pair.sam" "c1#clip.sam" "ce#2.sam" "ce#tag_padded.sam" "xx#blank.sam" "xx#rg.sam" "c1#pad1.sam" "ce#5.sam" "ce#unmap.sam" "xx#large_aux.sam" "xx#triplet.sam" "c1#pad2.sam" "ce#5b.sam" "ce#unmap1.sam" "xx#large_aux2.sam" "xx#unsorted.sam")

echo "Running samstat tests:";


for file in "${files[@]}"
do
	error=$( ../src/samstat -l ${testdatafiledir}/$file   2>&1 )
	status=$?
	if [[ $status -eq 0 ]]; then
		printf "%10s%20s%10s\n"  samstat $file SUCCESS;
	else
		printf "%10s%20s%10s\n"  samstat $file FAILED;
	printf "with ERROR $status and Message:\n\n$error\n\n";
	exit 1;
fi

done




