#!/bin/bash

SOURCE_DIR=$1
TARGET_DIR=$2
HG19_DIR="/home/cilga/Desktop/VPC/ref/bowtie_indexes/hg19"

echo 'Going to $SOURCE_DIR'
cd $SOURCE_DIR

FILE_ARRAY=/$(ls *.fastq.gz)
partial_file_found=false

for file in $FILE_ARRAY; do
	if [ "$partial_file_found" = true ]; then
		partial_file_found=false
		continue
	fi

	echo "File: $file"

	serial_number=${file%_pass*}
	echo "Working on $serial_number"
	mkdir -p $TARGET_DIR/$serial_number
	if [[ $file = *"_pass_"* ]]; then
		echo "Partial file ..."
		partial_file_found=true

		part_one="$serial_number""_pass_1.fastq.gz"
		part_two="$serial_number""_pass_2.fastq.gz"
		echo "Parts: $part_one & $part_two"

		tophat2 -o $TARGET_DIR/$serial_number $HG19_DIR $SOURCE_DIR/$part_one 			$SOURCE_DIR/$part_two
	else
		echo "Single file ..."
		tophat2 -o $TARGET_DIR/$serial_number $HG19_DIR $SOURCE_DIR/$file
	fi

	sleep 2
done

