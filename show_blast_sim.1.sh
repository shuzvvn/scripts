#!/bin/bash
read -p "sim file: " sim_file
read -p "locus_tag prefix: " locus_tag_prefix
while [ 1 ]
do
echo -e "\n\n\n\n\n\n"
read -p "Locus_tag#: " number
grep "${locus_tag_prefix}${number}" ${sim_file}
done
