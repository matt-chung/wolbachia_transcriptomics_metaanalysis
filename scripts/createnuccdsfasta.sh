#!/bin/bash

TEMP=`getopt -o n:g:f:i: --long nuc-fasta:,gff:,feature-id:,id-attribute: -n 'createnuccdsfasta.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -n|--nuc-fasta)
            case "$2" in
                "") shift 2 ;;
                *) nuc_fasta=$2 ; shift 2 ;;
            esac ;;
        -g|--gff)
            case "$2" in
                "") shift 2 ;;
                *) gff=$2 ; shift 2 ;;
            esac ;;
        -f|--feature-id)
            case "$2" in
                "") shift 2 ;;
                *) feature_id=$2 ; shift 2 ;;
            esac ;;
        -i|--id-attribute)
	    case "$2" in
	        "") shift 2 ;;
                *) id_attribute=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error" ; exit 1 ;;
    esac
done

samtools_bin_dir="/usr/local/packages/samtools-1.9/bin"

#nuc_fasta="/local/projects-t3/EBMAL/mchung_dir/fadu/geark/Arkansas.fsa"
#gff="/local/projects-t3/EBMAL/mchung_dir/fadu/geark/Arkansas.gff3"
#feature_id="CDS"
#id_attribute="ID"
#output_file=./temp

IFS=""

"$samtools_bin_dir"/samtools faidx "$nuc_fasta"

while read Line
do

if [ "$(echo $Line | awk -F "\t" '{print $3}')" = "$feature_id" ]
then
contig=$(echo $Line | awk -F "\t" '{print $1}')
start=$(echo $Line | awk -F "\t" '{print $4}')
stop=$(echo $Line | awk -F "\t" '{print $5}')
strand=$(echo $Line | awk -F "\t" '{print $7}')
id=$(echo $Line | awk -F "\t" '{print $9}' | sed "s/.*"$id_attribute"=//g" | sed "s/;.*//g")

if [ "$strand" = "+" ]
then
"$samtools_bin_dir"/samtools faidx "$nuc_fasta" "$contig":"$start"-"$stop" | sed -e "s/>.*/>"$id"/g"
else
"$samtools_bin_dir"/samtools faidx -i "$nuc_fasta" "$contig":"$start"-"$stop" | sed -e "s/>.*/>"$id"/g"
fi
fi
done < "$gff"
