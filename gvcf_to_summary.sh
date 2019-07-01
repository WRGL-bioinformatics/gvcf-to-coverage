#!/usr/bin/env bash

# get the VCF file and BED file from the arguments for now
# This line could easily be added to analysis scripts instead and use those params
# Or just load the .config and .variables files here anyway?

# 20x minimum depth is required
# DEV: maybe spin out into a config file for clarity??
mindp=20

# panel bed and coverage file are passed via args, so they can be easily changed.
bed=$1
vcf=$2
# output filenames are derived from VCF file name
out="$( basename $vcf .vcf )".summary.txt
depthout="$( basename $vcf .vcf )".coverage.txt
# quickly open the VCF file to get the sample ID
sampleid=$( grep -v "^##" "$vcf" | grep "^#" | cut -f 10 )

# Check inputs
# See if there's a help flag
if [ "$panel" = "-h" -o "$panel" = "--help" ]; then
	usage
	exit 1
fi

# check that the files exist
if [ ! -f "$bed" ]; then
	echo 'ERROR: ROI file '"$bed"' does not exist or could not be opened'
	usage
	exit 1
fi
if [ ! -f "$vcf" ]; then
	echo 'ERROR: Coverage file '"$vcf"' does not exist or could not be opened'
	usage
	exit 1
fi


echo "## Coverage report for: $sampleid" > "$out"
echo "## Minimum depth of coverage: $mindp" >> "$out"
echo -e "#GENE\tLENGTH\tCOVERED\tCOVERAGE" >> "$out"

# Convert the BP resolution gVCF to a whole-gene summary
cat "$vcf" | ./gvcf_to_bed.py "$mindp" | bedtools map -a "$bed" -b stdin -o count | ./gene_summariser.py >> "$out"

# Also create a samtools depth style report for pipeline dowload
# Set minimum depth to 0 to ensure depth is output for all target positions
cat "$vcf" | ./gvcf_to_bed.py 0 | bedtools intersect -a stdin -b "$bed" | cut -f 1,3,5 > "$depthout"

usage(){
	echo 'USAGE: gvcf_to_summary.sh <ROI file> <Coverage.txt file>'
}