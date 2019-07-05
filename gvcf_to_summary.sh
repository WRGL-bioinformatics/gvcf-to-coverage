#!/usr/bin/env bash

module load bedtools/2.21.0
module load python/3.5.1

# usage prompt for user help
usage(){
	echo 'USAGE: gvcf_to_summary.sh <ROI file> <Coverage.txt file>'
}

# Check inputs
# See if there's a help flag
if [ "$1" = "-h" -o "$1" = "--help" ]; then
	usage
	exit 1
fi


# 20x minimum depth is required
# DEV: maybe spin out into a config file for clarity??
mindp=20

bed=$1
vcf=$2

# get the path where the script is located, NOT the pwd of the user panel bed and coverage file are passed via args, so they can be easily changed.
# tried $BASH_SOURCE without success, so using $0 even though it's not recommendedbed=$1
script_path="${0%/*}"/python_scripts

# Check that the files exist
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

# output filenames are derived from VCF file name
out="$( basename $vcf .vcf )".summary.txt
depthout="$( basename $vcf .vcf )".coverage.txt
# quickly open the VCF file to get the sample ID
sampleid=$( grep -v "^##" "$vcf" | grep "^#" | cut -f 10 )

echo "## Coverage report for: $sampleid" > "$out"
echo "## Minimum depth of coverage: $mindp" >> "$out"
echo -e "#GENE\tLENGTH\tCOVERED\tCOVERAGE" >> "$out"

# Convert the BP resolution gVCF to a whole-gene summary
cat "$vcf" | \
python3 "$script_path"/gvcf_to_bed.py "$mindp" | \
bedtools map -g /scratch/WRGL/REFERENCE_FILES/REFERENCE_GENOME/GRCh37_no_gl000201.genome -a "$bed" -b stdin -c 5 -o count -null 0 | \
python3 "$script_path"/gene_summariser.py >> "$out"

# Also create a samtools depth style report for pipeline dowload
# Set minimum depth to 0 to ensure depth is output for all target positions
cat "$vcf" | \
python3 "$script_path"/gvcf_to_bed.py 0 | \
bedtools intersect -sorted -g /scratch/WRGL/REFERENCE_FILES/REFERENCE_GENOME/GRCh37_no_gl000201.genome -a stdin -b "$bed" | \
cut -f 1,3,5 > "$depthout"
