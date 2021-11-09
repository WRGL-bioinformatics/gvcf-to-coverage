#!/bin/bash -e

# Iridis5 has no bedtools module (yet - we could request it)
# So instead just use an alias to the manually installed version.
#module load bedtools/2.21.0
# DEV: This setting is needed to allow alias in scripts
shopt -s expand_aliases
alias bedtools="/scratch/bs5n14/software/bedtools/2.30.0/bedtools"

# TODO: Will need to add a Python 3 module load here once updated. Currently works
#       as the default version is 2.7, but no default Python 3 install is available.
#       Currently the scripts are explicitly calling python2, both at the command
#       and the shebang within the scripts (could remove python command if set up
#       correctly, and just treat them like a standard script or executable).

# usage prompt for user help
usage(){
	echo 'USAGE: gvcf_to_summary.sh <ROI file> <genome file> <gVCF file>'
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

bed="$1"
genome="$2"
vcf="$3"

# get the path where the script is located.
script_path="$( dirname "$(readlink -f "$0")" )"/python_scripts

# Check that the files exist
if [ ! -f "$bed" ]; then
	>&2 echo 'ERROR: ROI file '"$bed"' does not exist or could not be opened'
	usage
	exit 1
fi

if [ ! -f "$genome" ]; then
    >&2 echo 'ERROR: Genome file '"$genome"' does not exist of could not be opened'
    usage
    exit 1
fi

if [ ! -f "$vcf" ]; then
	>&2 echo 'ERROR: Coverage file '"$vcf"' does not exist or could not be opened'
	usage
	exit 1
fi

# output filenames are derived from VCF file name
out="$( basename $vcf .gvcf.gz )"
out="$( basename $out .vcf )".summary.txt
depthout="$( basename $vcf .gvcf.gz )"
depthout="$( basename $depthout .vcf )".coverage.txt
# quickly open the VCF file to get the sample ID
sampleid=$( grep -v "^##" "$vcf" | grep "^#" | cut -f 10 )

echo "## Coverage report for: $sampleid" > "$out"
echo "## Minimum depth of coverage: $mindp" >> "$out"
echo -e "#GENE\tLENGTH\tCOVERED\tCOVERAGE" >> "$out"

# Convert the BP resolution gVCF to a whole-gene summary
if (file "$vcf" | grep -q compressed ) ; then
    zcat "$vcf" | \
    python2 "$script_path"/gvcf_to_bed.py "$mindp" | \
    bedtools map -g "$genome" -a "$bed" -b stdin -c 5 -o count -null 0 | \
    python2 "$script_path"/gene_summariser.py >> "$out"
else
    cat "$vcf" | \
    python2 "$script_path"/gvcf_to_bed.py "$mindp" | \
    bedtools map -g "$genome" -a "$bed" -b stdin -c 5 -o count -null 0 | \
    python2 "$script_path"/gene_summariser.py >> "$out"
fi

# Also create a samtools depth style report for pipeline dowload
# Set minimum depth to 0 to ensure depth is output for all target positions
if (file "$vcf" | grep -q compressed ) ; then
    zcat "$vcf" | \
    python2 "$script_path"/gvcf_to_bed.py 0 | \
    bedtools intersect -u -sorted -g "$genome" -a stdin -b "$bed" | \
    cut -f 1,3,5 > "$depthout"
else
    cat "$vcf" | \
    python2 "$script_path"/gvcf_to_bed.py 0 | \
    bedtools intersect -u -sorted -g "$genome" -a stdin -b "$bed" | \
    cut -f 1,3,5 > "$depthout"
fi
