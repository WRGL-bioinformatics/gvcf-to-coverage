# Run the GVCF to BED Python script only

# get the path where the script is located.
scriptpath="$( dirname "$(readlink -f "$0")" )"/python_scripts

# check if bedtools module is already loaded (e.g. from pipeline)
bedtools --version > /dev/null 2> /dev/null
if [[  $? -ne 0 ]]; then
    # Non-zero exit means bedtools isn't loadedy
    module load bedtools
fi

usage(){
    >&2 echo "USAGE:"
    >&2 echo "      gvcf-to-bed <gVCF> <minDP (optional. default=20)"
    exit 1
}

gvcf="$1"
mindp="${2:-20}"

if [ ! -f "$gvcf" ]; then
    >&2 echo "ERROR: Cannot open gVCF file $gvcf"
    usage
fi

"$scriptpath"/gvcf_to_bed.py "$gvcf" "$mindp" | bedtools merge
