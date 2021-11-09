# Run the GVCF to BED Python script only

# get the path where the script is located.
script_path="$( dirname "$(readlink -f "$0")" )"/python_scripts

usage(){
    >&2 echo "USAGE:"
    >&2 echo "      gvcf-to-bed <gVCF> <BED> <opt: min DP (default=20)>"
    exit 1
}

gvcf="$1"
bed="$2"
mindp="${2:-20}"

if [ ! -f "$gvcf" ]; then
    >&2 echo "ERROR: Cannot open gVCF file $gvcf"
    usage
fi

if [ ! -f "$bed" ]; then
    >&2 echo "ERROR: Cannot open BED file $bed"
    usage
fi

"$scriptpath"/gvcf_to_bed.py

