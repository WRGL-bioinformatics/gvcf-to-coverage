# Run the GVCF to BED Python script only

# get the path where the script is located.
scriptpath="$( dirname "$(readlink -f "$0")" )"/python_scripts

# Iridis5 has no bedtools module (yet - we could request it)
# So instead just use an alias to the manually installed version.
#module load bedtools/2.21.0
# DEV: This setting is needed to allow alias in scripts
shopt -s expand_aliases
alias bedtools="/scratch/bs5n14/software/bedtools/2.30.0/bedtools"

usage(){
    >&2 echo "USAGE:"
    >&2 echo "      gvcf-to-bed <gVCF> <BED> <opt: min DP (default=20)>"
    exit 1
}

gvcf="$1"
mindp="${2:-20}"

if [ ! -f "$gvcf" ]; then
    >&2 echo "ERROR: Cannot open gVCF file $gvcf"
    usage
fi

"$scriptpath"/gvcf_to_bed.py "$gvcf" "$mindp" | bedtools merge
