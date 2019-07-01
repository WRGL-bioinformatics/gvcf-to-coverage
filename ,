#!/usr/bin/env bash

# get the VCF file and BED file from the arguments for now
# This line could easily be added to analysis scripts instead and use those params
# Or just load the .config and .variables files here anyway?

vcf="181112_M00321_0257_000000000-C5KD2_240240121.vcf"
bed="PanelROI_WRGL_52893_b37_v1.bed"
out="$( basename $vcf .vcf )".summary.txt
depthout="$( basename $vcf .vcf )".coverage.txt
mindp=20

sampleid=$( grep -v "^##" "$vcf" | grep "^#" | cut -f 10 )

echo "## Coverage report for: $sampleid" > "$out"
echo "## Minimum depth of coverage: $mindp" >> "$out"
echo -e "#GENE\tLENGTH\tCOVERED\tCOVERAGE" >> "$out"

# Convert the BP resolution gVCF to a whole-gene summary
cat "$vcf" | ./gvcf_to_bed.py "$mindp" | bedtools map -a "$bed" -b stdin -o count | ./gene_summariser.py >> "$out"

# Also create a samtools depth style report for pipeline dowload
# Set minimum depth to 0 to ensure depth is output for all target positions
cat "$vcf" | ./gvcf_to_bed.py 0 | bedtools intersect -a stdin -b "$bed" | cut -f 1,3,5 > "$depthout"

