GVCF COVERAGE SCRIPTS
=====================

Author: ben.sanders@nhs.net

In some circumstances, the current Samtools depth coverage output is inadequate - it
does not accurately represent the depth of coverage used for variant calling. In at
least one case this has lead to a missed variant - depth for calling was <20x, and so
would have triggered a Sanger infill that would have detected the variant, however 
Samtools depth suggested the coverage was >70x.

These scripts produce a per-gene summary, and also mimic the Samtools depth output used
by the MiSeq pipeline, using the GVCF produced by HaplotypeCaller. The depth information 
in this file shows the depth of reads actually considered by the variant caller, and so 
is a more accurate reflection of the true coverage.

Usage
-----

The scripts are written into the anlaysis scripts, but are simple to call on thier own
as part of a basic pipeline:

per-gene coverage:
cat <vcf> | \
./gvcf_to_bed.py <mindp DEFAULT=20> | \
bedtools map -a <roi bed> -b stdin -o count | \
./gene_summariser.py 

samtools depth output:
cat <vcf> | \
./gvcf_to_bed.py 0 | \
bedtools intersect -a stdin -b <roi bed> | \
cut -f 1,3,5