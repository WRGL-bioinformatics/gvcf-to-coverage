#!/usr/bin/env python

import sys
import collections

def summarise_gene(INFILE=sys.stdin):
    sys.stderr.write("INFO: Summarising per-gene coverage...\n")

    # Use an OrderedDict to keep the genes in chromsomal order
    genedict = collections.OrderedDict()

    for line in INFILE:
        line = line.rstrip()
        # Ignore any VCF header lines
        if line.startswith("#"):
            pass
        # Now process the data lines
        else:
            line = line.split()
            # Length is just end - start
            ampliconlength = int(line[2]) - int(line[1])
            ampliconcoverage = int(line[4])
            # By WRGL convention, gene name and amplicon/exon number
            # are separated by an underscore. So we can just split on
            # this and take the first part.
            genename = line[3].split("_")[0]
            # add each region to the appropriate gene in the dict
            try:
                genedict[genename]['covered'] += ampliconcoverage
                genedict[genename]['length'] += ampliconlength
            except KeyError:
                genedict[genename] = {'covered': ampliconcoverage, 'length': ampliconlength}

    # print the gene summaries, calculating %age coverage
    for k,v in genedict.items():
        print("%s\t%d\t%d\t%.2f" % (k, v['length'], v['covered'], (float(v['covered'])/v['length'])*100))

    # close file streams
    try:
        sys.stdin.close()
    except:
        INFILE.close()
    sys.stderr.close()


if __name__ == "__main__":
    try:
        summarise_gene(INFILE=sys.argv[1])
    except IndexError:
        summarise_gene()
