#!/usr/bin/env python3

import sys

print("INFO: Summarising per-gene coverage...", file=sys.stderr)

CURRENT_GENE = {"name": None, "length": None, "coverage": None}
f = sys.stdin
for line in f:
    line = line.rstrip()
    # Ignore any VCF header lines
    if line.startswith("#"):
        pass
    # Now process the data lines
    else:
        line = line.split()

        # DEV
        # Do a simple ordering check?
        # e.g. if the same chrom, is current start > previous start
        # could do be creating an amplicon object/dict and retaining
        # the previous amplicon for comparison

        # Length is just end - start
        ampliconlength = int(line[2]) - int(line[1])
        ampliconcoverage = int(line[4])
        # By WRGL convention, gene name and amplicon/exon number
        # are separated by an underscore. So we can just split on
        # this and take the first part.
        genename = line[3].split("_")[0]

        # When the gene name changes, print details and reset
        if genename != CURRENT_GENE["name"]:
            # "name" is initialised to None, so don't print as there's
            # no actual information yet.
            if CURRENT_GENE["name"] is not None:
                assert CURRENT_GENE["coverage"] <= CURRENT_GENE["length"], "ERROR: Total coverage appears greater than gene length"
                print("{}\t{}\t{}\t{:.2f}".format(CURRENT_GENE["name"], CURRENT_GENE["length"], CURRENT_GENE["coverage"], (CURRENT_GENE["coverage"]/CURRENT_GENE["length"])*100))
            CURRENT_GENE["name"] = genename
            CURRENT_GENE["length"] = 0
            CURRENT_GENE["coverage"] = 0
        # Update CURRENT_GENE with current amplion details
        CURRENT_GENE["length"] += ampliconlength
        CURRENT_GENE["coverage"] += ampliconcoverage

# print the last gene
assert CURRENT_GENE["coverage"] <= CURRENT_GENE["length"], "ERROR: Total coverage appears greater than gene length"
print("{}\t{}\t{}\t{:.2f}".format(CURRENT_GENE["name"], CURRENT_GENE["length"], CURRENT_GENE["coverage"], (CURRENT_GENE["coverage"]/CURRENT_GENE["length"])*100))

