#!/usr/bin/env python

import sys

def summarise_gene(INFILE=sys.stdin):
    sys.stderr.write("INFO: Summarising per-gene coverage...\n")

    CURRENT_GENE = {"name": None, "length": None, "coverage": None}

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

            # When the gene name changes, print details and reset
            if genename != CURRENT_GENE["name"]:
                # "name" is initialised to None, so don't print as there's
                # no actual information yet.
                if CURRENT_GENE["name"] is not None:
                    #assert CURRENT_GENE["coverage"] <= CURRENT_GENE["length"], "ERROR: Total coverage appears greater than gene length"
                    print "%s\t%d\t%d\t%.2f" % (CURRENT_GENE["name"], CURRENT_GENE["length"], CURRENT_GENE["coverage"], (float(CURRENT_GENE["coverage"])/CURRENT_GENE["length"])*100)
                CURRENT_GENE["name"] = genename
                CURRENT_GENE["length"] = 0
                CURRENT_GENE["coverage"] = 0
            # Update CURRENT_GENE with current amplion details
            CURRENT_GENE["length"] += ampliconlength
            CURRENT_GENE["coverage"] += ampliconcoverage

    # print the last gene
#    assert CURRENT_GENE["coverage"] <= CURRENT_GENE["length"], "ERROR: Total coverage appears greater than gene length"
    print "%s\t%d\t%d\t%.2f" % (CURRENT_GENE["name"], CURRENT_GENE["length"], CURRENT_GENE["coverage"], (CURRENT_GENE["coverage"]/CURRENT_GENE["length"])*100)

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
