#!/usr/bin/env python3

import sys

# MINDP is the mininum depth threshold
# Default to 20, but allow to be overriden from the cmd line
try:
    MINDP = int(sys.argv[1])
except IndexError:
    MINDP = 20

print("INFO: Converting gVCF to BED for coverage...", file=sys.stderr)
print("INFO: Minimum depth = {}".format(MINDP), file=sys.stderr)

# Print an out-of-ROI BED line to ensure bedtools doesn't throw a fit if
# there is no coverage (map expects at least one line of BED file!)
print("1\t0\t1\t{}\t0".format("This is a dummy line. See code for details."), file=sys.stdout)

# Track progress
PASSING_DP = 0
TOTAL_POS = 0

f = sys.stdin
for line in f:
    line = line.rstrip()
    # ignore the VCF header lines
    if line.startswith("##"):
        pass
    # get the actual column headers
    elif line.startswith("#"):
        line = line.lstrip("#")
        headers = line.split()
    # now process the data lines
    else:
        line = line.split()
        line = dict(zip(headers, line))
        # Join the FORMAT and INFO fields into a dict
        # This is flexible, and can account for differences between lines
        sampleheaders = line["FORMAT"].split(":")
        # Sample data (GT, DP, etc) is under the column headed by the sample ID
        # Since we don't know this, we have to use the fact that it's the last
        # column and get the last item of the original headers list to access it.
        sampledata = line[headers[-1]].split(":")
        sampledict = dict(zip(sampleheaders, sampledata))
        TOTAL_POS += 1
        if int(sampledict['DP']) >= MINDP:
            PASSING_DP += 1
            # Output in BED format (subtract 1 from POS for the correct indexed start position)
            print("{}\t{}\t{}\t.\t{}".format(line['CHROM'],
                                             int(line['POS'])-1,
                                             line['POS'],
                                             sampledict['DP']), file=sys.stdout)

if PASSING_DP == 0:
    print("WARNING: 0 positions passing minimum depth", file=sys.stderr)
print("INFO: {} of {} ({:.2f}%) positions above minimum depth".format(PASSING_DP, TOTAL_POS, (PASSING_DP/TOTAL_POS)*100), file=sys.stderr)
