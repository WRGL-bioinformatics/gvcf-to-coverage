#!/usr/bin/env python2

import sys

class ProcessVcf(object):
    """
    Generates BED-format output from an input gVCF file.
    For more accurate coverage (related to depth used for calling) compared to 
    current smtools depth based methods.
    """

    def __init__(self, mindp=20, infile=sys.stdin):
        # MINDP is the mininum depth threshold
        self.MINDP = mindp
        # Track progress
        self.PASSING_DP = 0
        self.TOTAL_POS = 0
        # To store the VCF column headers
        self.HEADERS = None
        self.INFILE = infile

        # run the analysis
        self.process_file()

    def process_file(self):
        """
        Opens the INFILE and reads each line. This should support either stdin
        or a file input.

        Ignores the VCF file header lines (starting "##"), and stores the
        column headers (single "#") as these are useful as dictionary keys to
        simplify access to this data.
        """

        # Print basic progress info to STDERR
        # Sample metadata recording is handled by gvcf_to_summary.sh
        sys.stderr.write("INFO: Converting gVCF to BED for coverage...\n")
        sys.stderr.write("INFO: Minimum depth = %d\n" % (self.MINDP))
        # Print an out-of-ROI BED line to ensure bedtools doesn't throw a fit if
        # there is no coverage (map expects at least one line of BED file!)
        print "1\t0\t1\t%s\t0".format("This is a dummy line. See code for details.")

        # Identify and appropriately process file header, column header, and data lines
        for line in self.INFILE:
            line = line.rstrip()
            # ignore the VCF header lines
            if line.startswith("##"):
                pass
            # get the actual column headers
            elif line.startswith("#"):
                line = line.lstrip("#")
                self.HEADERS = line.split()
            # now process the data lines
            else:
                self.process_line(line)

        # print an output summary
        if self.PASSING_DP == 0:
            sys.stderr.write("WARNING: 0 positions passing minimum depth\n")
        sys.stderr.write("INFO: %d of %d (%.2f) positions above minimum depth\n" % (self.PASSING_DP, self.TOTAL_POS, (float(self.PASSING_DP)/self.TOTAL_POS)*100))

        try:
            sys.stdin.close()
        except:
            self.INFILE.close()
        sys.stderr.close()



    def process_line(self, line):
        """
        Process a line of gVCF input to a BED format output.

        Split the line, and then zip mergeinto a dict with the column headers for easy
        access. The individual sample details column must also be split and zipped with
        the FORMAT column

        If the DP value for the line is >= MINDP it can be printed
        """

        line = dict(zip(self.HEADERS, line.split()))
        # Join the FORMAT and INFO fields into a dict
        # This is flexible, and can account for differences between lines
        # Sample data (GT, DP, etc) is under the column headed by the sample ID
        # Since we don't know this, we have to use the fact that it's the last
        # column and get the last item of the original headers list to access it.
        sampledict = dict(zip(line["FORMAT"].split(":"),
                              line[self.HEADERS[-1]].split(":")))
        self.TOTAL_POS += 1
        try:
            if int(sampledict['DP']) >= self.MINDP:
                self.PASSING_DP += 1
                # Print the BED format output
                self.print_line(line['CHROM'], line['POS'], sampledict['DP'])
        except KeyError:
            # No DP data found. Check if DP=0 is in the INFO field, otherwise fail
            sys.stderr.write("WARNING: No DP field found at position %s:%s. Assuming 0 depth. VCF INFO %s\n" % (line['CHROM'], line['POS'], line['INFO']))
            # Print the BED format output as 0 depth only if MINDP is set to report ALL base positions
            if self.MINDP == 0:
                self.print_line(line['CHROM'], line['POS'], 0)

    @staticmethod
    def print_line(chrom, pos, depth):
        """
        Print BED format output from the gVCF data.

        Need to subtract 1 from the VCF position to get the start of a single base
        bed interval.
        """

        print "%s\t%d\t%s\t.\t%s" % (chrom, int(pos)-1, pos, depth)

if __name__ == "__main__":
    """
    Handle the various possible combinations of inputs, either directly passing a file name
    or using stdin as part of a pipeline. Allow file and depth to be specified in either
    order. Also have to allow depth to be specified in a pipe.
    """
    if len(sys.argv) == 3:
        try:
            # min depth and target file user specified
            vcfreader = ProcessVcf(mindp=int(sys.argv[1]), infile=open(sys.argv[2]))
        except ValueError:
            # min depth and target file specified but reversed
            vcfreader = ProcessVcf(infile=open(sys.argv[1]), mindp=int(sys.argv[2]))
    elif len(sys.argv) == 2:
        try:
            # depth specified, target file = stdin
            vcfreader = ProcessVcf(mindp=int(sys.argv[1]))
        except ValueError:
            # target file specified, depth = default
            vcfreader = ProcessVcf(infile=open(sys.argv[1]))
    else:
        try:
            # depth = default, target file = stdin
            vcfreader = ProcessVcf()
        except:
            sys.stderr.write("ERROR: There was a problem creating the ProcessVcf reader\n")
