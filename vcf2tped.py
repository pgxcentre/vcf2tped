#!/usr/bin/env python

import os
import re
import sys
import gzip
import argparse
from shutil import copyfile
from collections import defaultdict


# The version of the script
prog_version = 0.1

def main():
    """The main function.

    """
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    # Summarize the options used
    print "# Command used:"
    print "{} \\".format(sys.argv[0])
    print "    --vcf {} \\".format(args.vcf)
    print "    --out {} \\".format(args.out)

    # Converting
    convert_vcf(args.vcf, args.out)


def convert_vcf(i_filename, o_prefix):
    """Reads a VCF and creates a TPED/TFAM from it.

    :param i_filename: the name of the VCF file (might be gzip).
    :param o_prefix: the prefix of the output files.

    :type i_filename: string
    :type o_prefix: string

    """
    # The open function
    open_f = open
    if re.search("\.gz$", i_filename):
        open_f = gzip.open

    # Reading the file
    with open_f(i_filename, 'rb') as i_file:
        line = i_file.readline()

        # We read until the header
        while re.search("^##", line) is not None:
            line = i_file.readline()

        # This should be the header
        if re.search("^##CHROM", line) is not None:
            msg = "{}: no header".format(i_filename)
            raise ProgramError(msg)

        # Creating the header
        row = line.rstrip("\r\n").split("\t")
        header = {name:  i for i, name in enumerate(row)}

        # Checking some names
        for name in ["#CHROM", "POS", "ID", "REF", "ALT", "FORMAT"]:
            if name not in header:
                msg = "{}: no column named {}".format(i_filename, name)
                raise ProgramError(msg)

        # Printing the TFAMs
        tfam_names = ["{}.snv.2_alleles.tfam".format(o_prefix),
                      "{}.snv.n_alleles.tfam".format(o_prefix),
                      "{}.indel.2_alleles.tfam".format(o_prefix),
                      "{}.indel.n_alleles.tfam".format(o_prefix)]
        print_same_tfams(tfam_names, row[header["FORMAT"]+1:])

        # Those positions have been seen
        seen_pos = defaultdict(int)

        # The output files
        tped_snv_2 = None
        tped_snv_n = None
        tped_indel_2 = None
        tped_indel_n = None
        snv_ref = None
        try:
            tped_snv_2 = open("{}.snv.2_alleles.tped".format(o_prefix), 'w')
            tped_snv_n = open("{}.snv.n_alleles.tped".format(o_prefix), 'w')
            tped_indel_2 = open("{}.indel.2_alleles.tped".format(o_prefix), "w")
            tped_indel_n = open("{}.indel.n_alleles.tped".format(o_prefix), "w")
            snv_ref = open("{}.snv.ref".format(o_prefix), "w")
            indel_ref = open("{}.indel.ref".format(o_prefix), "w")
        except IOError:
            msg = "couldn't write output files"
            raise ProgramError(msg)

        # Reading the rest of the data
        for line in i_file:
            row = line.rstrip("\r\n").split("\t")

            # Getting the information
            chrom = encode_chr(row[header["#CHROM"]])
            pos = row[header["POS"]]
            name = row[header["ID"]]
            alleles = [row[header["REF"]]] + row[header["ALT"]].split(",")
            g_format = row[header["FORMAT"]].split(":")
            g_format = {name: i for i, name in enumerate(g_format)}
            genotypes = [i.split(":")[g_format["GT"]]
                            for i in row[header["FORMAT"]+1:]]

            # Getting rid of the "." (so that it becomes "./.")
            genotypes = [re.sub("^\.$", "./.", i) for i in genotypes]

            # Is this an INDEL?
            indel = False
            for allele in alleles:
                if len(allele) > 1:
                    indel = True
                    break

            # The output file to choose from (default is SNV and bi-allelic)
            o_file = tped_snv_2
            o_ref_file = snv_ref
            if len(alleles) > 2:
                o_file = tped_snv_n

            # Constructing the genotype map
            g_map = {str(i): a for i, a in enumerate(alleles)}
            if indel:
                # This is a new genotype map, only for INDELs
                g_map = {str(i): str(i+1) for i in xrange(len(alleles))}

                # These are the required files
                o_ref_file = indel_ref
                o_file = tped_indel_2
                if len(alleles) > 2:
                    o_file = tped_indel_n

            # Adding the unknown genotype
            g_map["."] = "0"

            # Checking if the position have a name
            if name == ".":
                name = "{}:{}".format(chrom, pos)

            # We increase the number of time we saw this name, and check if we
            # saw it more than once
            seen_pos[name] += 1
            if seen_pos[name] > 1:
                name = "{}-{}".format(name, seen_pos[name])

            # The first part of the TPED line
            first_part = [chrom, name, "0", pos]

            genotypes = [i.split("/") for i in genotypes]
            genotypes = [recode_genotype(i, g_map, chrom, pos)
                                                for i in genotypes]
            print >>o_file, "\t".join(first_part + genotypes)

            # Saving the alleles
            print >>o_ref_file, "\t".join([chrom, pos, name, alleles[0],
                                           ",".join(alleles[1:])])

        # Closing the output files
        tped_snv_2.close()
        tped_snv_n.close()
        tped_indel_2.close()
        tped_indel_n.close()
        snv_ref.close()
        indel_ref.close()


def recode_genotype(genotype, encoding, chromosome, position):
    """Encodes the genotypes.

    :param genotype: the genotypes (list of alleles)
    :param encoding: the allele encoding
    :param chromosome: the chromosome on which the marker is
    :param position: the position of the marker

    :type genotype: list of str
    :type encoding: map
    :type chromosome: str
    :type position: str

    :returns: a string with the two alleles separated by a space.

    """
    if len(genotype) == 1:
        # This is an haploid genotype
        if chromosome == "23" or chromosome == "24":
            return " ".join(genotype * 2)
        else:
            return " ".join([encoding["."]] * 2)

    return "{} {}".format(encoding[genotype[0]], encoding[genotype[1]])


def print_same_tfams(file_names, sample_names):
    """Print the same TFAM in different files.

    :param file_names: list of output file names
    :param sample_names: the list of samples

    :type file_names: list of string
    :type sample_names: list of string

    """
    try:
        # Printing only the first TFAM
        with open(file_names[0], 'w') as o_file:
            for name in sample_names:
                print >>o_file, "\t".join([name] * 2 + ["0"] * 3 + ["-9"])

        # Copying the first TFAM into the others
        for name in file_names[1:]:
            copyfile(file_names[0], name)
    except IOError:
        msg = "couldn't write output TFAMs"
        raise ProgramError(msg)


def encode_chr(chrom):
    """Encodes a chromosome.

    :param chrom: the chromosome to encode
    :type chrom: string

    Remove the leading "chr" if necessary, and encodes sex chromosomes and other
    to their numeric counterpart.

        - X    X chromosome                    -> 23
        - Y    Y chromosome                    -> 24
        - XY   Pseudo-autosomal region of X    -> 25
        - MT   Mitochondrial                   -> 26

    .. doctest::

        >>> [encode_chr(str(i)) for i in range(0, 11)]
        ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
        >>> [encode_chr("chr{}".format((i))) for i in range(11, 21)]
        ["11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
        >>> [encode_chr(str(i)) for i in range(21, 27)]
        ["21", "22", "23", "24", "25", "26"]
        >>> [encode_chr(i) for i in ["X", "Y", "XY", "MT"]]
        ["23", "24", "25", "26"]
        >>> [encode_chr("chr{}".format(i) for i in range(21, 27)]
        ["23", "24", "25", "26"]
        >>> encode_chr("27")
        Traceback (most recent call last):
            ...
        ProgramError: 27: not a valid chromosome
        >>> encode_chr("0")
        Traceback (most recent call last):
            ...
        ProgramError: 0: not a valid chromosome
        >>> encode_chr("chrXX")
        Traceback (most recent call last):
            ...
        ProgramError: chrXX: not a valid chromosome

    """
    # Removing the leading "chr"
    chrom = re.sub("^chr", "", chrom, flags=re.I)

    # Checking for special chromosome
    if re.match("^x$", chrom, flags=re.I):
        return "23"
    if re.match("^y$", chrom, flags=re.I):
        return "24"
    if re.match("^xy$", chrom, flags=re.I):
        return "25"
    if re.match("^mt$", chrom, flags=re.I):
        return "26"

    try:
        numerical_chrom = int(chrom)
        if numerical_chrom >= 1 and numerical_chrom <= 26:
            return chrom
    except ValueError:
        pass

    msg = "{}: not a valid chromosome".format(chrom)
    raise ProgramError(msg)


def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    if not os.path.isfile(args.vcf):
        msg = "{}: no such file".format(args.vcf)
        raise ProgramError(msg)

    if not (re.search("\.vcf$", args.vcf) or re.search("\.vcf\.gz$", args.vcf)):
        msg = "{}: not a vcf file".format(args.vcf)
        raise ProgramError(msg)

    return True


def parse_args():
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the different
              options.

    ===============   =======  ================================================
        Options        Type                      Description
    ===============   =======  ================================================
    ===============   =======  ================================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    return parser.parse_args()


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.
    
    :param msg: the message to print to the user before exiting.
    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


# The parser object
desc = """Convert a VCF to a TPED (version {}).""".format(prog_version)
parser = argparse.ArgumentParser(description=desc)

# The input file
group = parser.add_argument_group("Input Files")
group.add_argument("--vcf", type=str, metavar="FILE", required=True,
                   help="The VCF file.")

# The output file
group = parser.add_argument_group("Output Files")
group.add_argument("-o", "--out", type=str, metavar="STR",
                   default="vcf_to_tped",
                   help=("The suffix of the output file. "
                         "[Default: %(default)s]"))


# Calling the main, if necessery
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
