"""VCF Comparison Tool.
Created on Thu Oct 24 11:18:19 2019
@author: medinatd
"""
import gzip
import io
import sys
import os
import argparse


class BED:
    def __init__(self, BEDfile):
        self.intervals = {}
        self.read_from_file(BEDfile)

    def read_from_file(self, BEDfile):
        with open(BEDfile, "r") as infile:
            bed_lines = infile.readlines()
        bed_lines = [line.strip().split("\t") for line in bed_lines if not line.startswith("#")]
        chrom_list = set([bed_line[0] for bed_line in bed_lines])

        self.intervals = {chrom: [] for chrom in chrom_list}

        for bed_line in bed_lines:
            self.intervals[bed_line[0]].append(range(int(bed_line[1]) + 1, int(bed_line[2]) + 1))
        return

    def sort_intervals(self):
        for chrom in self.intervals:
            self.intervals[chrom].sort(key=lambda x: x.start)
        return


class Variant:
    def __init__(self, chrom,
                 pos,
                 rsid,
                 ref,
                 alt,
                 qual,
                 filter_status,
                 info,
                 format_field,
                 genotypes):
        self.chr = chrom
        self.pos = pos
        self.id = rsid
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_status
        self.info = info
        self.format = format_field
        self.genotypes = genotypes

    def __eq__(self, other):
        """Override the default Equals behavior."""
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __repr__(self):
        to_str = [
            self.chr.decode(),
            str(self.pos),
            self.id.decode(),
            self.ref.decode(),
            ",".join([x.decode() for x in self.alt]),
            self.qual.decode(),
            ",".join([x.decode() for x in self.filter]),
            # ";".join([f"{x.decode()}={self.info[x].decode()}" for x in self.info]),
            ":".join([x.decode() for x in self.format]),
            "\t".join([":".join(x.decode() for x in y.values()) for y in self.genotypes.values()])
            ]
        return "\t".join(to_str)

    def __hash__(self):
        return hash(self.__dict__.__str__())

    def raw(self):
        infos = []
        for key, value in self.info.items():
            if isinstance(value, list):
                new_value = b"|".join(value).decode()
            else:
                new_value = value.decode()
            infos.append(f"{key.decode()}={new_value}")
        infos = ";".join(infos)

        to_str = [
            self.chr.decode(),
            str(self.pos),
            self.id.decode(),
            self.ref.decode(),
            ",".join([x.decode() for x in self.alt]),
            self.qual.decode(),
            ",".join([x.decode() for x in self.filter]),
            infos,
            ":".join([x.decode() for x in self.format]),
            "\t".join([":".join(x.decode() for x in y.values()) for y in self.genotypes.values()])
            ]

        return "\t".join(to_str)

    def calculate_variant_length(self):
        self.span = max([len(allele) for allele in self.alt + [self.ref]])
        return

    def count_alt_alleles(self):
        return len(self.alt)

    def detect_homo_ref_genos(self):
        for sample in self.genotypes:
            if self.genotypes[sample][b"GT"] == b"0/0":
                print(f"{self.genotypes[sample].decode()}: {self.genotypes[sample]['GT'].decode()}")


class VCF:
    def __init__(self, vcf_file=""):
        self.vcf_file = vcf_file
        self.header = []
        self.fields = []
        self.samples = []
        self.variants = []
        self.span_set = False

    def __hash__(self):
        return hash(self.__dict__.__str__())

    def read_header_from_file(self, vcf=None):
        if vcf is None:
            vcf = self.vcf_file
        infile = io.BufferedReader(gzip.open(vcf))
        for line in infile:
            if line.startswith(b"#"):
                self.header.append(line.strip())
            else:
                break
        infile.close()
        self.fields = self.header[-1][1:].split(b"\t")
        self.samples = self.fields[9:]
        return

    def read_body_from_file(self, vcf=None):
        if vcf is None:
            vcf = self.vcf_file
        infile = io.BufferedReader(gzip.open(vcf))
        body = infile.readlines()
        infile.close()
        body = [line.strip() for line in body if not line.startswith(b"#")]
        for line in body:
            self.variants.append(self.convert_string_variant_to_dict(line))
        return

    def read_from_file(self, vcf=None):
        if vcf is None:
            vcf = self.vcf_file
        self.read_header_from_file(vcf)
        self.read_body_from_file(vcf)
        self.convert_dict_variants_to_objs()
        return

    def convert_string_variant_to_dict(self, var_str):
        variant = var_str.split(b"\t")
        field_count = len(variant)
        variant[1] = int(variant[1])
        variant[4] = variant[4].split(b",")
        variant[6] = variant[6].split(b";")
        variant[7] = variant[7].split(b";")
        info_field = {}
        for info in variant[7]:
            if b"=" in info:
                info_split = info.split(b"=")
                info_field[info_split[0]] = b"".join(info_split[1:])
            else:
                info_field[info] = info
        if b"ANN" in info_field:
            info_field[b"ANN"] = info_field[b"ANN"].split(b"|")
        variant[7] = info_field

        if field_count < 9:
            variant.append([b"."])
            variant.append({b".": {}})
            return dict(zip(self.fields[:9] + [b"FORMAT", b"GENOTYPE"], variant))

        variant[8] = variant[8].split(b":")
        if field_count < 10:
            variant.append({b".": {}})
        else:
            variant[9] = [geno.split(b":") for geno in variant[9:]]
            variant[9] = [dict(zip(variant[8], geno)) for geno in variant[9]]
            variant[9] = dict(zip(self.samples, variant[9]))
        return dict(zip(self.fields[:9] + [b"GENOTYPE"], variant))

    def convert_dict_variants_to_objs(self):
        self.variants = [Variant(*variant.values()) for variant in self.variants]
        return

    def calculate_variant_lengths(self):
        for variant in self.variants:
            variant.calculate_variant_length()
        self.span_set = True
        return


class VCF_Comparison:
    def __init__(self, VCF1, VCF2):
        self.VCF1 = VCF1
        self.VCF2 = VCF2

        self.chrom_set = self.make_chrom_set()
        self.VCF1_chrom_dict = self.split_per_chrom(VCF1)
        self.VCF2_chrom_dict = self.split_per_chrom(VCF2)

        self.shared_pos_vars = {"VCF1": [], "VCF2": []}
        self.unshared_pos_vars = {"VCF1": [], "VCF2": []}

        self.shared_allele_vars = {"VCF1": [], "VCF2": []}
        self.unshared_allele_vars = {"VCF1": [], "VCF2": []}

        self.shared_genotype_vars = {"VCF1": [], "VCF2": []}
        self.unshared_genotype_vars = {"VCF1": [], "VCF2": []}

        self.shared_filter_vars = {"VCF1": [], "VCF2": []}
        self.unshared_filter_vars = {"VCF1": [], "VCF2": []}

    def __hash__(self):
        return hash(self.__dict__.__str__())

    def check_span_set(self, VCF):
        if VCF.span_set is True:
            return
        else:
            VCF.calculate_variant_lengths()
        return

    def make_chrom_set(self):
        chrom_set = set()
        for variant in self.VCF1.variants + self.VCF2.variants:
            chrom_set.add(variant.chr)
        return chrom_set

    def split_per_chrom(self, VCF):
        self.check_span_set(VCF)
        chrom_dict = {chrom: [] for chrom in self.chrom_set}
        for variant in VCF.variants:
            chrom_dict[variant.chr].append(variant)
        return chrom_dict

    def compare_by_pos(self):
        self.shared_pos_vars = {"VCF1": [], "VCF2": []}
        self.unshared_pos_vars = {"VCF1": [], "VCF2": []}
        for chromosome in self.chrom_set:
            pos_list1 = set([variant.pos for variant in self.VCF1_chrom_dict[chromosome]])
            pos_list2 = set([variant.pos for variant in self.VCF2_chrom_dict[chromosome]])
            shared_pos = pos_list1 & pos_list2
            for VCF_by_chrom, VCF in zip([self.VCF1_chrom_dict[chromosome], self.VCF2_chrom_dict[chromosome]],
                                         ["VCF1", "VCF2"]):
                for variant in VCF_by_chrom:
                    if variant.pos in shared_pos:
                        self.shared_pos_vars[VCF].append(variant)
                    else:
                        self.unshared_pos_vars[VCF].append(variant)
        return

    def compare_by_alleles(self):
        self.shared_allele_vars = {"VCF1": [], "VCF2": []}
        self.unshared_allele_vars = {"VCF1": [], "VCF2": []}
        for var1, var2 in zip(self.shared_pos_vars["VCF1"], self.shared_pos_vars["VCF2"]):
            if sorted(var1.alt) == sorted(var2.alt):
                self.shared_allele_vars["VCF1"].append(var1)
                self.shared_allele_vars["VCF2"].append(var2)
            else:
                self.unshared_allele_vars["VCF1"].append(var1)
                self.unshared_allele_vars["VCF2"].append(var2)
        return

    def compare_by_genotypes(self, shared_list=None):
        if shared_list is None:
            shared_list = self.shared_allele_vars
        self.shared_genotype_vars = {"VCF1": [], "VCF2": []}
        self.unshared_genotype_vars = {"VCF1": [], "VCF2": []}
        # samples1 = self.VCF1.samples
        # samples2 = self.VCF2.samples
        if len(self.VCF1.samples) != len(self.VCF2.samples):
            print("Different numbers of samples in each VCF. Cannot compare genotypes.")
            return
        if sorted(self.VCF1.samples) != sorted(self.VCF2.samples):
            if len(self.VCF1.samples) != 1:
                print("Sample name mismatch in multisample VCFs. Cannot compare genotypes.")
                return
            elif len(self.VCF1.samples) == 1:
                print("Sample name mismatch in single sample VCFs. Comparing genotypes anyway.")
                sample_pairs = list(zip(self.VCF1.samples, self.VCF2.samples))
        elif sorted(self.VCF1.samples) == sorted(self.VCF2.samples):
            sample_pairs = list(zip(self.VCF1.samples, self.VCF2.samples))
        for var1, var2 in zip(shared_list["VCF1"], shared_list["VCF2"]):
            geno_mismatches = []
            for sample_pair in sample_pairs:
                geno1 = self.decode_genotype(var1, sample_pair[0])
                geno2 = self.decode_genotype(var2, sample_pair[1])

                if sorted(geno1) != sorted(geno2):
                    geno_mismatches.append(sample_pair[0])

            if len(geno_mismatches) == 0:
                self.shared_genotype_vars["VCF1"].append(var1)
                self.shared_genotype_vars["VCF2"].append(var2)
            elif len(geno_mismatches) > 0:
                self.unshared_genotype_vars["VCF1"].append([var1, geno_mismatches])
                self.unshared_genotype_vars["VCF2"].append([var2, geno_mismatches])
        return

    def compare_by_filter(self, shared_list=None):
        if shared_list is None:
            shared_list = self.shared_genotype_vars
        self.shared_filter_vars = {"VCF1": [], "VCF2": []}
        self.unshared_filter_vars = {"VCF1": [], "VCF2": []}
        for var1, var2 in zip(shared_list["VCF1"], shared_list["VCF2"]):
            if sorted(var1.filter) == sorted(var2.filter):
                self.shared_filter_vars["VCF1"].append(var1)
                self.shared_filter_vars["VCF2"].append(var2)
            else:
                self.unshared_filter_vars["VCF1"].append(var1)
                self.unshared_filter_vars["VCF2"].append(var2)
        return

    def compare_all(self):
        self.compare_by_pos()
        self.compare_by_alleles()
        self.compare_by_genotypes()
        self.compare_by_filter()
        return

    @staticmethod
    def decode_genotype(variant, sample):
        alleles = [variant.ref] + variant.alt
        geno = variant.genotypes[sample][b"GT"].decode()
        if "/" in geno:
            geno = geno.split("/")
        elif "|" in geno:
            geno = geno.split("|")
        if "." in geno:
            return geno
        return [alleles[int(i)] for i in geno]

    def count_shared_by_pos(self):
        return len(self.shared_pos_vars["VCF1"])

    def count_unique_by_pos(self):
        return [len(self.unshared_pos_vars["VCF1"]), len(self.unshared_pos_vars["VCF2"])]

    def count_variants(self):
        return [len(self.VCF1.variants), len(self.VCF2.variants)]

    def count_shared_by_alleles(self):
        return len(self.shared_allele_vars["VCF1"])

    def count_unique_by_alleles(self):
        return [len(self.unshared_allele_vars["VCF1"]), len(self.unshared_allele_vars["VCF2"])]

    def count_shared_by_genotype(self):
        return len(self.shared_genotype_vars["VCF1"])

    def count_unique_by_genotype(self):
        return [len(self.unshared_genotype_vars["VCF1"]), len(self.unshared_genotype_vars["VCF2"])]

    def count_shared_by_filter(self):
        return len(self.shared_filter_vars["VCF1"])

    def count_unique_by_filter(self):
        return [len(self.unshared_filter_vars["VCF1"]), len(self.unshared_filter_vars["VCF2"])]

    def summary_count(self):
        total = self.count_variants()
        if 0 in total:
            print(f"Summary:\tVCF1\tVCF2\n"
                  f"Total calls:\t{total[0]}\t{total[1]}\n"
                  f"\n"
                  f"Empty VCF found. No comparison performed.")
            return
        pos_match = self.count_shared_by_pos()
        pos_diff = self.count_unique_by_pos()
        allele_match = self.count_shared_by_alleles()
        allele_diff = self.count_unique_by_alleles()
        geno_match = self.count_shared_by_genotype()
        geno_diff = self.count_unique_by_genotype()
        filter_match = self.count_shared_by_filter()
        filter_diff = self.count_unique_by_filter()
        return (f"Summary:\tVCF1\tVCF2\n"
                f"\n"
                f"Total calls:\t{total[0]}\t{total[1]}\n"
                f"...of which match by position:\t{pos_match}\t{100*pos_match/total[0]:6.3f} | {100*pos_match/total[1]:6.3f}\n"
                f"...of which match by alleles:\t{allele_match}\t{100*allele_match/total[0]:6.3f} | {100*allele_match/total[1]:6.3f}\n"
                f"...of which match by genotype:\t{geno_match}\t{100*geno_match/total[0]:6.3f} | {100*geno_match/total[1]:6.3f}\n"
                f"...of which match by filter status:\t{filter_match}\t{100*filter_match/total[0]:6.3f} | {100*filter_match/total[1]:6.3f}\n"
                f"\n"
                f"Unique calls:\tVCF1\tVCF2\n"
                f"By positions:\t{pos_diff[0]}\t{pos_diff[1]}\n"
                f"By alleles:\t{allele_diff[0]}\t{allele_diff[1]}\n"
                f"By genotypes:\t{geno_diff[0]}\t{geno_diff[1]}\n"
                f"By filter status:\t{filter_diff[0]}\t{filter_diff[1]}")

    def show_unshared_alleles(self):
        sort_order = dict()
        for chrom in self.chrom_set:
            try:
                sort_order[chrom.decode()] = int(chrom)
            except ValueError:
                if chrom == b"X":
                    sort_order[chrom] = 23
                elif chrom == b"Y":
                    sort_order[chrom] = 24
                elif chrom == b"MT":
                    sort_order[chrom] = 25
                else:
                    sort_order[chrom] = 26
        header = "CHROM:POS\tREF\tALT1 | ALT2\n"
        printer = []
        for variant1, variant2 in zip(self.unshared_allele_vars["VCF1"], self.unshared_allele_vars["VCF2"]):
            alts1 = ",".join([alt.decode() for alt in variant1.alt])
            alts2 = ",".join([alt.decode() for alt in variant2.alt])
            genos1 = ",".join([variant1.genotypes[sample][b"GT"].decode() for sample in self.VCF1.samples])
            genos2 = ",".join([variant2.genotypes[sample][b"GT"].decode() for sample in self.VCF2.samples])
            printer.append(
                    [variant1.chr.decode(), variant1.pos,
                     variant1.ref.decode(), alts1, genos1,
                     variant2.ref.decode(), alts2, genos2]
                    )
        printer.sort(key=lambda x: x[1])
        printer.sort(key=lambda x: sort_order[x[0]])
        printer = [f"{x[0]}:{x[1]}\t{x[2]}\t{x[3]}\t{x[[4]]}\n" +
                   f"\t\t{x[5]}\t{x[6]}\t{x[7]}" for x in printer]
        printer = header + "\n".join(printer)
        return printer

    @staticmethod
    def get_variant_relative_bed_info(BED_object, var_list):
        relatives = []
        for variant in var_list:
            for interval in BED_object.intervals[variant.chr.decode()]:
                if variant.pos in interval:
                    relatives.append([variant.chr.decode(), variant.pos,
                                      interval.start-1, interval.stop-1,
                                      len(interval),
                                      (interval.index(variant.pos)/len(interval))])
                    break
        return relatives


def get_params():
    """Define, receive and return command line parameter values

    Returns
    -------
    dict
        Command line parameter values
    """
    ygor_args = argparse.ArgumentParser()
    ygor_args.add_argument("-1", "--vcf1", dest="vcf1", required=True, type=str, help="Path to the first VCF file")
    ygor_args.add_argument("-2", "--vcf2", dest="vcf2", required=True, type=str, help="Path to the second VCF file")
    ygor_args.add_argument("-o", "--outdir", dest="outdir", required=True, type=str, help="Path to folder to write output files to")
    ygor_args.add_argument("-op", "--out-prefix", dest="outprefix", type=str, help="Prefix to use for output files")
    ygor_args.add_argument("-s", "--summary", dest="summary", action="store_true", type=bool, help="Display the comparison summary?")
    return vars(ygor_args.parse_args())


def params_ok(ygor_arg_values):
    """Check whether the required parameters are ok.

    Parameters
    ----------
    ygor_arg_values : dict
        Ygor command line parameter values.

    Returns
    -------
    params_ok : bool
        True if required parameters are ok, False if not.
    """
    vcf_ext = (".vcf.gz", ".VCF.GZ")
    params_ok = True
    if not os.path.isfile(ygor_arg_values["vcf1"]) and not os.path.isfile(ygor_arg_values["vcf2"]):
        params_ok = False
        print("One (or both) of the supplied VCF files does not exist")
    else:
        if not ygor_arg_values["vcf1"].endswith(vcf_ext) or not ygor_arg_values["vcf2"].endswith(vcf_ext):
            params_ok = False
            print("One (or both) of the supplied is not a gzipped VCF file")
    if not os.path.isdir(ygor_arg_values["outdir"]):
        params_ok = False
        print("The supllied output directory does not exist")
    return params_ok


def main(vcf1, vcf2, allo_check=False):
    vcf1_obj = VCF(vcf1)
    vcf2_obj = VCF(vcf2)
    vcf1_obj.read_from_file()
    vcf2_obj.read_from_file()

    comparison = VCF_Comparison(vcf1_obj, vcf2_obj)
    if not allo_check:
        comparison.compare_by_pos()
        comparison.compare_by_genotypes(comparison.shared_pos_vars)
        comparison.compare_by_filter()
    else:
        comparison.compare_all()
    print(comparison.summary_count())
    print(comparison.unshared_filter_vars)
    return comparison


def main2(allo_check=False):
    """Perform the main comparison functionality.

    Parameters
    ----------
    allo_check : bool
        Perform all checks if True, only position, genotype and filter if False.
    """
    ygor_params = get_params()
    if params_ok(ygor_params):
        vcf1_obj = VCF(ygor_params["vcf1"])
        vcf2_obj = VCF(ygor_params["vcf2"])
        vcf1_obj.read_from_file()
        vcf2_obj.read_from_file()

        comparison = VCF_Comparison(vcf1_obj, vcf2_obj)
        if not allo_check:
            comparison.compare_by_pos()
            comparison.compare_by_genotypes(comparison.shared_pos_vars)
            comparison.compare_by_filter()
        else:
            comparison.compare_all()

        if ygor_params["summary"]:
            print(comparison.summary_count())

        outputdir = ygor_params["outdir"]
        # Start writing shared variant output files
        write_comparison_data(comparison.shared_pos_vars, f"{outputdir}/shared_pos_vars")
        write_comparison_data(comparison.shared_genotype_vars, f"{outputdir}/shared_genotype_vars")
        write_comparison_data(comparison.shared_filter_vars, f"{outputdir}/shared_filter_vars")

        # Start writing unshared variant output files
        write_comparison_data(comparison.unshared_pos_vars, f"{outputdir}/unshared_pos_vars")
        write_comparison_data(comparison.unshared_genotype_vars, f"{outputdir}/unshared_genotype_vars")
        write_comparison_data(comparison.unshared_filter_vars, f"{outputdir}/unshared_filter_vars")

        # Start writing the differences in genotype and filter values between the VCF1 and VCF2 variant on the same position
        write_differences(f"{outputdir}/genotype_differences.txt", comparison.unshared_genotype_vars, "genotype")
        write_differences(f"{outputdir}/filter_differences.txt", comparison.unshared_filter_vars, "filter")

        return comparison
    else:
        print("Not all parameters are ok :(")
        return None


def write_comparison_data(variantdata, outputlocation):
    """Write shared/unshared variants

    Parameters
    ----------
    variantdata : dict
        Shared/Unshared variant data for the two VCF files
    outputlocation : str
        Directory to write output files to
    """
    vcf1_out = f"{outputlocation}_vcf1.txt"
    vcf2_out = f"{outputlocation}_vcf2.txt"
    write_comparison_outfile(variantdata["VCF1"], vcf1_out)
    write_comparison_outfile(variantdata["VCF2"], vcf2_out)


def write_comparison_outfile(singlevars, outfilelocation):
    """Write shared/unshared comparison variants for VCF1 or VCF2.

    Parameters
    ----------
    singlevars : list of str
        Shared/Unshared variants for a single (VCF1/VCF2) vcf file
    outfilelocation : str
        Path to write output file to
    """
    try:
        with open(outfilelocation, 'w') as outfile:
            for variantline in singlevars:
                outfile.write(f"{variantline}\n")
    except IOError:
        print(f"Could not write data to {outfilelocation}")


def write_differences(outfileloc, differing_variants, difference_field):
    """Write variants that differ between VCF1 and VCF2

    Parameters
    ----------
    outfileloc : str
        Path to write differences output file to
    differing_variants : dict
        Variants for VCF1 and VCF2 that differ by genotype or filter
    difference_field : int
        Which variant data field to select

    Returns
    -------
    write_success : bool
        True if the differences file has been succesfully written, False if not
    """
    write_success = False
    header_field = get_difference_field(difference_field)
    vcf1_vars = differing_variants["VCF1"]
    vcf2_vars = differing_variants["VCF2"]
    num_of_iterations = min([len(vcf1_vars), len(vcf2_vars)])

    try:
        with open(outfileloc, 'w') as diff_file:
            diff_file.write(f"Chrom\tPos\tIdentifier\t{difference_field.title()}\n")

            # Start writing the differences to the outoput file
            for noi in range(0, num_of_iterations):
                vcf1_variant = vcf1_vars[noi][0].__repr__()
                vcf2_variant = vcf2_vars[noi][0].__repr__()
                difference_line = get_difference_line(vcf1_variant, vcf2_variant, header_field)
                diff_file.write(f"{difference_line}\n")
        write_success = True
    except IOError:
        print(f"Could not write differences to {outfileloc}")
    finally:
        return write_success


def get_difference_field(differencetype):
    """Determine which VCF line field to write to differences output file.

    Parameters
    ----------
    differencetype : str
        Which difference to write (filter or genotype)

    Returns
    -------
    int
        Field number to write to differencs output file
    """
    difference_types = ["filter", "genotype"]
    if differencetype == "filter":
        return 6
    elif differencetype == "genotype":
        return 4


def get_difference_line(vcf1_variant, vcf2_variant, diff_field):
    """Form andf return the fileline to write to the difference output file.   

    Parameters
    ----------
    vcf1_variant : str
        Variant line for the variant in VCF1
    vcf2_variant : str
        Variant line for the variant in VCF2
    diff_field
        Index of which variant data element to select

    Returns
    -------
    str or None
        Difference file line to write or None if the variants are not the same
    """
    vcf1_variantdata = vcf1_variant.split("\t")
    vcf2_variantdata = vcf2_variant.split("\t")

    if vcf1_variantdata[0] == vcf2_variantdata[0] and vcf1_variantdata[1] == vcf2_variantdata[1]:
        return f"{vcf1_variantdata[0]}\t{vcf1_variantdata[1]}\t{vcf1_variantdata[2]}\t{vcf1_variantdata[diff_field]}|{vcf2_variantdata[diff_field]}"
    else:
        return None


if __name__ == "__main__":
    # Comparator = main(sys.argv[1], sys.argv[2])
    Comparator = main2()
