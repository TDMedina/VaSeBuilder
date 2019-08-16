#!/usr/bin/env python

# Import necessary modules.
import logging
import sys
import uuid
import argparse
import pysam
import time

# Import VaSe classes.
from ParamChecker import ParamChecker
from VcfBamScanner import VcfBamScanner
from VaSeBuilder import VaSeBuilder


class VaSe:
    # Performs the check that VaSe is run with Python 3.x
    def __init__(self):
        assert (sys.version_info[0] >= 3 and sys.version_info[1] >= 6), "Please run this program in Python 3.6 or " \
                                                                        "higher"
        assert (int(pysam.version.__version__.split(".")[0]) >= 0 and int(pysam.version.__version__.split(".")[1]) >=
                15), "Please run this program with Pysam 0.15 or higher"
        self.valid_runmodes = ["A", "AC", "D", "DC", "F", "FC", "P", "PC", "X", "XC"]

    # Runs the program.
    def main(self):
        # Parse the command line parameters and check their validity.
        vase_arg_list = self.get_vase_parameters()
        pmc = ParamChecker()
        # Start the logger and initialize this run with an ID number.
        self.vaselogger = self.start_logger(pmc, vase_arg_list["log"], vase_arg_list["debug"])

        vase_called_command = " ".join(sys.argv)
        self.vaselogger.info(f"python {vase_called_command}")
        vase_b = VaSeBuilder(uuid.uuid4().hex)

        # Check whether a configuration file was supplied than all required command line parameters.
        if vase_arg_list["configfile"] is not None:
            configfileloc = vase_arg_list["configfile"]
            vase_arg_list = self.read_config_file(configfileloc)    # Set the read config file as the parameter list
            # Check optional parameters and set those missing to None
            vase_arg_list = pmc.optional_parameters_set(vase_arg_list)

        # Exit if not all of the required parameters have been set
        if not pmc.required_parameters_set(vase_arg_list["runmode"], vase_arg_list):
            self.vaselogger.critical("Not all required parameters have been set")
            exit()

        # Exit if the supplied parameters are incorrect.
        if not pmc.check_required_runmode_parameters(vase_arg_list["runmode"], vase_arg_list):
            self.vaselogger.critical("Not all required parameters are correct. Please check log for more info.")
            exit()

        if "A" in pmc.runmode:
            donor_fastq_files = self.read_donor_fastq_list_file(pmc.get_donorfqlist())
            r1_dfqs = [dfq[0] for dfq in donor_fastq_files]
            r2_dfqs = [dfq[1] for dfq in donor_fastq_files]
            vase_b.build_validation_from_donor_fastqs(pmc.get_first_fastq_in_location(),
                                                      pmc.get_first_fastq_in_location(),
                                                      r1_dfqs, r2_dfqs, pmc.varconin, pmc.get_fastq_out_location())
        else:
            # Scan the variant and alignment files in the provided lists.
            vbscan = VcfBamScanner()
            vcf_file_map = vbscan.scan_vcf_files(vase_arg_list["donorvcf"])
            bam_file_map = vbscan.scan_bamcram_files(vase_arg_list["donorbam"])
            sample_id_list = vbscan.get_complete_sample_ids()

            # Exit if no valid pairs of variant/alignment files are found.
            if not sample_id_list:
                self.vaselogger.critical("No valid samples available to "
                                         "create new validation set")
                exit()

            variantfilter = None
            if pmc.get_variant_list_location() != "":
                variantfilter = self.read_variant_list(pmc.get_variant_list_location())

            if "C" not in pmc.runmode:
                vase_b.build_varcon_set(sample_id_list,
                                        vcf_file_map, bam_file_map,
                                        pmc.get_acceptor_bam(),
                                        pmc.get_out_dir_location(),
                                        pmc.get_reference_file_location(),
                                        pmc.get_variant_context_out_location(),
                                        variantfilter)
            elif "C" in pmc.runmode:
                vase_b.build_donor_from_varcon(pmc.varconin,
                                               bam_file_map,
                                               pmc.get_reference_file_location())
            if "X" not in pmc.runmode:
                vase_b.build_validation_set(pmc.runmode,
                                            pmc.get_acceptor_bam(),
                                            pmc.get_first_fastq_in_location(),
                                            pmc.get_second_fastq_in_location(),
                                            pmc.get_fastq_out_location())

        self.vaselogger.info("VaSeBuilder run completed succesfully.")
        elapsed = time.strftime(
                "%Hh:%Mm:%Ss",
                time.gmtime(time.time() - vase_b.creation_time.timestamp())
                )
        self.vaselogger.info(f"Elapsed time: {elapsed}.")

    # Method that creates the logger that will write the log to stdout
    # and a log file.
    def start_logger(self, paramcheck, logloc, debug_mode=False):
        vaselogger = logging.getLogger("VaSe_Logger")
        if debug_mode:
            vaselogger.setLevel(logging.DEBUG)
        else:
            vaselogger.setLevel(logging.INFO)
        vaselog_format = logging.Formatter("%(asctime)s	%(name)s	%(levelname)s	%(message)s")

        # Add the log stream to stdout.
        vase_cli_handler = logging.StreamHandler(sys.stdout)
        if debug_mode:
            vase_cli_handler.setLevel(logging.DEBUG)
        else:
            vase_cli_handler.setLevel(logging.INFO)
        vase_cli_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_cli_handler)

        # Create the log stream to log file.
        logloc = paramcheck.check_log(logloc)
        if logloc == "":
            logloc = "VaSeBuilder.log"
        vase_file_handler = logging.FileHandler(logloc)

        if debug_mode:
            vase_file_handler.setLevel(logging.DEBUG)
        else:
            vase_file_handler.setLevel(logging.INFO)
        vase_file_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_file_handler)
        return vaselogger

    # Returns the vase Parameters.
    def get_vase_parameters(self):
        # Set the VaSe parameters for the program.
        vase_argpars = argparse.ArgumentParser()
        vase_argpars.add_argument("-m", "--runmode", dest="runmode", default="F", choices=self.valid_runmodes, help="RUNMODE HELP")
        vase_argpars.add_argument("-v", "--donorvcf", dest="donorvcf", help="File containing a list of VCF/VCF.GZ/BCF files.")
        vase_argpars.add_argument("-b", "--donorbam", dest="donorbam", help="File containing a list of BAM/CRAM files.")
        vase_argpars.add_argument("-a", "--acceptorbam", dest="acceptorbam", help="BAM file for identifying acceptor reads to exclude.")
        vase_argpars.add_argument("-1", "--templatefq1", dest="templatefq1", nargs="*", help="Location and name of the first fastq in file.")
        vase_argpars.add_argument("-2", "--templatefq2", dest="templatefq2", nargs="*", help="Location and name of the second fastq in file.")
        vase_argpars.add_argument("-o", "--out", dest="out", help="Directory to write output files to.")
        vase_argpars.add_argument("-r", "--reference", dest="reference", help="Location of the reference genome. This reference genome should be used by all VCF/BCF and BAM/CRAM files.")
        vase_argpars.add_argument("-of", "--fastqout", dest="fastqout", help="Name for the two FastQ files to be produced.")
        vase_argpars.add_argument("-ov", "--varcon", dest="varcon", help="File name to write variants and their contexts to.")
        vase_argpars.add_argument("-l", "--log", dest="log", help="Location to write log files to (will write to working directory if not used).")
        vase_argpars.add_argument("-!", "--debug", dest="debug", action="store_true", help="Run the program in debug mode")
        vase_argpars.add_argument("-vl", "--variantlist", dest="variantlist", help="File containing a list of variants to use. Will only use these variants if provided. Will use all variants if no list is provided.")
        vase_argpars.add_argument("-iv", "--varconin", dest="varconin", help="Provide a Vasebuilder output variant context file to build a validation set.")
        vase_argpars.add_argument("-dq", "--donorfastqs", dest="donorfastqs", help="Location to donor fastq list file")
        vase_argpars.add_argument("-c", "--config", dest="configfile", help="Supply a config file")
        vase_args = vars(vase_argpars.parse_args())
        return vase_args

    # Reads the variant list. Assumes that sampleId, chromosome and startpos are columns 1,2 and 3
    def read_variant_list(self, variantlistloc):
        variant_filter_list = {}
        try:
            with open(variantlistloc) as variantlistfile:
                next(variantlistfile)    # Skip the header line
                for fileline in variantlistfile:
                    filelinedata = fileline.strip().split("\t")
                    if filelinedata[0] not in variant_filter_list:
                        variant_filter_list[filelinedata[0]] = []
                    variant_filter_list[filelinedata[0]].append((filelinedata[1], int(filelinedata[2])))
        except IOError:
            self.vaselogger.critical(f"Could not open variant list file {variantlistloc}")
        finally:
            return variant_filter_list

    # Reads the config file with the settings
    def read_config_file(self, configfileloc):
        configdata = {}
        try:
            with open(configfileloc, "r") as configfile:
                for fileline in configfile:
                    fileline = fileline.strip()
                    if not fileline.startswith("#"):
                        configentry = fileline.split("=")
                        if len(configentry) == 2:
                            parameter_name = configentry[0].strip().lower()
                            parameter_value = configentry[1].strip()

                            # Check whether the current parameter equals either 'templatefq1' or 'templatefq2'
                            if parameter_name == "templatefq1" or parameter_name == "templatefq2":
                                template_files = parameter_value.split(",")
                                configdata[parameter_name] = [tmplfile.strip() for tmplfile in template_files]
                            else:
                                configdata[parameter_name] = parameter_value.strip()
        except IOError:
            self.vaselogger.critical(f"Could not read configuration file: {configfileloc}")
        return configdata

    # Reads a list of donor fastq files in the format (R1.fq\tR2.fq)
    def read_donor_fastq_list_file(self, donorfq_listfileloc):
        donor_fastqs = []
        try:
            with open(donorfq_listfileloc) as donorfqlistfile:
                for fileline in donorfqlistfile:
                    donor_fastqs.append(fileline.strip().split("\t"))
        except IOError:
            self.vaselogger.warning(f"Could not read donor fastq list file {donorfq_listfileloc}")
        finally:
            return donor_fastqs


# Run the program.
vaseb = VaSe()
vaseb.main()
