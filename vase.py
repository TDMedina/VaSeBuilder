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
from VariantContextFile import VariantContextFile


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

        # Check whether a configuration file was supplied than all required command line parameters.
        if vase_arg_list["configfile"] is not None:
            configfileloc = vase_arg_list["configfile"]
            vase_arg_list = self.read_config_file(configfileloc)    # Set the read config file as the parameter list
            # Check whether the DEBUG parameters has been set or not
            if "debug" not in vase_arg_list:
                vase_arg_list["debug"] = False
            # Check optional parameters and set those missing to None
            vase_arg_list = pmc.optional_parameters_set(vase_arg_list)

        # Start the logger and initialize this run with an ID number.
        self.vaselogger = self.start_logger(pmc, vase_arg_list["log"], vase_arg_list["debug"])

        vase_called_command = " ".join(sys.argv)
        self.vaselogger.info(f"python {vase_called_command}")
        vase_b = VaSeBuilder(uuid.uuid4().hex)

        # Exit if not all of the required parameters have been set
        if not pmc.required_parameters_set(vase_arg_list["runmode"], vase_arg_list):
            self.vaselogger.critical("Not all required parameters have been set")
            exit()

        # Exit if the supplied parameters are incorrect.
        if not pmc.check_required_runmode_parameters(vase_arg_list["runmode"], vase_arg_list):
            self.vaselogger.critical("Not all required parameters are correct. Please check log for more info.")
            exit()

        # Check if a variantfilter has been set.
        variantfilter = None
        if pmc.get_variant_list_location() != "":
            variantfilter = self.read_variant_list(pmc.get_variant_list_location())

        # Run the selected mode.
        self.run_selected_mode(pmc.get_runmode(), vase_b, pmc, variantfilter)

        self.vaselogger.info("VaSeBuilder run completed succesfully.")
        elapsed = time.strftime(
                "%Hh:%Mm:%Ss",
                time.gmtime(time.time() - vase_b.creation_time.timestamp())
                )
        self.vaselogger.info(f"Elapsed time: {elapsed}.")

    # Method that creates the logger that will write the log to stdout and a log file.
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
        debug_param_vals = ["True", "1", "T"]
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
                            # Check if the parameter is the debug parameter
                            elif parameter_name == "debug":
                                configdata[parameter_name] = parameter_value.title() in debug_param_vals
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

    def run_selected_mode(self, runmode, vaseb, paramcheck, variantfilter):
        vbscan = VcfBamScanner()
        varconfile = None    # Declare the varconfile variable so we can use it in C and no-C.

        if "X" in runmode:
            vcf_file_map = vbscan.scan_vcf_files(paramcheck.get_valid_vcf_filelist())
            bam_file_map = vbscan.scan_bamcram_files(paramcheck.get_valid_bam_filelist())
            sample_id_list = vbscan.get_complete_sample_ids()
            vaseb.run_x_mode(sample_id_list, vcf_file_map, bam_file_map, paramcheck.get_acceptor_bam(),
                             paramcheck.get_reference_file_location(), paramcheck.get_out_dir_location(),
                             paramcheck.get_variant_context_out_location(), variantfilter)
        else:
            if "C" in runmode:
                varconfile = VariantContextFile(paramcheck.get_variantcontext_infile())
                if "A" in runmode:
                    donor_fastq_files = self.read_donor_fastq_list_file(paramcheck.get_donorfqlist())
                    vaseb.run_ac_mode(paramcheck.get_first_fastq_in_location(),
                                      paramcheck.get_second_fastq_in_location(),
                                      donor_fastq_files, varconfile, paramcheck.get_out_dir_location())
            else:
                # Scan the variant and alignment files in the provided lists.
                vcf_file_map = vbscan.scan_vcf_files(paramcheck.get_valid_vcf_filelist())
                bam_file_map = vbscan.scan_bamcram_files(paramcheck.get_valid_bam_filelist())
                sample_id_list = vbscan.get_complete_sample_ids()
                varconfile = vaseb.build_varcon_set(sample_id_list, vcf_file_map, bam_file_map,
                                                    paramcheck.get_acceptor_bam(), paramcheck.get_out_dir_location(),
                                                    paramcheck.get_reference_file_location(),
                                                    paramcheck.get_variant_context_out_location(), variantfilter)
            # Check for modes D,F,P
            if "D" in runmode:
                vaseb.run_d_mode(varconfile, paramcheck.get_fastq_out_location())
            if "F" in runmode:
                vaseb.run_f_mode(varconfile, paramcheck.get_first_fastq_in_location(),
                                 paramcheck.get_second_fastq_in_location(), paramcheck.get_fastq_out_location())
            if "P" in runmode:
                vaseb.run_p_mode(varconfile, paramcheck.get_fastq_out_location())


# Run the program.
vaseb = VaSe()
vaseb.main()
