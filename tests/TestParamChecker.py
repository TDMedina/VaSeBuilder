#!/usr/bin/env python
# Import required modules/libraries 
import unittest
import os
import sys

# Import required class
sys.path.insert(0, './../')
from ParamChecker import ParamChecker


# Unittest class for the ParamChecker test.
class TestParamChecker(unittest.TestCase):
    # Set up required 
    def setUp(self):
        self.param_check = ParamChecker()
        self.param_list = {"donorvcf": "testdata/vcfDir/vcflistfile.txt",
                           "donorbam": "testdata/bamDir/bamlistfile.txt",
                           "acceptorbam": "testdata/valbam/SRR1039513.bam",
                           "templatefq1": "testdata/fqDir/SRR1039513_1.fastq.gz",
                           "templatefq2": "testdata/fqDir/SRR1039513_2.fastq.gz",
                           "out": "testdata/outDir",
                           "reference": "testdata/ref/reference.fa",
                           'fastqout': "testdata/outDir",
                           'varcon': "testdata/outDir/varcon.txt",
                           "variantlist": "testdata/variantlist.txt"
                           }
        self.default_log_answer = "VaSeBuilder.log"
        self.default_varcon_answer = "varcon.txt"
        self.required_parameters = {"donorvcf": "testdata/vcfDir/vcflistfile.txt",
                                    "donorbam": "testdata/bamDir/bamlistfile.txt",
                                    "acceptorbam": "testdata/valbam/SRR1039513.bam",
                                    "templatefq1": "testdata/fqDir/SRR1039513_1.fastq.gz",
                                    "templatefq2": "testdata/fqDir/SRR1039513_2.fastq.gz",
                                    "out": "testdata/outDir",
                                    "reference": "testdata/ref/reference.fa",
                                    "varconin": "testdata/varcon.txt"
                                    }

    # Tests that the --log will be written to the log file.
    def test_check_log_pos(self):
        log_pos_answer = "testdata/outDir/vaseutest.log"
        self.assertEqual(self.param_check.check_log(log_pos_answer), log_pos_answer,
                         f"The log location should have been {log_pos_answer}")

    # Tests that the --log parameter has not been set and the log will therfore be written to VaSeBuilder.log 
    def test_check_log_no_param_set(self):
        self.assertEqual(self.param_check.check_log(None), self.default_log_answer,
                         f"The log file location should have been {self.default_log_answer}")

    # Tests that the --log parameter has been set but the filename does not end with .log and will therefore be written
    # to VaSeBuilder.log
    def test_check_log_no_log(self):
        invalid_log_name = "testdata/outDir/logfile.file"
        self.assertEqual(self.param_check.check_log(invalid_log_name), self.default_log_answer,
                         f"The log location should have been {self.default_log_answer}")

    # Test that the --log parameter has been set but the filename does not end with .txt  and will therefore be written
    # to VaSeBuilder.log
    def test_check_log_no_txt(self):
        invalid_log_name = "testdata/outDir/logfile.file"
        self.assertEqual(self.param_check.check_log(invalid_log_name), self.default_log_answer,
                         f"The log location should have been {self.default_log_answer}")

    # Tests that the povided testdata folder 'bamDir' indeed contains two BAM files.
    def test_check_folder_contents_pos(self):
        self.assertEqual(self.param_check.check_folder_contents("testdata/bamDir", "bam"), 2,
                         "Two bam files should have been found")

    # Tests that the provided testdata folder 'bamDir' contains no VCF files
    def test_check_folder_contents_neg(self):
        self.assertEqual(self.param_check.check_folder_contents("testdata/bamDir", "vcf"), 0,
                         "No VCF files should have been found")

    # Tests that the provided testdata folder 'bamDir' exists and contains BAM files.
    def test_check_folders_exist_pos(self):
        self.assertListEqual(self.param_check.check_folders_exist(['testdata/bamDir'], 'bam'), ['testdata/bamDir'],
                             "The folder should have been testdata/bamDir")

    # Tests that the provided testdata folder 'bamDir' exists and does not contain VCF files.
    def test_check_folders_exist_neg(self):
        self.assertListEqual(self.param_check.check_folders_exist(['testdata/bamDir'], 'vcf'), [],
                             "There should be no folder")

    # Tests that a provided file indeed exists.
    def test_check_file_exists_pos(self):
        existing_file = "testdata/bamDir/SRR1039508.bam"
        self.assertTrue(self.param_check.check_file_exists(existing_file), f"{existing_file} should exist")

    # Tests whether a non existing file indeed does not exist.
    def test_check_file_exists_neg(self):
        nonexisting_file = ""
        self.assertFalse(self.param_check.check_file_exists(nonexisting_file),
                         f"Non existing file {nonexisting_file} should not exist")

    # Test that the output location for out.txt is valid.
    def test_is_valid_output_location_pos(self):
        valid_out_loc = "testdata/outDir/out.txt"
        self.assertTrue(self.param_check.is_valid_output_location(valid_out_loc),
                        f"{valid_out_loc} should have been a valid output location")

    # Test that the output location for out.txt is invalid.
    def test_is_valid_output_location_neg(self):
        invalid_out_loc = "testdata/doesnotexist/out.txt"
        self.assertFalse(self.param_check.is_valid_output_location(invalid_out_loc),
                         f"{invalid_out_loc} should have been an invalid output location")

    # Tests that the foldername a file is located in gets returned properly.
    def test_get_folder_name_pos(self):
        valid_path = "testdata/noSampleDir/noSampleBam.bam"
        valid_path_answer = "testdata/noSampleDir"
        self.assertEqual(self.param_check.get_folder_name(valid_path), valid_path_answer,
                         f"Folder name should have been {valid_path_answer}")

    # Test that the foldername is empty.
    def test_get_folder_name_neg(self):
        self.assertEqual(self.param_check.get_folder_name(""), "", "Folder names should both be empty.")

    # Tests that all provided parameters are ok.
    # First test should return True, all others False due to one parameter missing
    def test_check_parameters_pos(self):
        self.assertTrue(self.param_check.check_parameters(self.param_list))

    def test_check_parameters_no_donor_vcf(self):
        par_list = self.param_list.copy()
        par_list['donorvcf'] = 'testdata/doesnotexist.txt'    # Set donorvcf parameter tp nonexisting vcf list file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_donor_bam(self):
        par_list = self.param_list.copy()
        par_list['donorbam'] = 'testdata/doesnotexist.txt'    # Set donorbam parameter to nonexisting bam list file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_acceptor_bam(self):
        par_list = self.param_list.copy()
        par_list["acceptorbam"] = "testdata/doesnotexist.bam"    # Set acceptorbam parameter to nonexisting bam file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_acceptor_fq1(self):
        par_list = self.param_list.copy()
        par_list["templatefq1"] = "testdata/doesnotexist.fq"    # Set templatefq1 parameter to nonexisting fastq file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_acceptor_fq2(self):
        par_list = self.param_list.copy()
        par_list['templatefq2'] = 'testdata/doesnotexist.fq'    # Set templatefq2 parameter to nonexisting fastq file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_fq_out(self):
        par_list = self.param_list.copy()
        par_list["fastqout"] = "testdata/doesnotexist/faqOut.fq"    # Set fastqout parameter to nonexisting location
        self.assertTrue(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_varcon(self):
        par_list = self.param_list.copy()
        par_list["varcon"] = "testdata/doesnotexist/varcon.txt"    # Set varcon parameter to nonexisting location.
        self.assertTrue(self.param_check.check_parameters(par_list))

    # Tests that the VCF file list location is saved correctly in ParamChecker
    def test_get_valid_vcf_filelist(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_valid_vcf_filelist(), self.param_list["donorvcf"],
                         "The saved donor vcf list file location should have been " + self.param_list["donorvcf"])

    # Tests that the BAM file list location is saved correctly in ParamChecker
    def test_get_valid_bam_filelist(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_valid_bam_filelist(), self.param_list["donorbam"],
                         "The saved donor bam list file location should have been " + self.param_list["donorbam"])

    # Tests that the acceptor BAM/CRAM location is saved correctly in ParamChecker
    def test_get_acceptor_bam(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_acceptor_bam(), self.param_list["acceptorbam"],
                         "The saved acceptor BAM file should have been " + self.param_list["acceptorbam"])

    # Tests that the out directory location is saved correctly in ParamChecker
    def test_get_out_dir_location(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_out_dir_location(), self.param_list["out"],
                         "The saved out directory location should have been " + self.param_list["out"])

    # Tests that the reference is saved correctly in ParamChecker
    def test_get_reference_file_location(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_reference_file_location(), self.param_list["reference"],
                         "The saved reference file location should have been " + self.param_list["reference"])

    # Tests that the varcon out location is saved correctly in ParamChecker
    def test_get_variant_context_out_location(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_variant_context_out_location(), self.param_list["varcon"],
                         "The saved varcon out location should have been " + self.param_list["varcon"])

    # Tests that the default log file location is saved correctly in ParamChecker
    def test_get_log_file_location_default(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_log_file_location(), self.default_log_answer,
                         f"The saved default log file location should have been {self.default_log_answer}")

    # Tests that a provided log file location is saved correctly in ParamChecker
    def test_get_log_file_location(self):
        logloc = "testdata/outDir/testlog.log"
        self.assertEqual(self.param_check.check_log(logloc), logloc,
                         f"The saved log file location should have been {logloc}")

    # Tests that the variant list location is saved correctly in ParamChecker
    def test_get_variant_list_location(self):
        self.param_check.check_parameters(self.param_list)
        self.assertEqual(self.param_check.get_variant_list_location(), self.param_list["variantlist"],
                         "The saved variant list file location should have been " + self.param_list["variantlist"])

    # Tests that all required prameters have been set for runmode 'A'
    def test_required_parameters_set_amode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("donorvcf", None)
        set_parameters.pop("donorbam", None)
        set_parameters.pop("acceptorbam", None)
        set_parameters.pop("reference", None)
        self.assertTrue(self.param_check.required_parameters_set("A", set_parameters),
                        "All required parameters for A mode should have been ok and therefore return True")

    # Tests that all required paramneters have been set for runmode 'D'
    def test_required_parameters_set_dmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("templatefq1", None)
        set_parameters.pop("templatefq2", None)
        set_parameters.pop("varconin", None)
        self.assertTrue(self.param_check.required_parameters_set("D", set_parameters),
                        "All required parameters for D mode should have been ok and therefore return True")

    # Tests that all required parmeters have been set for runmode 'DC'
    def test_required_parameters_set_dcmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("templatefq1", None)
        set_parameters.pop("templatefq2", None)
        self.assertTrue(self.param_check.required_parameters_set("DC", set_parameters),
                        "All required parameters for DC mode should have been ok and therefore return True")

    # Tests that the required parameters have been set for runmode 'F'
    def test_required_parameters_set_fmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("varconin", None)
        self.assertTrue(self.param_check.required_parameters_set("F", set_parameters),
                        "All required parameters for F mode should have been ok and therefore return True")

    # Tests that the required parameters have been set for runmode 'FC'
    def test_required_parameters_set_fcmode(self):
        set_parameters = self.required_parameters.copy()
        self.assertTrue(self.param_check.required_parameters_set("FC", set_parameters),
                        "All required parameters for FC mode should have been ok and therefore return True")

    # Tests that the required parameters have been set for runmode 'P'
    def test_required_parameters_set_pmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("templatefq1", None)
        set_parameters.pop("templatefq2", None)
        set_parameters.pop("varconin", None)
        self.assertTrue(self.param_check.required_parameters_set("P", set_parameters),
                        "All required parameters for P mode should have been ok and therefore return True")

    # Tests that the required parameters have been set for runmode 'PC'
    def test_required_parameters_set_pcmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("templatefq1", None)
        set_parameters.pop("templatefq2", None)
        self.assertTrue(self.param_check.required_parameters_set("PC", set_parameters),
                        "All required parameters for PC mode should have been ok and therefore return True")

    # Tests that the required parameters have been set for runmode 'X'
    def test_required_parameters_set_xmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("templatefq1", None)
        set_parameters.pop("templatefq2", None)
        set_parameters.pop("varconin", None)
        self.assertTrue(self.param_check.required_parameters_set("X", set_parameters),
                        "All required parameters for X mode should have been ok and therefore return True")

    # Tests that the required parameters have been set for runmode 'XC'
    def test_required_parameters_set_xcmode(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("templatefq1", None)
        set_parameters.pop("templatefq2", None)
        self.assertTrue(self.param_check.required_parameters_set("XC", set_parameters),
                        "All required parameters for XC mode should have been ok and therefore return True")

    # Tests that False is returned when an incorrect runmode is specified
    def test_required_parameters_set_invalidmode(self):
        runmode = "G"
        self.assertFalse(self.param_check.required_parameters_set(runmode, {}), "False should have been returned")

    # Tests that False is returned when not all required parameters are set for F mode
    def test_required_parameters_set_notsetamodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("donorfastqs", None) # Removing the required donorfastqs parameter
        self.assertFalse(self.param_check.required_parameters_set("A", set_parameters),
                         "Not all required parameters for A mode should have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for AC mode is not set
    def test_required_parameters_set_notsetacmodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("varconin", None)    # Remove the required varconin parameter
        self.assertFalse(self.param_check.required_parameters_set("AC", set_parameters),
                         "The required parameters for AC mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for D mode is not set
    def test_required_parameters_set_notsetdmodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("", None)    # Removing the required parameter
        self.assertFalse(self.param_check.required_parameters_set("D", set_parameters),
                         "The required parameters for D mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for DC mode in not set
    def test_required_parameters_set_notsetdcmodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("varconin", None)    # Removing the required varconin parameter
        self.assertFalse(self.param_check.required_parameters_set("DC", set_parameters),
                         "The required parameters for DC mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for F mode is not set
    def test_required_parameters_set_notestfmodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("donorvcf", None)
        self.assertFalse(self.param_check.required_parameters_set("F", set_parameters),
                         "The required parameters for F mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for FC mode is not set
    def test_required_parameters_set_notsetfcmodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("varconin", None)    # Removing the required varconin parameter
        self.assertFalse(self.param_check.required_parameters_set("FC", set_parameters),
                         "The required parameters for FC mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for P mode is not set
    def test_required_parameters_set_invalidpmodeparam(self):
        set_parameters = self.required_parameters.copy()
        set_parameters.pop("", None)
        self.assertFalse(self.param_check.required_parameters_set("P", set_parameters),
                         "The required parameters for P mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for PC mode is not set
    def test_required_parameters_set_invalidpcmodeparam(self):
        set_parameters = {}
        self.assertFalse(self.param_check.required_parameters_set("PC", set_parameters),
                         "The required parameters for PC mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for X mode is not set
    def test_required_parameters_set_invalidxmodeparam(self):
        set_parameters = {}
        self.assertFalse(self.param_check.required_parameters_set("X", set_parameters),
                         "The required parameters for X mode should not have been ok and therefore return False")

    # Tests that False is returned when one of the required parameters for XC mode is not set
    def test_required_parameters_set_invalidxcmodeparam(self):
        set_parameters = {}
        self.assertFalse(self.param_check.required_parameters_set("XC", set_parameters),
                         "The required parameters for XC mode should not have been ok and therefore return False")

    # =====TEST SOME OTHER METHODS=====
    def test_get_required_runmode_parameters_acmode(self):
        ac_params_answer = ["runmode", "templatefq1", "templatefq2", "donorfastqs", "varconin", "out"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("AC"), ac_params_answer,
                             f"The list of required AC-mode parameters should have been {ac_params_answer}")

    def test_get_required_runmode_parameters_dmode(self):
        d_params_answer = ["runmode", "donorvcf", "donorbam", "acceptorbam", "out", "reference"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("D"), d_params_answer,
                             f"The list of required D-mode parameters should have been {d_params_answer}")

    def test_get_required_runmode_parameters_dcmode(self):
        dc_params_answer = ["runmode", "donorvcf", "donorbam", "out", "reference", "varconin"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters(""), dc_params_answer,
                             f"The list of required DC-mode parameters should have been {dc_params_answer}")

    def test_get_required_runmode_parameters_fmode(self):
        f_params_answer = ["runmode", "donorvcf", "donorbam", "acceptorbam", "templatefq1", "templatefq2", "out",
                           "reference"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("F"), f_params_answer,
                             f"The list of required F-mode parameters should have been {f_params_answer}")

    def test_get_required_runmode_parameters_fcmode(self):
        fc_params_answer = ["runmode", "donorvcf", "donorbam", "templatefq1", "templatefq2", "out", "reference",
                            "varconin"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("FC"), fc_params_answer,
                             f"the list of required FC-mode parameters should have been {fc_params_answer}")

    def test_get_required_runmode_parameters_pmode(self):
        p_params_answer = ["runmode", "donorvcf", "donorbam", "acceptorbam", "out", "reference"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("P"), p_params_answer,
                             f"The list of required parameters should have been {p_params_answer}")

    def test_get_required_runmode_parameters_pcmode(self):
        pc_params_answer = ["runmode", "donorvcf", "donorbam", "out", "reference", "varconin"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("PC"), pc_params_answer,
                             f"The list of required parameters should have been {pc_params_answer}")

    def test_get_required_runmode_parameters_xmode(self):
        x_params_answer = ["runmode", "donorvcf", "donorbam", "acceptorbam", "out", "reference"]
        self.assertListEqual(self.param_check.get_required_runmode_parameters("X"), x_params_answer,
                             f"The list of required parameters should have been {x_params_answer}")

    def test_get_required_runmode_parameters_nonexistent(self):
        self.assertListEqual(self.param_check.get_required_runmode_parameters("T"), [],
                             f"The list of required parameters for a non existent mode should have been empty.")

    def test_get_optional_parameters(self):
        opt_params_answer = ["fastqout", "varcon", "variantlist", "seed"]
        self.assertListEqual(self.param_check.get_optional_parameters(), opt_params_answer,
                             f"The list of optional parameters should have been {opt_params_answer}")
