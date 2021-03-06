# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 17:47:12 2020

@author: tdmedina
"""

import logging
import os
import sys
import subprocess
import pysam


class VcfBamScanner:
    """VcfBamScanner offers functionality to scan variant (VCF/BCF) and alignment (BAM/CRAM) files.

    The VcfBamScanner can be used to scan folders for variant and alignments files, check whether variant and alignment
    files have sample names and extract them. A list file containing paths to either variant or alignment files can
    also be.scanned. Files within the list files will be a saved per sample name in an identifier.

    Attributes
    ----------
    vcf_sample_map : dict
        Dictionary with variant files per sample
    bam_sample_map : dict
        Dictionary with alignment files per sample
    valid_file_types : list of str
        Valid variant and alignment file types
    """

    def __init__(self):
        """Obtain the vase logger and creates two empty dictionaries to save the variant and alignment files."""
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_sample_map = {}
        self.bam_sample_map = {}
        self.valid_file_types = ["VCF", "BCF", "BAM", "CRAM"]

    # ===METHODS TO GET SAVED DATA FROM THE VCFBAMSCANNER======================
    def get_vcf_sample_map(self):
        """Return a dictionary with VCF/BCF files saved per sample.

        Returns
        -------
        vcf_sample_map : dict
            Paths to variant files per sample name
        """
        return self.vcf_sample_map

    def get_bam_sample_map(self):
        """Return a dictionary with BAM/CRAM files saved per sample.

        Returns
        -------
        bam_sample_map : dict
            Paths to alignment files per sample name
        """
        return self.bam_sample_map

    def get_vcf_to_bam_map(self):
        """Return a dictionary with a variant file linked to an alignment file of the same sample.

        For each sample, the associated variant file is linked to the associated.alignment file are linked via a
        dictionary. The variant and alignment file are saved as key and value respectively.

        Returns
        -------
        vcf_to_bam_map : dict
            Variant files linked to alignment files
        """
        vcf_to_bam_map = {}
        for sampleid in self.vcf_sample_map:
            if sampleid in self.bam_sample_map:
                vcf_to_bam_map[self.vcf_sample_map[sampleid]] = self.bam_sample_map[sampleid]
        return vcf_to_bam_map

    # ===METHODS TO SCAN VCF/BAM FOLDERS AND FILES=============================
    # Scans the folders containing VCF files and returns a map that
    # links sample ids with vcf files.
    def scan_vcf_folders(self, vcffolders):
        """Scan a list of folders containing variant files and return a dictionary with valid files per sample name.

        Iterates over a list of paths to folders containing VCF/BCF files. A non-existing folder is skipped. For each
        variant file in a valid folder, the sample name is extracted.and the file is saved under its sample name in a
        dictionary.

        Parameters
        ----------
        vcffolders : list
            List of paths to folders contain alignment files

        Returns
        -------
        vcf_sample_map : dict
            String paths to variant files per sample name

        """
        self.vaselogger.info("Start scannig VCF files")

        for vcffolder in vcffolders:
            if not os.path.isdir(vcffolder):
                self.vaselogger.info(f"Folder {vcffolder} does not exist.")
                continue
            self.vaselogger.info(f"Scanning VCF files in {vcffolder}")
            for vcf_filename in os.listdir(vcffolder):
                if vcf_filename.endswith((".vcf", ".vcf.gz", ".bcf")):
                    # Do something to scan the file for the sample.
                    try:
                        vcffile = pysam.VariantFile(
                            f"{vcffolder}/{vcf_filename}",
                            "r"
                            )
                        self.vaselogger.debug(
                            "Scanning VCF file "
                            f"{vcffolder}/{vcf_filename}"
                            )
                        if len(vcffile.header.samples) > 0:
                            # This is the sample identifier.
                            sampleid = vcffile.header.samples[0]
                            self.vcf_sample_map[sampleid] = (
                                f"{vcffolder}/{vcf_filename}"
                                )
                        vcffile.close()
                    except IOError:
                        self.vaselogger.warning(
                            f"VCF file {vcffolder}/{vcf_filename} "
                            "does not exist"
                            )
        self.vaselogger.info("Finished scanning VCF files")
        return self.vcf_sample_map

    def scan_bam_folders(self, bamfolders):
        """Scan a list of folders containing alignment files and return a dictionary with valid files per sample name.

        Iterates over a list of paths to folders containing BAM/CRAM files. Non-existing folders are skipped. For each
        alignment file in a valid folder, the sample name is extracted.and the file is saved under its sample name in a
        dictionary.

        Parameters
        ----------
        bamfolders : list
            List of String paths to folders containing alignment files

        Returns
        -------
        bam_sample_map : dict
            String paths to alignment files per sample name
        """
        self.vaselogger.info("Start scanning BAM files")

        # Scan BAM files in all provided folders.
        for bamfolder in bamfolders:
            if not os.path.isdir(bamfolder):
                self.vaselogger.info(f"Folder {bamfolder} does not exist.")
                continue
            self.vaselogger.info(f"Scanning BAM files in {bamfolder}")
            for bamfilename in os.listdir(bamfolder):
                if bamfilename.endswith(".bam"):
                    try:
                        bamfile = pysam.AlignmentFile(
                            f"{bamfolder}/{bamfilename}",
                            "rb"
                            )
                        self.vaselogger.debug(
                            "Scanning BAM file "
                            f"{bamfolder}/{bamfilename}"
                            )
                        if self.bam_has_sample_name(bamfile):
                            # The sample identifier.
                            sampleid = bamfile.header["RG"][0]["SM"]
                            self.bam_sample_map[sampleid] = (
                                f"{bamfolder}/{bamfilename}"
                                )
                        bamfile.close()
                    except IOError:
                        self.vaselogger.warning(
                            f"BAM file {bamfolder}/{bamfilename} "
                            "does not exist"
                            )
        self.vaselogger.info("Finished scanning BAM files")
        return self.bam_sample_map

    # Scans the provided BAM/CRAM files from a provided donor list file
    def scan_bamcram_files(self, bamcramlistfile):
        """Read a list file containing the locations of BAM/CRAM files and saves them per sample.

        Parameters
        ----------
        bamcramlistfile : str
            String path to list file with paths to alignment files

        Returns
        -------
        bam_sample_map : dict
            Paths to alignment files per sample name
        """
        self.bam_sample_map = self.read_donorlistfile(bamcramlistfile)
        return self.bam_sample_map

    # Scan the provided VCF files from a provided donor list file
    def scan_vcf_files(self, vcflistfile):
        """Read a list file containing the locations of VCF/BCF files and saves them per sample name.

        Parameters
        ----------
        vcflistfile : str
            String path to list file with paths to variant files

        Returns
        -------
        vcf_sample_map : dict
            Variant files per sample name
        """
        self.vcf_sample_map = self.read_donorlistfile(vcflistfile, "v")
        return self.vcf_sample_map

    def read_donorlistfile(self, listfileloc, listtype="a"):
        """Read a list file containing paths to donor files and return the donor files per sample name.

        Parameters
        ----------
        listfileloc : str
            Path to the list file containing paths to donor files
        listtype : str
            Type of variant files in the list file (a=alignment ; v=variant)

        Returns
        -------
        donor_sample_files : dict
            Paths to donor variant files per sample name
        """
        donor_sample_files = {}
        try:
            with open(listfileloc, "r") as listfile:
                for fileline in listfile:
                    if os.path.isfile(fileline.strip()):
                        sampleid = self.get_sample_id(fileline.strip(), listtype)

                        if sampleid is not None:
                            if isinstance(sampleid, list):
                                for sid in sampleid:
                                    donor_sample_files[sid] = fileline.strip()
                            else:
                                donor_sample_files[sampleid] = fileline.strip()
        except IOError:
            self.vaselogger.critical(f"Could not open donor list file {listfileloc}")
            sys.exit()
        return donor_sample_files

    @staticmethod
    def bam_has_sample_name(bamfile):
        """Check and return whether a BAM/CRAM file has a sample name/identifier.

        Parameters
        ----------
        bamfile : pysam AlignmentFile
            Already opened pysam AlignmentFile

        Returns
        -------
        bool
            True if file has a sample name, otherwise False
        """
        if "RG" in bamfile.header:
            if len(bamfile.header["RG"]) > 0:
                if "SM" in bamfile.header["RG"][0]:
                    return True
        return False

    @staticmethod
    def vcf_has_sample_name(vcffile):
        """Check and return whether the VCF/BCF file has a sample name/identifier.

        Returns
        -------
        bool
            True if the variant file has sample names, otherwise False
        """
        if vcffile.header.samples[0] != "":
            return True
        return False

    # Returns the first sample identifier of the VCF file (as for now we only work with one sample per file
    def get_vcf_sample_name(self, vcffileloc):
        """Extract and return the list of sample names from a specified VCF file.

        Parameters
        ----------
        vcffileloc : str
            Path to a VCF file

        Returns
        -------
        list or None
            VCF Sample name if present, otherwise None
        """
        vcffile = pysam.VariantFile(vcffileloc, "r")
        if self.vcf_has_sample_name(vcffile):
            return list(vcffile.header.samples)
        return None

    def get_bam_sample_name(self, bamfileloc):
        """Extract and return the sample name/identifier from a BAM file.

        The sample name is extracted by first checking whether the BAM file contains a sample name via the
        VcfBamScanner method bam_has_sample_name(). If so, the sample name is extracted from the BAM file by returning
        the value of the 'SM' field from the header line starting with '@RG'. None is returned if the BAM file has no
        sample name.

        Parameters
        ----------
        bamfileloc : str
            Path to a BAM file

        Returns
        -------
        str or None
            Sample name if present, otherwise None
        """
        try:
            bamfile = pysam.AlignmentFile(bamfileloc, "rb")
            if self.bam_has_sample_name(bamfile):
                return bamfile.header["RG"][0]["SM"]
            self.vaselogger.debug(f"BAM file {bamfileloc} has no sample identifier")
            return None
        except IOError:
            self.vaselogger.warning(f"Could not open {bamfileloc} to extract sample identifier")

    def get_cram_sample_name(self, cramfileloc):
        """Extract and return the sample name/identifier from a specified CRAM file.

        The sample name is extracted by first checking whether the CRAM file contains a sample name via the
        VcfBamScanner method cram_has_sample_name(). If so, the sample name is extracted from the CRAM file by
        returning the value of the 'SM' field from the header line starting with '@RG'. None is returned if the CRAM
        file has no sample name.

        Parameters
        ----------
        cramfileloc : str
            Path to a CRAM file

        Returns
        -------
        str or None
            Sample name if present, otherwise None
        """
        try:
            cramfile = pysam.AlignmentFile(cramfileloc, "rc")
            if self.bam_has_sample_name(cramfile):
                return cramfile.header["RG"][0]["SM"]
            return None
        except IOError:
            self.vaselogger.warning(f"Could not open {cramfileloc} to extract sample identifier")

    # General methods that returns one or more sample ids (in case of a VCF)
    def get_sample_id(self, donorfileloc, donorlisttype):
        """Return the sample name/identifier for a specified alignment or variant file.

        Parameters
        ----------
        donorfileloc : str
            Path to a donor variant/alignment file
        donorlisttype : str
            Type of donor file (a=alignment ; v=variant)

        Returns
        -------
        str
            File type of provided donor file
        """
        donor_file_type = self.get_donor_file_type(donorfileloc)
        if donorlisttype == "a":
            if donor_file_type[1] == "BAM":
                return self.get_bam_sample_name(donorfileloc)
            if donor_file_type[0] == "CRAM":
                return self.get_cram_sample_name(donorfileloc)
        if donorlisttype == "v":
            if donor_file_type[3] == "(VCF)" or donor_file_type[3] == "(BCF)":
                return self.get_vcf_sample_name(donorfileloc)

    # Returns the filetype of a donor file
    @staticmethod
    def get_donor_file_type(donorfileloc):
        """Determine and return the file type of a specified donor file.

        To determine the file type, the linux command 'file' is used which returns a String containing the file type
        info. The received String is split on spaces and the resulting list is returned.

        Returns
        -------
        list of str
            List with file type data
        """
        filetype_proc = subprocess.Popen(["file", "-b", "-z", donorfileloc], stdout=subprocess.PIPE)
        filetype_data, filetype_err = filetype_proc.communicate()
        return filetype_data.decode().split(" ")

    # Returns a list of sample identifiers that have both a VCF and a BAM/CRAM
    def get_complete_sample_ids(self):
        """Determine which samples have a variant and alignment file and return the list of identifiers.

        Returns
        -------
        set
            Set of sample names with a variant and alignment file
        """
        return set(self.vcf_sample_map.keys()) & set(self.bam_sample_map.keys())

    # Extracts the sequence names from an already opened BAM/CRAM file
    @staticmethod
    def get_alignment_sequence_names(alignment_file):
        """Extract and return the chromosome names from an alignment file.

        Chromosome names are extracted from a specified BAM or CRAM file by checking the headers for the 'SQ' tag. If
        present the SQ entry is checked for the 'SN' field which contains the specific chromosome names.

        Parameters
        ----------
        alignment_file : pysam.AlignmentFile
            Already opened pysam AlignmentFile

        Returns
        -------
        sequence_names : list
            List with extracted chromosome names
        """
        sequence_names = set()
        if "SQ" in alignment_file.header:
            for sn_entry in alignment_file.header["SQ"]:
                sequence_names.add(sn_entry["SN"])
        return sequence_names

    @staticmethod
    def get_variant_sequence_names(variant_file):
        """Extract and return the chromosome names from a variant file.

        Chromosome names are extracted from a specified VCF or BCF file by obtaining the header contig names via the
        pysam method .header.contigs(). The provided variant file is expected to be an already opened pysam
        VariantFile.

        Parameters
        ----------
        variant_file : pysam.VariantFile
            Already opened pysam VariantFile

        Returns
        -------
        sequence_names : list of str
            List with extracted chromosome names
        """
        sequence_names = set()
        if len(variant_file.header.contigs) > 0:
            for seqname in list(variant_file.header.contigs):
                sequence_names.add(seqname)
        return sequence_names