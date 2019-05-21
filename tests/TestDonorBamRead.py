import unittest
from DonorBamRead import DonorBamRead

class TestDonorBamRead(unittest.TestCase):
    # Creates the answer variables and the DonorBamRead used for testing
    def setUp(self):
        self.readIdAnswer = 'HHKY2CCXX160108:1:2122:24160:2522'
        self.readPnAnswer = '1'
        self.readChromAnswer = '21'
        self.readPosAnswer = 9411193
        self.readLenAnswer = 151
        self.readSeqAnswer = 'AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGACTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA'
        self.readQualsAnswer = '><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>>==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>'
        self.readMapQAnswer = 40
        self.dbamRead = DonorBamRead(self.readIdAnswer, self.readPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
    
    
    
    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_getBamReadId(self):
        self.assertEqual(self.dbamRead.get_bam_read_id(), self.readIdAnswer, f"Both BAM read identifiers should have been {self.readIdAnswer}")
    
    def test_getBamReadPairNumber(self):
        self.assertEquals(self.dbamRead.get_bam_read_pair_number(), self.readPnAnswer, f"Both read pair numbers should have been {self.readPnAnswer}")
    
    def test_getBamReadChrom(self):
        self.assertEquals(self.dbamRead.get_bam_read_chrom(), self.readChromAnswer, f"Both read chromosomes should have been {self.readChromAnswer}")
    
    def test_getBamReadRefPos(self):
        self.assertEqual(self.dbamRead.get_bam_read_ref_pos(), self.readPosAnswer, f"Both read positions should have been {self.readPosAnswer}")
    
    def test_getBamReadLength(self):
        self.assertEqual(self.dbamRead.get_bam_read_length(), self.readLenAnswer, f"Both read lengths should have been {self.readLenAnswer}")
    
    def test_getBamReadRefEnd(self):
        readEndAnswer = 9411344
        self.assertEqual(self.dbamRead.get_bam_read_ref_end(), readEndAnswer, f"Both read end positions should have been {readEndAnswer}")
    
    def test_getBamReadQual(self):
        self.assertEqual(self.dbamRead.get_bam_read_qual(), self.readQualsAnswer, f"Both read qualities should have been {self.readQualsAnswer}")
    
    def test_getBamReadQScores(self):
        qscoresAnswer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29, 30, 30, 30, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
        self.assertEqual(self.dbamRead.get_bam_read_q_scores(), qscoresAnswer, f"Both Q-scores should have been {qscoresAnswer}")
    
    def test_getBamReadMapQ(self):
        self.assertEqual(self.dbamRead.get_mapping_qual(), self.readMapQAnswer, f"Both MapQ values should have been {self.readMapQAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR THE STATISTICS METHODS====================
    def test_getAverageQScore(self):
        avgQScoreAnswer = 28.490066225165563
        self.assertEqual(self.dbamRead.get_average_qscore(), avgQScoreAnswer, f"Both average Q-scores should have been {avgQScoreAnswer}")
    
    def test_getMedianQScore(self):
        medQScoreAnswer = 28
        self.assertEqual(self.dbamRead.get_median_q_score(), medQScoreAnswer, f"Both median Q-scores should have been {medQScoreAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR CHECKING WHETHER THE READ IS R1 OR R2====================
    def test_isRead1(self):
        self.assertTrue(self.dbamRead.is_read1(), "Donor BAM read should have been read 1")
    
    def test_isRead2(self):
        self.assertFalse(self.dbamRead.is_read2(), "Donor BAM read should not have been read 2")
    
    
    
    # ====================PERFORM THE TESTS FOR RETURNING BAM READ STRING REPRESENTATIONS====================
    def test_toString(self):
        toStringAnswer = f"{self.readIdAnswer}\t{self.readPnAnswer}\t{self.readChromAnswer}\t{self.readPosAnswer}\t{self.readLenAnswer}\t{self.readSeqAnswer}\t{self.readQualsAnswer}\t{self.readMapQAnswer}"
        self.assertEqual(self.dbamRead.to_string(), toStringAnswer, f"Both answers should have been {toStringAnswer}")
    
    def test_getAsFastQSeq_pairnum(self):
        fastqPnAnswer = f"@{self.readIdAnswer}/{self.readPnAnswer}\n{self.readSeqAnswer}\n+\n{self.readQualsAnswer}\n"
        self.assertEqual(self.dbamRead.get_as_fastq_seq(True), fastqPnAnswer, f"Both answers should have been {fastqPnAnswer}")
    
    def test_getAsFastQSeq_nopairnum(self):
        fastqAnswer = f"@{self.readIdAnswer}\n{self.readSeqAnswer}\n+\n{self.readQualsAnswer}\n"
        self.assertEqual(self.dbamRead.get_as_fastq_seq(), fastqAnswer, f"Both answers should have been {fastqAnswer}")
