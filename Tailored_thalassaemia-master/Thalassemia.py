#!/usr/bin/env python
'''
Tailored thalassaemia mutation detection pipeline
******************
Rescue the multiply-reads and generation of rescued BAM files
Input: BAM files go through bwa alignment, followed by PICARD remove duplication
Output: Ready BAM file for Thalassamia mutation detection.
******************
'''
from BamOPR import *
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-bam", "--bamfile",
                    help="/Path/to/input_Bam_file.bam, example:/home/data/input.bam")
parser.add_argument(
    "-o", "--output", help="/Path/to/output_Bam_file.bam, example:/home/data/out.bam")
parser.add_argument("-ref", "--reference", default="hg19",
                    help="reference of human genome assebly, hg19(Default) or GRCh38")
args = parser.parse_args()

if __name__ == "__main__":

    path2bam = args.bamfile
    path2outputfile = args.output
    assembly_version = args.reference

    if assembly_version == "hg19":
        chr16_pos1 = 200000
        chr16_pos2 = 230000
        alpha_1_peak_pos1 = 226000
        alpha_1_peak_pos2 = 226500
        chr11_pos1 = 5269000
        chr11_pos2 = 5278000
    elif assembly_version == "GRCh38":
        chr16_pos1 = 150001
        chr16_pos2 = 180001
        alpha_1_peak_pos1 = 176001
        alpha_1_peak_pos2 = 176501
        chr11_pos1 = 5247770
        chr11_pos2 = 5256770

    Mean_InsertSize, Std_InsertSize = EstimateInsertSize(path2bam)

    Inferred_read_length = InferReadLength(path2bam)

    if InferReadLength <= 100:
        alpha_rescue_set = RescueSet_v2(
            path2bam, "chr16", chr16_pos1, chr16_pos2)
        alpha_rescue_set_3 = RescueSet_v3(
            path2bam, "chr16", chr16_pos1, chr16_pos2)
        alpha_1_peak_rescue_set = RescueSet_v2(
            path2bam, "chr16", alpha_1_peak_pos1, alpha_1_peak_pos2, MAQ_L=30, MAQ_H=60)

        beta_rescue_set = RescueSet_v2(
            path2bam, "chr11", chr11_pos1, chr11_pos2)
        beta_rescue_set_3 = RescueSet_v3(
            path2bam, "chr11", chr11_pos1, chr11_pos2)

        hemoglobin_rescue_set = alpha_rescue_set | beta_rescue_set | alpha_rescue_set_3 | beta_rescue_set_3 | alpha_1_peak_rescue_set

    elif InferReadLength > 100:
        alpha_rescue_set = RescueSet_v2(
            path2bam, "chr16", chr16_pos1, chr16_pos2)
        alpha_rescue_set_3 = RescueSet_v3(
            path2bam, "chr16", chr16_pos1, chr16_pos2)

        beta_rescue_set = RescueSet_v2(
            path2bam, "chr11", chr11_pos1, chr11_pos2)
        beta_rescue_set_3 = RescueSet_v3(
            path2bam, "chr11", chr11_pos1, chr11_pos2)

        hemoglobin_rescue_set = alpha_rescue_set | beta_rescue_set | alpha_rescue_set_3 | beta_rescue_set_3

    probability = 0.5
    RescueMARReads(path2bam, path2outputfile,
                   probability, hemoglobin_rescue_set)
else:
    print "Error! Please check the read length!"
