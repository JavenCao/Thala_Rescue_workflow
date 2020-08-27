#!/usr/bin/env python
'''
Module for NGS4THAL - a one-stop molecular diagnosis and screening tool for thalassaemia
******************
Rescue the multiply-reads and generation of rescued BAM files
Input: bwa-mem BAM files
Output: Re-aligned BAM file ready for SNV and InDels callings
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
        chr11_pos1 = 5269000
        chr11_pos2 = 5278000
    elif assembly_version == "GRCh38":
        chr16_pos1 = 150001
        chr16_pos2 = 180001
        chr11_pos1 = 5247770
        chr11_pos2 = 5256770

    alpha_rescue_set = RescueSet(
        path2bam, "chr16", chr16_pos1, chr16_pos2)

    beta_rescue_set = RescueSet(
        path2bam, "chr11", chr11_pos1, chr11_pos2)

    hemoglobin_rescue_set = alpha_rescue_set | beta_rescue_set

    probability = 0.5
    RescueMARReads(path2bam, path2outputfile,
                   probability, hemoglobin_rescue_set)
else:
    pass
