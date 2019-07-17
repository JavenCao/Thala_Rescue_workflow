# update on July 8, 2019
# This file includes all supporting funcions
# June 1, 2018
# Cao Yujie(boxyjcao@gmail.com)
"""
Supporting functions for Thalassemia.py, which is a tailored thalassaemia mutation detection pipeline
**************
supporting function lists:
def EstimateInsertSize()
def RescueSet()
def RescueSet_v2()
def RescueSet_v3()
def RescueMARReads()
def ExtractUnmappedSingleEndReads()
def ExtractUnmapped()
def GivenNameSetGetSeq()
def ChangeMAQOnly()
def InferReadLength()
**************
"""
import pysam
import random
import numpy as np
from BamOPR import *


def EstimateInsertSize(path2bam, chromosome="chr16", start_pos=135000, stop_pos=200000, MAQ=60):
    """
    Description: Use high quality alignments to estimate insert size mean and SD.
    Input:
         path2bam: path to alignment file(.bam)
         chromosome(string), start_pos(int),stop_pos(int): explicitly give a region to use, default is hg19 chr16:135000-200000
         MAQ(int): mapping quality cut-off, default is 60
    Return:
         a turple contains (Mean_InsertSize(float), Std_InsertSize(float))
    """

    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(chromosome, start=start_pos, end=stop_pos)
    if samfile.count(contig=chromosome, start=start_pos, end=stop_pos) < 100:
        return int(300), int(80)
    else:
        isize = []

        for x in chr_sam:
            if x.is_paired and x.is_proper_pair and x.is_read1 and x.mapping_quality >= MAQ and abs(x.tlen) < 10000:
                isize.append(abs(x.tlen))

        Mean_InsertSize = np.mean(isize)
        Std_InsertSize = np.std(isize)
    return int(Mean_InsertSize), int(Std_InsertSize)


def RescueSet(path2bam, chromosome, start_pos, stop_pos):
    """
    Description: Return a set of reads name that needs to be rescued, becasue of mutliple alignment
    Input:
         path2bam(string): path to bwa mapping file(.bam)
    Return:
         rescue_set: a set of reads name that needs to be rescued, becasue of mutliple alignment
    Here, we define mulitply reads as (1)MAQ=0; (2) 0<MAQ<30 and XA tag and NM <= 3
    """
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(chromosome, start=start_pos, end=stop_pos)
    rescue_set = set()
    for x in chr_sam:
        if x.mapping_quality >= 0 and x.mapping_quality < 50 and x.is_paired and x.is_proper_pair and x.is_read1:
            if x.mapping_quality == 0 and x.has_tag("XA"):
                rescue_set.add(x.query_name)
            elif x.mapping_quality > 0 and x.mapping_quality < 50 and x.has_tag("XA"):
                try:
                    NM_tag_line = x.get_tag("NM")
                    NM_NM = int(NM_tag_line)
                    XA_tag_line = x.get_tag("XA")
                    XA_NM = int(XA_tag_line.split(";")[0].split(",")[3])
                    if abs(NM_NM - XA_NM) <= 3:
                        rescue_set.add(x.query_name)
                    else:
                        continue
                except KeyError:
                    continue
            else:
                continue
    return rescue_set


def RescueSet_v2(path2bam, chromosome, start_pos, stop_pos, MAQ_L=0, MAQ_H=40):
    """
    Description: Return a set of reads name that needs to be rescued, becasue both the PE reads are with mutliple alignment
    Input:
         path2bam(string): path to bwa mapping file(.bam)
    Return:
         rescue_set: a set of reads name that needs to be rescued, becasue of mutliple alignment
    Here, we define mulitply reads as (1) 0<MAQ<=30; (2)R1 and R2 both have XA_tag showing the alternative tag on the genome
    """
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(chromosome, start=start_pos, end=stop_pos)
    rescue_set = set()
    PE_flag_dict = {}
    for x in chr_sam:
        if x.mapping_quality >= MAQ_L and x.mapping_quality <= MAQ_H and x.is_paired and (not x.is_secondary) and x.has_tag("XA"):
            if x.query_name in PE_flag_dict.keys():
                if x.is_read1:
                    PE_flag_dict[x.query_name] = PE_flag_dict[x.query_name] + "R1"
                if x.is_read2:
                    PE_flag_dict[x.query_name] = PE_flag_dict[x.query_name] + "R2"
            else:
                if x.is_read1:
                    PE_flag_dict[x.query_name] = "R1"
                if x.is_read2:
                    PE_flag_dict[x.query_name] = "R2"
        else:
            continue
    keys_list = PE_flag_dict.keys()
    for name in keys_list:
        if PE_flag_dict[name] == "R1R2" or PE_flag_dict[name] == "R2R1":
            rescue_set.add(name)
    return rescue_set


def RescueSet_v3(path2bam, chromosome, start_pos, stop_pos, MAQ_L=0, MAQ_H=60):
    """
    Description: Return a set of reads name that needs to be rescued, becasue one read is with high quality, the other is low quality
    Input:
         path2bam(string): path to bwa mapping file(.bam)
    Return:
         rescue_set: a set of reads name that needs to be rescued, becasue of mutliple alignment
    Here, we define mulitply reads as (1) 0<MAQ<=30; (2)R1 and R2 both have XA_tag showing the alternative tag on the genome
    """
    Mean_InsertSize, Std_InsertSize = EstimateInsertSize(path2bam)
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(chromosome, start=start_pos, end=stop_pos)
    rescue_set = set()
    PE_flag_dict = {}
    for x in chr_sam:
        if x.mapping_quality >= MAQ_L and x.mapping_quality <= MAQ_H and x.is_paired and (not x.is_secondary) and x.tlen >= Mean_InsertSize + 5 * Std_InsertSize:
            if x.query_name in PE_flag_dict.keys():
                if x.is_read1:
                    PE_flag_dict[x.query_name] = PE_flag_dict[x.query_name] + "R1"
                if x.is_read2:
                    PE_flag_dict[x.query_name] = PE_flag_dict[x.query_name] + "R2"
            else:
                if x.is_read1:
                    PE_flag_dict[x.query_name] = "R1"
                if x.is_read2:
                    PE_flag_dict[x.query_name] = "R2"
        else:
            continue
    keys_list = PE_flag_dict.keys()
    for name in keys_list:
        if PE_flag_dict[name] == "R1R2" or PE_flag_dict[name] == "R2R1":
            rescue_set.add(name)
    return rescue_set


def RescueMARReads(path2bam, path2newbam, probability, rescue_set):
    """
    Description: Rescue reads in MAR(multiple alignment region)
    Dependency: def EstimateInsertSize
    Input:
         path2bam: /path/to/bwa/mapping/file(.bam)
         path2newbam: /path/to/rescued/bam
         probability: assign probability
         rescue_set: a set of reads name that needs to be rescued(generated by RescueSet_v2 and RescueSet_v3)
    Return:
          1
    """
    Mean_InsertSize, Std_InsertSize = EstimateInsertSize(path2bam)
    samfile = pysam.AlignmentFile(path2bam, "rb")
# round 1 -----------------------rescue R1----------------------------------
    # rescue R1 now, and write to "R1_mid_outputsamfile"
    R1_mid_outputsamfile = path2newbam + ".R1.bam"
    fout = pysam.AlignmentFile(R1_mid_outputsamfile, 'wb', template=samfile)
    # scan all alignments
    chr_sam = samfile.fetch()
    # store R1 reads position for R2 use in the next round
    R1_pos_dict = {}
    R1_strand_dict = {}

    for x in chr_sam:
        if x.query_name in rescue_set and x.is_read1 and x.is_paired and (not x.is_secondary):
            try:
                if x.is_reverse:
                    strand_flag = "-"
                else:
                    strand_flag = "+"

                XA_tag_line = x.get_tag("XA")
                XA_chr = XA_tag_line.split(";")[0].split(",")[0]
                XA_chr_id = samfile.get_tid(XA_chr)
                # for hemoglobin region, ++, --
                XA_pos = abs(int(XA_tag_line.split(";")[0].split(",")[1])) - 1
                XA_cigar = XA_tag_line.split(";")[0].split(",")[2]
                XA_NM = int(XA_tag_line.split(";")[0].split(",")[3])

                rand = random.uniform(0, 1)
                if rand > probability:  # without changing the position. how about XA tag???? need to be tested
                    # needs to be considered again
                    x.mapping_quality = 60
                    x.next_reference_id = x.reference_id
                    fout.write(x)
                    R1_pos_dict[x.query_name] = x.reference_name + \
                        "_" + str(x.pos)  # chr16_215397
                    R1_strand_dict[x.query_name] = strand_flag
                else:  # needs to change the position and related parameters, such as POS, MAQ, CIGAR et al
                    x.reference_id = XA_chr_id
                    x.pos = XA_pos
                    x.mapping_quality = 60
                    x.cigarstring = XA_cigar
                    x.next_reference_id = x.reference_id
                    fout.write(x)
                    R1_pos_dict[x.query_name] = x.reference_name + \
                        "_" + str(x.pos)  # positive number in hg19(+)
                    R1_strand_dict[x.query_name] = strand_flag
            except KeyError:
                fout.write(x)
                continue
        else:
            fout.write(x)
    samfile.close()
    fout.close()
    #sort and index
    sort_R1 = R1_mid_outputsamfile + "_sort.bam"
    pysam.sort("-o", sort_R1, R1_mid_outputsamfile)
    pysam.index(sort_R1)
# round 2-----------------------rescue R2 and update PNEXT field of R1----------------------------------
    samfile = pysam.AlignmentFile(sort_R1, "rb")
    R2_mid_outputsamfile = sort_R1 + ".R2.bam"
    fout = pysam.AlignmentFile(R2_mid_outputsamfile, "wb", template=samfile)
    chr_sam_round2 = samfile.fetch()

    R2_pos_dict = {}
    R2_strand_dict = {}
    for x in chr_sam_round2:
        if x.query_name in rescue_set and x.is_read2 and x.is_paired and (not x.is_secondary):
            try:
                if x.is_reverse:
                    strand_flag = "-"
                else:
                    strand_flag = "+"

                XA_tag_line = x.get_tag("XA")
                XA_chr = XA_tag_line.split(";")[0].split(",")[0]
                XA_chr_id = samfile.get_tid(XA_chr)
                XA_pos = abs(int(XA_tag_line.split(";")[0].split(",")[1])) - 1
                XA_cigar = XA_tag_line.split(";")[0].split(",")[2]
                XA_NM = int(XA_tag_line.split(";")[0].split(",")[3])

                R1_pos = int(R1_pos_dict[x.query_name].split(
                    "_")[1])  # chr16_215397
                R2_lower = R1_pos + Mean_InsertSize - x.query_length - 10 * Std_InsertSize
                R2_upper = R1_pos + Mean_InsertSize - x.query_length + 10 * Std_InsertSize

                if x.pos >= R2_lower and x.pos <= R2_upper and (XA_pos < R2_lower or XA_pos > R2_upper):
                    x.mapping_quality = 60
                    x.next_reference_id = x.reference_id
                    x.next_reference_start = R1_pos
                    R2_pos_dict[x.query_name] = x.reference_name + \
                        "_" + str(x.pos)  # chr16_215397
                    R2_strand_dict[x.query_name] = strand_flag
                    # calculate tlen with sign
                    R2_pos = x.pos
                    absTlen = abs(R2_pos - R1_pos + x.query_length)
                    tlen = R2_strand_dict[x.query_name] + str(absTlen)
                    x.tlen = int(tlen)
                    fout.write(x)
                elif XA_pos >= R2_lower and XA_pos <= R2_upper and (x.pos < R2_lower or x.pos > R2_upper):
                    x.reference_id = XA_chr_id
                    x.pos = XA_pos
                    x.mapping_quality = 60
                    x.cigarstring = XA_cigar
                    x.next_reference_id = x.reference_id
                    x.next_reference_start = R1_pos
                    R2_pos_dict[x.query_name] = x.reference_name + \
                        "_" + str(x.pos)  # chr16_215397
                    R2_strand_dict[x.query_name] = strand_flag
                    # calculate tlen with sign
                    R2_pos = x.pos
                    absTlen = abs(R2_pos - R1_pos + x.query_length)
                    tlen = R2_strand_dict[x.query_name] + str(absTlen)
                    x.tlen = int(tlen)
                    fout.write(x)
                elif x.pos >= R2_lower and x.pos <= R2_upper and XA_pos >= R2_lower and XA_pos <= R2_upper:
                    # just keep the current hit
                    x.mapping_quality = 60
                    x.next_reference_id = x.reference_id
                    x.next_reference_start = R1_pos
                    R2_pos_dict[x.query_name] = x.reference_name + \
                        "_" + str(x.pos)  # chr16_215397
                    R2_strand_dict[x.query_name] = strand_flag
                    # calculate tlen with sign
                    R2_pos = x.pos
                    absTlen = abs(R2_pos - R1_pos + x.query_length)
                    tlen = R2_strand_dict[x.query_name] + str(absTlen)
                    x.tlen = int(tlen)
                    fout.write(x)

                else:
                    # ignore these reads
                    pass
            except KeyError:
                fout.write(x)
                continue
        else:
            fout.write(x)
    samfile.close()
    fout.close()
    #sort and index
    sort_R2 = R2_mid_outputsamfile + "_sort.bam"
    pysam.sort("-o", sort_R2, R2_mid_outputsamfile)
    pysam.index(sort_R2)
# Round 3---------------------next, update R2 information of rescued R1, and create the final rescued bam--------------
    samfile = pysam.AlignmentFile(sort_R2, "rb")
    fout = pysam.AlignmentFile(path2newbam, "wb", template=samfile)
    chr_sam_round3 = samfile.fetch()
    for x in chr_sam_round3:
        if x.query_name in rescue_set and x.is_read1 and x.is_paired and (not x.is_secondary):
            try:
                R2_pos = int(R2_pos_dict[x.query_name].split("_")[1])
                x.next_reference_start = R2_pos
                # calculate tlen with sign
                R1_pos = x.pos
                absTlen = abs(R2_pos - R1_pos + x.query_length)
                tlen = R1_strand_dict[x.query_name] + str(absTlen)
                x.tlen = int(tlen)
                fout.write(x)
            except KeyError:
                fout.write(x)
                continue
        else:
            fout.write(x)
    samfile.close()
    fout.close()
    #sort and index
    sort_rescue_bam = path2newbam + "_sort.bam"
    pysam.sort("-o", sort_rescue_bam, path2newbam)
    pysam.index(sort_rescue_bam)
    return 1


def ExtractUnmapped(path2bam, chromosome, start_pos, stop_pos):
    """
    Description: Extract a set of unmapped reads name
    Input:
         path2bam: path to the bam file
    Return: a set of sequence name
    """
    unmapped_name_set = set()
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(chromosome, start=start_pos,
                            end=stop_pos, until_eof=True)
    count_R1 = 0
    count_R2 = 0
    for x in chr_sam:
        if x.is_unmapped:
            if x.is_read1:
                count_R1 = count_R1 + 1
            if x.is_read2:
                count_R2 = count_R2 + 1
            unmapped_name_set.add(x.query_name)
    return unmapped_name_set


def ExtractUnmappedSingleEndReads(path2bam):
    """
    Description: Extract a set of unmapped reads name
    Input:
         path2bam: path to the bam file

    Return: two sets of unmapped reads name
    """
    unmapped_name_set_R1 = set()
    unmapped_name_set_R2 = set()
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(until_eof=True)
    for x in chr_sam:
        if x.is_unmapped:
            if x.is_read1:
                unmapped_name_set_R1.add(x.query_name)
            if x.is_read2:
                unmapped_name_set_R2.add(x.query_name)
    return unmapped_name_set_R1, unmapped_name_set_R2


def GivenNameSetGetSeq(name_set, path2RawR1, path2RawR2, prefix2NewFQ):
    """
    Descrition: given a set of sequence name, and fastq1, fastq2 file, extract paired-end reads from them
    Input:
        name_set: a set of name whose sequence will be extracted from the original R1 and R2 fastq files
        path2RawR1: path to raw R1.fastq file
        path2RawR2: path to raw R2.fastq file
        prefix2NewFQ: the prefix for new fastq file
    output:
        prefix2NewFQ.R1.fq and prefix2NewFQ.R2.fq contains raw PE reads
    """
    New_R1 = prefix2NewFQ + ".R1.fq"
    New_R2 = prefix2NewFQ + ".R2.fq"

    with pysam.FastxFile(path2RawR1) as p_Raw_R1, open(New_R1, mode='w') as fout:
        for x in p_Raw_R1:
            if x.name in name_set:
                fout.write(str(x) + "\n")
    with pysam.FastxFile(path2RawR2) as p_Raw_R2, open(New_R2, mode='w') as fout:
        for x in p_Raw_R2:
            if x.name in name_set:
                fout.write(str(x) + "\n")

    return 1


def ChangeMAQOnly(path2bam, path2newbam, chromosome, start_pos, stop_pos):
    """
    Description: only change the Mapping Qualtiy field into 60
    """
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch(chromosome, start=start_pos, end=stop_pos)
    temp_fout = path2newbam + "_temp.bam"
    fout = pysam.AlignmentFile(temp_fout, "wb", template=samfile)
    for x in chr_sam:
        x.mapping_quality = 60
        fout.write(x)
    pysam.sort("-o", temp_fout, path2newbam)
    return 1


def InferReadLength(path2bam):
    """
    Description: infer read length from CIGAR alignment.This method deduces the read length from the CIGAR alignment including hard-clipped bases.
    """
    samfile = pysam.AlignmentFile(path2bam, "rb")
    chr_sam = samfile.fetch()

    x = next(chr_sam)
    read_length = x.infer_read_length()

    while read_length is None:
        x = next(chr_sam)
        read_length = x.infer_read_length()
    return read_length
