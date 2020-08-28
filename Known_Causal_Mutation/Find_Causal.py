import sys
import os
import argparse
import pandas as pd
from Find_Causal import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--inputvcf")
parser.add_argument("-output", "--outputtxt")
args = parser.parse_args()


def F_Create_pseudo_vcf(inputvcf, outputvcf):
    alt_allele = [None] * 4
    with open(inputvcf, "r") as fp, open(outputvcf, "w") as fout:
        for line in fp:
            if line.startswith("##"):
                pass
            elif line.startswith("#"):
                div = line.rstrip("\n").split("\t")
                samplenamelists = line.strip('\n').split("\t")[9:]
                samplecount = len(samplenamelists)
                newline = div[0] + "\t" + div[1] + "\t" + \
                    div[3] + "\t" + div[4] + "\t" + \
                    "\t".join(samplenamelists) + "\n"
                fout.write(newline)
            else:
                div = line.rstrip("\n").split("\t")
                alt_num = len(line.strip('\n').split('\t')[4].split(','))
                chrom = div[0]
                pos = div[1]
                ref = div[3]
                alt = div[4]

                if alt_num == 1 and alt != "*":
                    newline = chrom + "\t" + pos + "\t" + ref + "\t" + alt
                    for i in range(9, samplecount + 9):
                        raw_GT = div[i].split(":")[0]
                        if raw_GT == "0/0" or raw_GT == "0|0":
                            new_GT = ref + "_" + ref
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "0/1" or raw_GT == "0|1":
                            new_GT = ref + "_" + alt
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "1/1" or raw_GT == "1|1":
                            new_GT = alt + "_" + alt
                            newline = newline + "\t" + new_GT
                        else:
                            new_GT = "Untyped"
                    newline = newline + "\n"
                    fout.write(newline)
                else:
                    newline = chrom + "\t" + pos + "\t" + ref + "\t" + alt
                    alt_allele[1] = alt.split(",")[0]
                    try:
                        alt_allele[2] = alt.split(",")[1]
                    except IndexError:
                        alt_allele[2] = "NA"
                    try:
                        alt_allele[3] = alt.split(",")[2]
                    except IndexError:
                        alt_allele[3] = "NA"

                    for i in range(9, samplecount + 9):
                        raw_GT = div[i].split(":")[0]
                        if raw_GT == "0/0" or raw_GT == "0|0":
                            new_GT = ref + "_" + ref
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "0/1" or raw_GT == "0|1":
                            new_GT = ref + "_" + alt_allele[1]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "0/2" or raw_GT == "0|2":
                            new_GT = ref + "_" + alt_allele[2]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "0/3" or raw_GT == "0|3":
                            new_GT = ref + "_" + alt_allele[3]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "1/1" or raw_GT == "1|1":
                            new_GT = alt_allele[1] + "_" + alt_allele[1]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "1/2" or raw_GT == "1|2":
                            new_GT = alt_allele[1] + "_" + alt_allele[2]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "1/3" or raw_GT == "1|3":
                            new_GT = alt_allele[1] + "_" + alt_allele[3]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "2/2" or raw_GT == "2|2":
                            new_GT = alt_allele[2] + "_" + alt_allele[2]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "2/3" or raw_GT == "2|3":
                            new_GT = alt_allele[2] + "_" + alt_allele[3]
                            newline = newline + "\t" + new_GT
                        elif raw_GT == "3/3" or raw_GT == "3|3":
                            new_GT = alt_allele[3] + "_" + alt_allele[3]
                            newline = newline + "\t" + new_GT
                        else:
                            new_GT = "Untyped"
                    newline = newline + "\n"
                    fout.write(newline)
    return 1


def F_split_pseudovcf_by_sample(vcffile, outputfolder):
    GT_matrix = pd.read_csv(vcffile, sep="\t", header=0)
    with open(vcffile, "r") as fp:
        line = fp.readline()
        div = line.rstrip("\n").split("\t")
        samplenamelists = line.strip('\n').split("\t")[4:]
        samplecount = len(samplenamelists)
    for sample in samplenamelists:
        outputfilename = str(sample) + ".txt"
        filepath = outputfolder + "/" + outputfilename
        df = GT_matrix.loc[:, ['#CHROM', 'POS', sample]]
        df.to_csv(filepath, index=False, sep="\t")
    return samplenamelists, samplecount


def F_match_causal(samplefile, knowncausal, outputfile):
    db_dict = {}
    with open(knowncausal, "r") as fp:
        for line in fp:
            if not line.startswith("#"):
                div = line.rstrip("\n").split("\t")
                key = div[0] + "_" + str(div[1]) + "_" + div[3] + "_" + div[4]
                value = div[7]
                db_dict[key] = value
    with open(samplefile, "r") as fp, open(outputfile, "w") as fout:
        for line in fp:
            if line.startswith("#"):
                fout.write(line)
            else:
                div = line.rstrip("\n").split("\t")
                current_key = div[0] + "_" + str(div[1]) + "_" + div[2]
                if current_key in db_dict.keys():
                    newline = line.rstrip("\n") + "\t" + \
                        db_dict[current_key] + "\n"
                    fout.write(newline)
    return 1


if __name__ == "__main__":
    # SNP
    code_F = os.getcwd()
    inputvcf = args.inputvcf
    wkd = code_F.replace("VCF_file/Joint", "")

    pseudo_vcf_file = wkd + "VCF_file/Joint/Pseudo_Candidate_SNP.recode.vcf"
    F_Create_pseudo_vcf(inputvcf=inputvcf, outputvcf=pseudo_vcf_file)
    Causal_SNV_F = wkd + "VCF_file/Joint/ind_vcf_SNP"
    os.mkdir(Causal_SNV_F)
    samplenamelists = []
    samplecount = []
    samplenamelists, samplecount = F_split_pseudovcf_by_sample(
        vcffile=pseudo_vcf_file, outputfolder=Causal_SNV_F)

    knowncausalSNV = wkd + \
        "Thala_Rescue_workflow/Known_Causal_Mutation/sorted_Causal_SNV_Thala_with_equivalent.vcf"
    for i in samplenamelists:
        samplefile = Causal_SNV_F + "/" + i + ".txt"
        outputfile = Causal_SNV_F + "/" + "pre." + i
        F_match_causal(samplefile=samplefile,
                       knowncausal=knowncausalSNV, outputfile=outputfile)
# indel
    pseudo_vcf_file = wkd + "VCF_file/Joint/Pseudo_Candidate_INDEL.recode.vcf"
    F_Create_pseudo_vcf(inputvcf=inputvcf, outputvcf=pseudo_vcf_file)
    Causal_INDEL_F = wkd + "VCF_file/Joint/ind_vcf_INDEL"
    os.mkdir(Causal_INDEL_F)
    samplenamelists = []
    samplecount = []
    samplenamelists, samplecount = F_split_pseudovcf_by_sample(
        vcffile=pseudo_vcf_file, outputfolder=Causal_INDEL_F)

    knowncausalINDEL = wkd + \
        "Thala_Rescue_workflow/Known_Causal_Mutation/sorted_normed_Causal_Indel_Thala_with_equivalent.vcf"
    for i in samplenamelists:
        samplefile = Causal_INDEL_F + "/" + i + ".txt"
        outputfile = Causal_INDEL_F + "/" + "pre." + i
        F_match_causal(samplefile=samplefile,
                       knowncausal=knowncausalINDEL, outputfile=outputfile)
else:
    pass
