#!/usr/bin/env python
'''
Find causal point mutations for thalassaemia
******************
python Find_causal.py -db 5_substitution.vcf -ind VN424_HardFiltering_SNP.recode.vcf -o VN424.causal.txt

Or

python Find_causal.py -db 5_substitution.vcf --bulkvcf /home/data/bulk.vcf --bulkbamf /home/data/project/Bam_file/ -o VN424.causal.txt

******************
'''

from Find_causal import *
import argparse


def Build_dict_error(path2file):
    singledb = {}
    orsingledb = {}
    with open(path2file, 'rb') as fp:
        for line in fp:
            if not line.startswith('#'):
                div = line.rstrip('\n').split("\t")
                info_string = div[7]
                key_string = div[0] + "_" + \
                    str(div[1]) + "_" + div[3] + "_" + div[4]
                value_string = info_string
                hgvs_name_string = info_string.split(";")[2].split("=")[1]

                if ('_or_' not in hgvs_name_string) & ('semicolon' not in hgvs_name_string):

                    singledb[key_string] = value_string

                elif ('_or_' in hgvs_name_string) & (hgvs_name_string.count('_or_') == 1) & ('semicolon' not in hgvs_name_string):

                    singledb[key_string] = value_string

                    next_line = fp.readline()
                    div = line.rstrip('\n').split("\t")
                    info_string = div[7]
                    key_string = div[0] + "_" + \
                        str(div[1]) + "_" + div[3] + "_" + div[4]
                    value_string = info_string
                    hgvs_name_string = info_string.split(";")[2].split("=")[1]
                    orsingledb[key_string] = value_string
                else:
                    pass
    return singledb, orsingledb


def Build_dict(path2file):
    singledb = {}
    orsingledb = {}
    with open(path2file, 'rb') as fp:
        lines = fp.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if not line.startswith('#'):
                div = line.rstrip('\n').split("\t")
                info_string = div[7]
                key_string = div[0] + "_" + \
                    str(div[1]) + "_" + div[3] + "_" + div[4]
                value_string = info_string
                hgvs_name_string = info_string.split(";")[2].split("=")[1]

                if ('_or_' not in hgvs_name_string) & ('semicolon' not in hgvs_name_string):

                    singledb[key_string] = value_string

                elif ('_or_' in hgvs_name_string) & (hgvs_name_string.count('_or_') == 1) & ('semicolon' not in hgvs_name_string):

                    singledb[key_string] = value_string
                    if (i < len(lines)):

                        next_line = lines[i + 1]
                        div = line.rstrip('\n').split("\t")
                        info_string = div[7]
                        key_string = div[0] + "_" + \
                            str(div[1]) + "_" + div[3] + "_" + div[4]
                        value_string = info_string
                        hgvs_name_string = info_string.split(";")[
                            2].split("=")[1]
                        orsingledb[key_string] = value_string
                else:
                    pass

    return singledb, orsingledb


def String_Scan_Dict(inputline, dbdict):
    flag = 0
    if inputline.startswith('#'):
        flag = 0
        return flag, "NA"

    div = inputline.rstrip("\n").split("\t")
    equry_key_string = div[0] + "_" + str(div[1]) + "_" + div[3] + "_" + div[4]

    if equry_key_string in dbdict.keys():
        flag = 1
        return flag, dbdict[equry_key_string]
    else:
        flag = 0
        return flag, dbdict.get(equry_key_string, "NA")


def Main_Find_causal(ind_vcf_file, db_vcf_file, outputfile):

    singledb = {}
    orsingledb = {}
    singledb, orsingledb = Build_dict(db_vcf_file)
    ind_hit = {}
    hit_flag_single = 0
    hit_flag_orsingle = 0
    dbvcf_info_single = ""
    dbvcf_info_orsingle = ""
    with open(ind_vcf_file, 'rb') as fp, open(outputfile, "wb") as fout:
        for line in fp:
            if not line.startswith("#"):
                hit_flag_single = 0
                hit_flag_orsingle = 0
                hit_flag_single, dbvcf_info_single = String_Scan_Dict(
                    line, singledb)
                hit_flag_orsingle, dbvcf_info_orsingle = String_Scan_Dict(
                    line, orsingledb)

                if (hit_flag_single == 1) & (hit_flag_orsingle == 0):
                    div = line.rstrip("\n").split("\t")
                    div[7] = dbvcf_info_single
                    newline = "\t".join(div) + "\n"
                    fout.write(newline)
                    hit_flag_single = 0
                    hit_flag_orsingle = 0
                elif hit_flag_single == 0 & hit_flag_orsingle == 1:
                    div = line.rstrip("\n").split("\t")
                    div[7] = dbvcf_info_orsingle
                    newline = "\t".join(div) + "\n"
                    fout.write(newline)
                    hit_flag_single = 0
                    hit_flag_orsingle = 0
    return 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-db", "--databasevcf",
                        help="/Path/to/databasevcf, example:/home/data/known_substitution.vcf")

    parser.add_argument(
        "-bulk", "--bulkvcf", help="/Path/to/vcf files with all samples listed, example:/home/data/bulk.vcf")

    parser.add_argument(
        "-bulkbf", "--bulkbamf", help="Path to Bam folder of the project,for example:/home/data/Bam_file. And the Bam_folder should be structured as described in Step0")

    parser.add_argument(
        "-ind", "--indvcf", help="/Path/to/individual.vcf, example:/home/data/individual.vcf")

    parser.add_argument(
        "-indb", "--indbam", help="For igv snapshot, /Path/to/individual.bam, example:/home/data/individual.bam")

    parser.add_argument(
        "-o", "--output", help="/Path/to/output_causal.pre, example:/home/data/output_causal.pre")

    args = parser.parse_args()

    ind_vcf_file = args.indvcf
    db_vcf_file = args.databasevcf
    outputfile = args.output
    Main_Find_causal(ind_vcf_file, db_vcf_file, outputfile)
