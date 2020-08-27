from Create_pseudo_vcf import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-eq", "--enquiry")
parser.add_argument("-o", "--output")
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


if __name__ == "__main__":
    eq_vcf_file = args.enquiry
    output_vcf_file = args.output
    F_Create_pseudo_vcf(inputvcf=eq_vcf_file, outputvcf=output_vcf_file)
else:
    pass
