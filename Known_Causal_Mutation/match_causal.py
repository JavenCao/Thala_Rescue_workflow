from match_causal import *
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-samplefile", "--samplefile")
parser.add_argument("-known", "--knownmutation")
parser.add_argument("-of", "--outputfile")
args = parser.parse_args()


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
    samplefile = args.samplefile
    knowncausal = args.knownmutation
    outputfile = args.outputfile
    F_match_causal(samplefile=samplefile,
                   knowncausal=knowncausal, outputfile=outputfile)
else:
    pass
