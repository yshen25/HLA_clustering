import os
import argparse

from parso import parse

def main(allele):
    for f in os.listdir("."):
        if f.endswith(".csv"):
            os.rename(f, f.split("_A")[0]+f"_{allele}.csv")

    return

parser = argparse.ArgumentParser()
parser.add_argument("allele_name")
arg = parser.parse_args()

main(arg.allele_name)
