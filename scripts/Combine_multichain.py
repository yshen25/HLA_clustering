#!usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def rename(full_name:str) -> str:
    return full_name.replace('*', '').replace(':','_')

def filename(A_name, B_name):
    A_new = rename(A_name)
    B_new = rename(B_name)

    return f"{A_new}_{B_new}.faa"

def A_extract(A_fasta, A_allele_file):
    """
    Read in one fasta file and corresponding gene names
    write each gene into seperate ffasta file
    """
    A_record = SeqIO.to_dict(SeqIO.parse(A_fasta, 'fasta'))

    with open(A_allele_file) as fh:
        A_list = [line.strip() for line in fh]

    for A_allele in A_list:
        if A_allele in A_record:
            A_rec = A_record[A_allele]
            SeqIO.write(A_rec, rename(A_allele)+'.faa', 'fasta')
        else:
            print(A_allele)

    return

def AB_combine(A_fasta, A_allele_file, B_fasta, B_allele_file):
    """
    Read in two fasta files containing A chain sequences and B chain sequences,
    combining both chains into single fasta files according to A and B allele files
    """
    A_record = SeqIO.to_dict(SeqIO.parse(A_fasta, 'fasta'))
    B_record = SeqIO.to_dict(SeqIO.parse(B_fasta, 'fasta'))

    with open(A_allele_file) as fh:
        A_list = [line.strip() for line in fh]

    with open(B_allele_file) as fh:
        B_list = [line.strip() for line in fh]

    for A_allele in A_list:
        if not A_allele in A_record:
            continue
        for B_allele in B_list:
            A_rec = A_record[A_allele]
            B_rec = B_record[B_allele]
            new_rec = SeqRecord(Seq(str(A_rec.seq)+":"+str(B_rec.seq)), id=f"{A_rec.id}_{B_rec.id}")
            SeqIO.write(new_rec, filename(A_allele, B_allele), 'fasta')

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="combine two chains")
    parser.add_argument("A_fas")
    parser.add_argument("A_names")
    parser.add_argument("B_fas", nargs='?', default=None)
    parser.add_argument("B_names", nargs='?', default=None)
    args = parser.parse_args()
    if args.B_fas:
        AB_combine(args.A_fas, args.A_names, args.B_fas, args.B_names)
    else:
        A_extract(args.A_fas, args.A_names)
