#!/usr/bin/env python3

import os
import re
import argparse
import hashlib
from Bio import SeqIO

def main():
    input_parser = argparse.ArgumentParser(description='gb2db: the script to parse GenBank files to load them in a database.')
    input_parser.add_argument('-f', metavar='input_file.gb', help='GenBank input file', required=True)
    input_parser.add_argument('-o', metavar='/path/to/output_dir', default='output', help='the path to the output directory', required=False)
    input_parser.add_argument('-c', action='store_true', help='extract only complete genomes', required=False)

    args = input_parser.parse_args()
    input_file = args.f
    comlete = args.c
    out_dir = args.o

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    file_desc = out_dir + '/description.csv'
    file_feat = out_dir + '/features.csv'

    with open(file_desc, 'w') as file1, open(file_feat, 'w') as file2:
        file1.write("\t".join(['id', 'name', 'source', 'organism', 'description', 'seq', 'seq_md5']) + "\n")
        file2.write("\t".join(['id', 'feature_number', 'feature_type', 'qualifier', 'value']) + "\n")
        for rec in SeqIO.parse(input_file, 'genbank'):
            if comlete and not re.search(r'complete.+genome', rec.description):
                continue
            name = rec.name
            id = rec.id
            seq = rec.seq
            seq_md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
            desc = rec.description
            source = ''
            organism = ''
            if 'source' in rec.annotations:
                source = rec.annotations["source"]
            if 'organism' in rec.annotations:
                organism = rec.annotations["organism"]
            file1.write("\t".join([id, name, source, organism, desc, str(seq), seq_md5]) + "\n")

            n_feat = 0
            for feat in rec.features:
                if feat.type == "CDS" or feat.type == "gene":
                    n_feat += 1
                    loc = feat.location
                    for qual in feat.qualifiers:
                        for item in feat.qualifiers[qual]:
                            tran_md5 = 0
                            if qual == 'translation':
                                tran_md5 = hashlib.md5(item.encode('utf-8')).hexdigest()
                                file2.write("\t".join([id, str(n_feat), feat.type, 'translation_md5sum', tran_md5]) + "\n")
                            file2.write("\t".join([id, str(n_feat), feat.type, qual, item]) + "\n")
# end of main()

if __name__ == '__main__':
    main()
    print("ok")
