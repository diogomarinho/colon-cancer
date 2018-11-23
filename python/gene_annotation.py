#!/usr/bin/python

import pandas as pd
import numpy as np
import argparse
import sys

'''
This script will generate a list of CpGs that belongs to a gene (one per row)
with a given windows size and a gene annotation file
'''
def cpg_annotation(illumina, gcoord, w, output):
    print('creating CpG list per gene')
    with open(output, 'w') as fout:
        for i, row in gcoord.iterrows():
            start = row.START - w
            if (start < 0):
                start += w
            end = row.END + w
            illumina_aux = illumina[illumina.CHR == row.CHR]
            downstream = illumina_aux.POS >= start
            upstream = illumina_aux.POS <= end
            which = (np.argwhere((downstream & upstream).values)).flatten()
            cpgs = illumina_aux.NAME.values[which].tolist()
            if (len(cpgs) > 2):
                str_format = row.GENE + ' ' + ' '.join(cpgs) + '\n'
                fout.write(str_format)


def create_baseline(prefix):
    print('Merge all annot data....')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate SNP annotations giving a sub set of genes')
    parser.add_argument('-c','--coordinate', action='store', required=True,
                        help='genes coordinate file')
    parser.add_argument('-i', '--illumina', action='store', required=True,
                     help='prefix of illumina file')
    parser.add_argument('-o', '--output', action='store', required=True,
                     help='prefix for output')
    parser.add_argument('-w', '--window', action='store', default='100000')
	# .
    args = parser.parse_args()

    # GENE    CHR START   END STRAND  HGNC
    try:
        gcoord = pd.read_csv(
            args.coordinate,
            dtype = {
                'ID':np.str,
                'start':np.int,
                'end':np.int,
                'chr':np.str
            }
        )
        gcoord.columns = 'GENE START END CHR'.split()
        gcoord = gcoord.dropna(axis=0)
        gcoord = gcoord.drop_duplicates()
        gcoord = gcoord.reset_index(drop=True)
    except Exception:
        sys.exit(0)

    try:
        illumina = pd.read_csv(
            args.illumina,
            usecols = [
                'chr',
                'pos',
                'Name'
            ]
        )
        illumina.columns = 'CHR POS NAME'.split()
    except Exception:
        sys.exit(0)
    try:
        w = int(args.window)
    except:
        w=100000
        print('window size invalid setting to 100000')
    # .
    output =  args.output
    chrom = str(illumina.iloc[0,0])

    # removing X Y chromosomes
    to_keep = ~illumina.CHR.isin('chrX chrY'.split()).values.reshape(-1)
    illumina = illumina.loc[to_keep, :]
    to_keep = ~gcoord.CHR.isin('X Y'.split()).values.reshape(-1)
    gcoord = gcoord.loc[to_keep, :]
    # . replacing chrA to A
    for i in range(1, 23):
        illumina.CHR = illumina.CHR.replace('chr{X}'.format(X=i), str(i))
    #
    annotation_output = cpg_annotation(illumina, gcoord, w, args.output)
