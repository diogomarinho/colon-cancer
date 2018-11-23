#!/usr/bin/python

#decompose gene expression data in PCS and write a CSV file
#to import this file inside the interactive mode type:
#  exec(open("./single_cell.py").read())

import sys
import pandas as pd
import numpy as np
from sklearn import decomposition
from sklearn.manifold import TSNE
#from ggplot import *
import argparse
from sklearn.decomposition import PCA

#import urllib
#import urllib2

class tsne_gene:
    def __init__(self, pcs, geneid):
        self.pcs = pcs
        self.geneid = geneid
    def applyTSNE(self):
        self.tsne = TSNE(n_components=3, verbose=0, perplexity=40, n_iter=500).fit_transform(self.pcs)
    # . . . . . . .

# use the pcs from the previous function to perform a tsne analysis
def perform_TSNE(pca_data):
    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=500)
    tsne_results = tsne.fit_transform(pca_data)
    return tsne_results

if __name__=='__main__':
    print('')
    #parser = argparse.ArgumentParser(description='Generates a dosage file with the wanted SNPs')
    #parser.add_argument('-s', '--sample', action='store', nargs='+', help='A set of snps to look for')
    #parser.add_argument('-b', '--beta_values', action='store', help='methylation file ')
    #parser.add_argument('-a', '--annotation', action='store', help='file with CpGs per gene')
    #parser.add_argument('-o', '--output', action='store', help='output pc1 and pc2 for the TSNE')
    # .
    #args = parser.parse_args(
    #beta_values = pd.read_csv(args.beta_values)
    #sample_annotation = pd.read_csv(args.sample)

    #args = parser.parse_args()
    beta_values = pd.read_csv('~/colon-cancer/processed_data/beta_values.csv', index_col=0)
    sample_annotation = pd.read_csv('~/colon-cancer/MethylationEPIC_Sample_Sheet_2018.csv')

    cpg_gene_mapping_file = '/home/people/diomar/colon-cancer/cpg_annotation.txt'
    pca = decomposition.PCA()
    thres = 0.9
    tsne_gene_list = []
    print('computing PCAs per gene')
    with open(cpg_gene_mapping_file) as fin:
        for lines in fin:
            tokz = lines.split()
            gene = tokz[0]
            snps = tokz[1:]
            beta_gene = beta_values.loc[snps, :].transpose()
            pca_transformed = pca.fit_transform(beta_gene)
            max_pcs = 0
            exp_var = 0.0
            for i in range(0, pca.explained_variance_ratio_.size):
                exp_var+=pca.explained_variance_ratio_[i]
                if (exp_var >= thres):
                    max_pcs = i
                    break
            if (max_pcs > 0):
                tsne_gene_list.append(tsne_gene(pca_transformed[:,0:max_pcs], gene))
    print('compute TSNE per gene')
    # compute tsne coordiantes
    for genes in tsne_gene_list:
        genes.applyTSNE()
        df = pd.DataFrame(genes.tsne).to_csv('/home/people/diomar/colon-cancer/genes/'+genes.geneid+'.csv', index=False, header=False)

    #print('Total SNPS:{S} --- Max PCS: {P}'.format(S=beta_gene.shape[1], P=max_pcs))
    # iterate through cpgs per gene file
    #with open(args.annotation) as f:
    #    pass
    #data =  pd.read_csv(args.input, index_col=0)
    #data = data.transpose()
    #labels = pd.read_csv(args.labels)
    #pca, pca_t = perform_pcs(data,2)
    #tsne_results = perform_TSNE(pca_t)
    #.
    #tsne_results_df = pd.DataFrame(tsne_results)
    #tsne_results_df.columns = ['pc1', 'pc2']
    #tsne_results_df['label'] = labels.label
    #tsne_results_df.to_csv(args.output, index=False)
