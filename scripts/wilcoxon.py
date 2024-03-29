from lib2to3.pytree import WildcardPattern
from operator import ge
import pandas as pd
import scipy as sp
import json
import os
from multiprocessing import Process, Manager, Queue
from scipy import stats

# /depmap_data
# ...CRISPR_gene_effect.csv
# ...CCLE_expression.csv
# /scripts
# ... wilcoxon.py


def load_dataset(filename):
    return pd.read_csv(filename)

def compute_essentiality_association(datapath, outpath=".", method="ISLE", n=10, tissue='None', bimodal=False):
    """
    datapath: path to depmap data (folder that includes expression and essentiality data)
    #   essentiality_path: path to depmap essentiality data
    #   gene_exp_data_path: Path to depmap gene expression data
    method: method to compute the essentiality association
    - "ISLE" (default): Wilcoxon tercile method
    n: number of genes
    - 10 (default): 10 of gene A x 10 of gene B
    - "all": runs the whole dataset
    tissue: filter by specific tissue type
    - "None" (default): do not filter by tissue type (take all tissues)
    - "pan": filter by Pancreatic Cancer tissues
    = "lung": filter by Lung Cancer tissues
    - "breast": filter by Bresat Cancer tissues
    bimodal: True to only consider bimodal genes in the gene expression data. False to consider all genes.
    """

    essentiality_path = os.path.join(datapath, 'CRISPR_gene_effect.csv')
    expression_path = os.path.join(datapath, 'CCLE_expression.csv')

    gene_effect = load_dataset(essentiality_path)
    expression = load_dataset(expression_path)

    # slicing two datasets to have same cell lines
    cell_lines_1 = gene_effect['DepMap_ID']
    cell_lines_2 = expression.iloc[:, 0]
    cell_lines = pd.Series(list(set(cell_lines_1) & set(cell_lines_2)))
    bimodals = pd.read_csv('../data/bimodal.csv')

    genes = gene_effect.columns
    bimodal_filtered=['DepMap_ID']
    bimodals = bimodals['Symbol']
    for gene in genes:
        spliced = gene.split(' (')[0]
        if spliced in list(bimodals):
            bimodal_filtered.append(gene)

    if bimodal:
        gene_effect = gene_effect[bimodal_filtered]

    gene_effect = gene_effect[gene_effect['DepMap_ID'].isin(cell_lines)]
    expression = expression[expression.iloc[:, 0].isin(cell_lines)]
    # they have same cell lines (1005 rows)

    if (tissue != "None"):
        genes_types = tissue_specific_info(datapath)
        if (tissue == 'pan'):
            print('pan')
            filtered_genes = genes_types[genes_types['primary_disease'] == 'Pancreatic Cancer']
        elif (tissue == 'lung'):
            filtered_genes = genes_types[genes_types['primary_disease'] == 'Lung Cancer']
            print('lung')
        elif (tissue == 'breast'):
            print('breast')
            filtered_genes = genes_types[genes_types['primary_disease'] == 'Breast Cancer']
        else:
            print('invalid tissue type! Available: "None"/"pan"/"lung"/"breast"')
            return

        filtered = pd.Series(list(filtered_genes.iloc[:,0]))

        gene_effect = gene_effect[gene_effect['DepMap_ID'].isin(filtered)]
        expression = expression[expression['Unnamed: 0'].isin(filtered)]

        if (method == "ISLE"):
            wilcoxon_results, df = perform_isle_method(gene_effect, expression, n)
            json_name = 'wilcoxon_results_' + tissue + '_n=' + str(n) + '.json'
            with open(os.path.join(outpath, json_name), 'w') as fp:
                json.dump(wilcoxon_results, fp)
            return df

    else:  # all tissues
        if (method == "ISLE"):
            wilcoxon_results, df = perform_isle_method(gene_effect, expression, n)
            json_name = 'wilcoxon_results_pancan_n=' + str(n) + '.json'
            with open(os.path.join(outpath, json_name), 'w') as fp:
                json.dump(wilcoxon_results, fp)
            return df

    print('done!')

    return


def perform_isle_method(gene_effect, expression, n, useParallel = False):
    # Assuming gene_effect and expression already has same cell lines
    total_cell_lines = gene_effect.shape[0]
    print(total_cell_lines)
    benchmark = total_cell_lines // 3


    if n == "all":
        total_gene_As = len(gene_effect.columns)
        total_gene_Bs = len(expression.columns)
    else: 
        total_gene_As = n+1
        total_gene_Bs = n+1

    if useParallel == 0:
      wilcoxon_results = {}
      for i in range(1, total_gene_As):

        if (i % 10 == 0):
            print(i)

        x = gene_effect.iloc[:, [0,i]]
        gene_A = list(x)[1]
        wilcoxon_results[gene_A] = {}
          
        for j in range(1, total_gene_Bs):
              y = expression.iloc[:, [0,j]]
              gene_B = list(y)[1]
  
              # get cell lines in bottom 1/3 and top 1/3 of expression for a gene B
              y = y.sort_values(by=list(y)[1:])
              bottom_third = y.head(benchmark).iloc[:, 0]
              top_third = y.tail(benchmark).iloc[:, 0]
              bottom_third = list(bottom_third)
              top_third = list(top_third)

              # get essentiality_bottom = essentiality_A[cell line is in bottom 1/3 of gene B]
              # and essentiality_top    = essentiality_A[cell line is in top 1/3 of gene B]
              essentiality_top = x[x['DepMap_ID'].isin(top_third)]
              essentiality_bottom = x[x['DepMap_ID'].isin(bottom_third)]
              essentiality_top = list(essentiality_top.iloc[:, 1])
              essentiality_bottom = list(essentiality_bottom.iloc[:,1])
 
              # perform the wilcoxon test
              wilcoxon = sp.stats.wilcoxon(essentiality_bottom, essentiality_top)
              # print(gene_A, gene_B)
              # print(wilcoxon)

              # saving in matrix of wilcoxon results
              wilcoxon_results[gene_A][gene_B] = wilcoxon
    
      wilcoxon_results_df = pd.DataFrame.from_dict(wilcoxon_results)
      return wilcoxon_results, wilcoxon_results_df
    
    else:   # Use parallel
        manager = Manager()
        wilcoxonDict = manager.dict()
     
        queueList = [Queue() for i in range(numProcs)]
        procList = [Process(target=isle_submethod, args=(queueList[i], FullModelResDict)) for i in range(numProcs)]

        for proc in procList:
          proc.start()

        for i in range(len(parGrid)):
          queueList[i%numProcs].put(parGrid[i])
 
        for queue in queueList:
          queue.put(-1)

        for proc in procList:
          proc.join()

        wilcoxon_results = pd.DataFrame.from_dict(wilcoxonDict._getvalue())
        return wilcoxon_results, wilcoxonDict._getvalue()

          

def isle_submethod(gene_effect, expression, i, j, resDict):
    y = expression.iloc[:, [0,j]]
    gene_B = list(y)[1]
    
    # get cell lines in bottom 1/3 and top 1/3 of expression for a gene B
    y = y.sort_values(by=list(y)[1:])
    bottom_third = y.head(benchmark).iloc[:, 0]
    top_third = y.tail(benchmark).iloc[:, 0]
    bottom_third = list(bottom_third)
    top_third = list(top_third)

    # get essentiality_bottom = essentiality_A[cell line is in bottom 1/3 of gene B]
    # and essentiality_top    = essentiality_A[cell line is in top 1/3 of gene B]
    essentiality_top = x[x['DepMap_ID'].isin(top_third)]
    essentiality_bottom = x[x['DepMap_ID'].isin(bottom_third)]
    essentiality_top = list(essentiality_top.iloc[:, 1])
    essentiality_bottom = list(essentiality_bottom.iloc[:,1])

    # perform the wilcoxon test
    wilcoxon = sp.stats.wilcoxon(essentiality_bottom, essentiality_top)

    # saving in matrix of wilcoxon results
    resDict[gene_A][gene_B] = wilcoxon


def tissue_specific_info(datapath):
    sample_info_path = os.path.join(datapath, 'sample_info.csv')
    info = load_dataset(sample_info_path)
    # umm = info['primary_disease'].unique()

    # for x in umm:
    #     genes_types[genes_types['primary_disease'] == 'Pancreatic Cancer']

    # umm = info['primary_disease']
    # print(umm.unique())
    #     ['Kidney Cancer' 'Leukemia' 'Lung Cancer' 'Non-Cancerous' 'Sarcoma'
    #  'Lymphoma' 'Colon/Colorectal Cancer' 'Pancreatic Cancer' 'Gastric Cancer'
    #  'Rhabdoid' 'Endometrial/Uterine Cancer' 'Esophageal Cancer'
    #  'Breast Cancer' 'Brain Cancer' 'Ovarian Cancer' 'Bone Cancer' 'Myeloma'
    #  'Head and Neck Cancer' 'Bladder Cancer' 'Skin Cancer' 'Bile Duct Cancer'
    #  'Prostate Cancer' 'Cervical Cancer' 'Thyroid Cancer' 'Neuroblastoma'
    #  'Eye Cancer' 'Liposarcoma' 'Gallbladder Cancer' 'Teratoma' 'Unknown'
    #  'Liver Cancer' 'Adrenal Cancer' 'Embryonal Cancer']
    genes_types = info[['DepMap_ID', 'primary_disease']]

    # for x in umm:
    #     print(x)
    #     print(genes_types[genes_types['primary_disease'] == x].shape)
    return genes_types


if __name__ == '__main__':

    # tissue_specific_info()


    datapath = '../depmap_data/'

    essentiality_path = datapath + 'CRISPR_gene_effect.csv'
    expression_path = datapath + 'CCLE_expression.csv'
    print(essentiality_path)
    compute_essentiality_association(datapath, "ISLE", 10, 'pan', True)
 
    # compute_essentiality_association(essentiality_path, expression_path, "ISLE", 10, 'breast')
    # compute_essentiality_association(essentiality_path, expression_path, "ISLE", 10, 'lung')

    
