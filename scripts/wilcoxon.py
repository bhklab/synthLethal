from lib2to3.pytree import WildcardPattern
import pandas as pd
import scipy as sp
import json

from scipy import stats

# /depmap_data
# ...CRISPR_gene_effect.csv
# ...CCLE_expression.csv
# /scripts
# ... wilcoxon.py
DATASETS = '../depmap_data/'

def load_dataset(filename):
    return pd.read_csv(filename)

def compute_essentiality_association(essentiality_data_path, gene_exp_data_path, method="ISLE", n=10):
    """
    essentiality_data_path: path to depmap essentiality data
    gene_exp_data_path: Path to depmap gene expression data
    method: method to compute the essentiality association
    - "ISLE" (default): Wilcoxon tercile method
    n: number of genes
    - 10 (default): 10 of gene A x 10 of gene B
    - "all": runs the whole dataset
    """
    gene_effect = load_dataset(essentiality_data_path)
    expression = load_dataset(gene_exp_data_path)

    # slicing two datasets to have same cell lines
    cell_lines_1 = gene_effect['DepMap_ID']
    cell_lines_2 = expression.iloc[:, 0]
    cell_lines = pd.Series(list(set(cell_lines_1) & set(cell_lines_2)))
    gene_effect = gene_effect[gene_effect['DepMap_ID'].isin(cell_lines)]
    expression = expression[expression.iloc[:, 0].isin(cell_lines)]
    # they have same cell lines (1005 rows)

    if (method == "ISLE"):
        perform_isle_method(gene_effect, expression, n)

    print('done!')

    return

def perform_isle_method(gene_effect, expression, n):
    # Assuming gene_effect and expression already has same cell lines
    total_cell_lines = gene_effect.shape[0]
    print(total_cell_lines)
    benchmark = total_cell_lines // 3

    wilcoxon_results = {}

    total_gene_As = n+1
    total_gene_Bs = n+1
    if n == "all":
        total_gene_As = len(gene_effect.columns)
        total_gene_Bs = len(expression.columns)

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

    with open('wilcoxon_results2.json', 'w') as fp:
        json.dump(wilcoxon_results, fp)

    return wilcoxon_results_df


if __name__ == '__main__':
    essentiality_path = DATASETS + 'CRISPR_gene_effect.csv'
    expression_path = DATASETS + 'CCLE_expression.csv'
    compute_essentiality_association(essentiality_path, expression_path, "ISLE", 10)


    