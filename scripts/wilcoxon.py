import pandas as pd
import scipy as sp

from scipy import stats

DATASETS = '../depmap_data/'

def load_dataset(filename):
    return pd.read_csv(filename)

if __name__ == '__main__':
    # /depmap_data
    # ...CRISPR_gene_effect.csv
    # ...CCLE_expression.csv
    # /scripts
    # ... wilcoxon.py
    gene_effect = load_dataset(DATASETS + 'CRISPR_gene_effect.csv')
    expression = load_dataset(DATASETS + 'CCLE_expression.csv')


    # slicing two datasets to have same cell lines
    cell_lines_1 = gene_effect['DepMap_ID']
    cell_lines_2 = expression.iloc[:, 0]
    cell_lines = pd.Series(list(set(cell_lines_1) & set(cell_lines_2)))
    gene_effect = gene_effect[gene_effect['DepMap_ID'].isin(cell_lines)]
    print(gene_effect)
    expression = expression[expression.iloc[:, 0].isin(cell_lines)]
    print(expression)
    # they have same cell lines (1005 rows)

    total_cell_lines = gene_effect.shape[0]
    print(total_cell_lines)
    benchmark = total_cell_lines // 3

    wilcoxon_results = {}

    for i in range(1, 5):
        x = gene_effect.iloc[:, [0,i]]
        gene_A = list(x)[1]

        wilcoxon_results[gene_A] = {}

        for j in range(1, 5):
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
            print(gene_A, gene_B)
            print(wilcoxon)

            # saving in matrix of wilcoxon results
            wilcoxon_results[gene_A][gene_B] = wilcoxon

    print(wilcoxon_results)

    