from scipy import stats

import pandas as pd
import matplotlib.pyplot as plt
import json
import scipy as sp
import numpy as np

# https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
# Going to use Benjamin Hochberg method
# https://www.researchgate.net/post/Can_anyone_explain_how_to_calculate_adjusted_p-values_q-values_following_Benjamini_Hochberg_correction
import statsmodels.stats.multitest as smt

def random_wilcoxons(n=1000):
    # 1000 random draws of wilcoxon
    sample_teststats = []
    sample_pvalues = []
    for i in range(1000):
        result = sp.stats.wilcoxon(np.random.rand(335), np.random.rand(335))
        sample_teststats.append(result.statistic)
        sample_pvalues.append(result.pvalue)

    plt.hist(sample_teststats, density=True, bins=32)
    plt.xlim(xmin=15000, xmax=35000)
    plt.ylabel('Frequency')
    plt.xlabel('Wilcoxon Test Statistics')
    plt.title('1000 random draws of wilcoxon with same list lengths (335)')
    plt.show()

    plt.hist(sample_pvalues, density=True, bins=32)
    plt.ylabel('Frequency')
    plt.xlabel('Wilcoxon p-values')
    plt.title('1000 random draws of wilcoxon with same list lengths (335)')
    plt.show()


if __name__ == '__main__':

    # random_wilcoxons()

    # Opening JSON file
    with open('wilcoxon_results.json') as json_file:
        data = json.load(json_file)
        i = 0
        test_stats = []
        p_values = []
        for geneA, scores in data.items():
            for geneB, stats in scores.items():
                test_stats.append(stats[0])
                p_values.append(stats[1])
        
        print('total stats:', len(test_stats))
        print(max(test_stats), max(p_values))

        # plt.hist(test_stats, density=True, bins=1005)
        # plt.xlim(xmin=15000, xmax=35000)
        # plt.ylabel('Frequency')
        # plt.xlabel('Wilcoxon Test Statistics')
        # plt.title('1005 genes x 1005 genes Wilcoxon results')
        # plt.show()

        # plt.hist(p_values, density=True, bins=1005)
        # plt.ylabel('Frequency')
        # plt.ylim(ymin=0, ymax=5)
        # plt.xlabel('Wilcoxon Test P-Values')
        # plt.title('1005 genes x 1005 genes Wilcoxon results')
        # plt.show()
        

        # https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
        # Going to use Benjamin Hochberg method
        # https://www.researchgate.net/post/Can_anyone_explain_how_to_calculate_adjusted_p-values_q-values_following_Benjamini_Hochberg_correction
        
        # p_values must not contain nan or result will all be nan
        print(np.isnan(p_values).any())
        p_values = np.nan_to_num(p_values)
        print(np.isnan(p_values).any())

        q_values = smt.multipletests(p_values, method='fdr_bh') #Benjamini/Hochberg method

        print(q_values)
        q_values_list = q_values[1].tolist()
        false_positives = 0
        for x in q_values[0]:
            if x == True:
                false_positives += 1
        print('false positives:', false_positives)

        valid_q_values = []
        print(len(q_values[0]) == len(q_values[1]))
        for i in range(len(q_values[0])):
            if (q_values[0][i]) == False:
                valid_q_values.append(q_values_list[i])

        plt.hist(valid_q_values, density=True, bins=1005)
        plt.ylabel('Frequency')
        plt.ylim(ymin=0, ymax=5)
        plt.xlabel('Wilcoxon Test Q-Values')
        plt.title('1005 genes x 1005 genes Wilcoxon results')
        plt.show()