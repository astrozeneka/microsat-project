
from itertools import combinations
from scipy.stats import ttest_ind

if __name__ == '__main__':
    data = {
        'a': [1, 1, 1, 1, 1],
        'b': [20, 20, 20],
        'c': [1, 2],
    }
    for list1, list2 in combinations(data.keys(), 2):
        t, p = ttest_ind(data[list1], data[list2])
        print(list1, list2, p)