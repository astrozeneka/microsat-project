
from scipy import stats
import numpy as np

if __name__ == '__main__':

    rvs = stats.norm.rvs(loc=5, scale=10, size=500)
    print(stats.ttest_1samp(rvs, 5.0))
    print()