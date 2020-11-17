import numpy as np
import glob, os
import csv

# extract sufficient statistics (seed, estimates, iterations) from scripts with pattern (i.e. a setting)
def extractSS(fn_pattern):
    results = []
    for fname in glob.glob(fn_pattern):
        info = get_info(fname)
        results.append(info)
    # total available results
    estimates = [item[1:] for item in results if len(item)>1]
    seeds = [item[0] for item in results if len(item)>1]
    if len(seeds) > 0:
        mean = np.mean(estimates, axis=0)
        se = np.sqrt(np.var(estimates, axis=0))
        print('    Mean: ', mean)
        print('    SE:', se)
    else:
        mean = [0]*6
        se = [0]*6
    print('    Average over %d results, '%len(seeds))
    return mean, se, estimates, seeds

# get seed, estimates, iterations from a single script
def get_info(filename):
    last_num = filename[-4::]
    last_num = int(last_num[last_num.index('.')+1:])
    seed = 13280 + last_num
    info = [seed]
    with open(filename) as file:
        for line in file.readlines():
            if 'Estimates (p, q, gamma, alpha, beta):' in line:
                estimates = [float(i) for i in line[38:].strip().split(' ')]
                # print(estimates, '\n')
                info = info + estimates
            if 'iteration:' in line:
                iteration = float(line[11:].strip())
                # print(iteration, '\n')
                info = info + [iteration]
    if len(info)<=1: print('Check file:', filename, '\n')
    return info

fn_patterns = ['N20M50_0.6_ER_BX1.o*', 'N20M50_0.6_PR_BX1.o*']
for fn_pattern in fn_patterns:
    print(fn_pattern, '\n')
    extractSS(fn_pattern)
