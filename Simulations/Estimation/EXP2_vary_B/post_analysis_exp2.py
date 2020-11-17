import numpy as np
import glob, os
import csv

# get mean, se of MLE for the setting N, M, B, rate_obs
def postAnalyze(N, M, rate_obs, B):
    mapp_B = {50000: 'BX1', 250000: 'BX5', 500000: 'BX10', 1000000: 'BX20'}
    fn_pattern1 = 'N'+str(N)+'M'+str(M)+'_'+mapp_B[B]+'.o*'
    # fn_pattern2 = 'N'+str(N)+'M'+str(M)+'_'+str(rate_obs)+'_PR_'+mapp_B[B]+'.o*'

    # get rate_detection for ER sequences
    mean, se, estimates, seeds  = extractSS(fn_pattern1, 'ER')
    res_all = [[N] + [M] +[round(float(M)/float(N)/rate_obs, ndigits=2)]+[B]+['ER']+i+[j] for i, j in zip(estimates, seeds)]
    res_summary = [[N, M, round(float(M)/float(N)/rate_obs, ndigits=2), B, 'ER']+ list(mean) + list(se) +[len(seeds)]]
    return res_all, res_summary

# extract sufficient statistics (seed, estimates, iterations) from scripts with pattern (i.e. a setting)
def extractSS(fn_pattern, process):
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

N = 20
rate_obs = 0.6
outfile = 'data_exp2.csv'
outfile2 = 'data_exp2_all.csv'
print('N=%d, rate_obs=%.1f'%(N, rate_obs))
data = [['N', 'M', 'scaled_t_obs', 'B', 'process', 'mean_p', 'mean_q', 'mean_gamma', 'mean_alpha', 'mean_beta', 'mean_iter', 'se_p', 'se_q', 'se_gamma', 'se_alpha', 'se_beta', 'se_iter', 'total']]
data_all = [['N', 'M', 'scaled_t_obs', 'B', 'process', 'p_est', 'q_est', 'gamma_est', 'alpha_est', 'beta_est', 'iter', 'seed']]

for B in [50000*1, 50000*5, 50000*10, 50000*20]:
    print ('->B=%d'%B)
    for M in [50]:
        print('-->M=%d'%M)
        res_all, res_summary = postAnalyze(N, M, rate_obs, B)
        data = data + res_summary
        data_all = data_all + res_all
    print('\n')

with open(outfile, "w+") as f:
    wr = csv.writer(f)
    wr.writerows(data)

with open(outfile2, "w+") as f:
    wr = csv.writer(f)
    wr.writerows(data_all)
