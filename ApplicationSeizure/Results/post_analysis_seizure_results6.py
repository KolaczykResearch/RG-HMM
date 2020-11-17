import numpy as np
import glob, os
import csv

fn_pattern = 'P1S1_BX10_1_0.5_178_209_left.o*'
# fn_pattern = 'P1S2_BX10_1_0.5_237_255_left.o*'
# fn_pattern = 'P1S3_BX10_1_0.5_98_164_left.o*'
# fn_pattern = 'P1S3_BX10_1_0.5_183_193_left.o*'

seizure_name = 'S3'

# get bf from single run in one script
def get_est_one(filename):
    info = []
    with open(filename) as file:
        for line in file.readlines():
            if 'para_est_ER=' in line:
                estimates_ER = [float(i) for i in line[13:].strip().split(' ')]
                # print(estimates_ER, '\n')
                info = info + estimates_ER
            if 'para_est_PR=' in line:
                estimates_PR = [float(i) for i in line[13:].strip().split(' ')]
                # print(estimates_PR, '\n')
                info = info + estimates_PR
            if 'Bayes Factor is: ' in line:
                BF = [float(line[17:].strip())]
                info = info + BF
            if '...-' in line:
                string = line[line.index('-'):]
                loglik = [float(string.strip())]
                info = info + loglik
    if len(info)<13: print('Check file:', filename, '\n')
    return info
    
# get estimates from all scripts and take average
def avg_estimates(fn_pattern):
    mat_est = []
    for fname in glob.glob(fn_pattern):
        # print(fname)
        est = get_est_one(fname)
        if len(est)==13:
            mat_est.append(est)
    n_rep = len(mat_est)
    mat_est = np.array(mat_est)
    avg = np.mean(mat_est, axis=0)
    sd = np.sqrt(np.var(mat_est, axis=0))  
    print('Mean of est_para_est_ER: ', avg[0:5], '\n')
    print('                     sd: ', sd[0:5], '\n')
    print('Mean of est_para_est_PR: ', avg[5:10], '\n')
    print('                     sd: ', sd[5:10], '\n')
    print('Mean of loglik_ER: ', avg[10], '\n')
    print('               sd: ', sd[10], '\n')
    print('Mean of loglik_PR: ', avg[11], '\n')
    print('               sd: ', sd[11], '\n')
    print('Mean of BF: ', avg[12], '\n')
    print('        sd: ', sd[12], '\n')
#    print('Mean of BF: ', avg[10], '\n')
#    print('        sd: ', sd[10], '\n')
    print('replicates: ', n_rep)

print('Avg for ', seizure_name, '\n')
print(fn_pattern, '\n')
avg_estimates(fn_pattern)
