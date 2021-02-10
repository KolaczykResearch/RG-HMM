import numpy as np
import glob, os
import csv

fn_patterns = ['P1S1_BX10_1_0.5_126_134_left.o*', 'P1S1_BX10_1_0.5_453_463_left.o*', 'P1S1_BX10_1_0.5_609_624_left.o*']
# fn_patterns = ['P1S2_BX10_1_0.5_236_253_left.o*', 'P1S2_BX10_1_0.5_295_318_left.o*', 'P1S2_BX10_1_0.5_627_645_left.o*']
# fn_patterns = ['P1S3_BX10_1_0.5_229_238_left.o*', 'P1S3_BX10_1_0.5_471_487_left.o*', 'P1S3_BX10_1_0.5_256_269_left.o*', 'P1S3_BX10_1_0.5_309_316_left.o*', 'P1S3_BX10_1_0.5_352_375_left.o*', 'P1S3_BX10_1_0.5_385_394_left.o*']

seizure_name = 'S1'
# seizure_name = 'S2'
# seizure_name = 'S3'

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
    print('Mean of est_para_est_ER: ', np.round(avg[0:5], 3), '\n')
    print('                     sd: ', np.round(sd[0:5], 3), '\n')
    print('Mean of est_para_est_PR: ', np.round(avg[5:10], 3), '\n')
    print('                     sd: ', np.round(sd[5:10], 3), '\n')
    print('Mean of loglik_ER: ', np.round(avg[10], 3), '\n')
    print('               sd: ', np.round(sd[10], 3), '\n')
    print('Mean of loglik_PR: ', np.round(avg[11], 3), '\n')
    print('               sd: ', np.round(sd[11], 3), '\n')
    print('Mean of BF: ', np.round(avg[12], 3), '\n')
    print('        sd: ', np.round(sd[12], 3), '\n')   
    print('replicates: ', n_rep)

np.set_printoptions(suppress=True)
print('Avg for ', seizure_name, '\n')
for fn_pattern in fn_patterns:
    print(fn_pattern, '\n')
    avg_estimates(fn_pattern)
