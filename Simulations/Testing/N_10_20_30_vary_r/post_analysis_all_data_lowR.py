import numpy as np
import glob, os
import csv

# get rate of detection (se), and table of 0/1 results for the setting N, M, B, rate_obs
def postAnalyze(N, M, rate_obs, B):
    arr = []
    test_res = []
    mapp_B = {50000: 'BX1', 250000: 'BX5', 500000: 'BX10'}
    # fn_pattern1 = folder +'/N'+str(N)+'M'+str(M)+'_'+str(rate_obs)+'_ER_'+mapp_B[B]+'.o*'
    # fn_pattern2 = folder + '/N'+str(N)+'M'+str(M)+'_'+str(rate_obs)+'_PR_'+mapp_B[B]+'.o*'
    fn_pattern1 = 'N'+str(N)+'M'+str(M)+'_'+str(rate_obs)+'_ER_'+mapp_B[B]+'.o*'
    fn_pattern2 = 'N'+str(N)+'M'+str(M)+'_'+str(rate_obs)+'_PR_'+mapp_B[B]+'.o*'
    
    # get rate_detection for ER sequences
    rate_ER, n_total, binary, seeds, diffs = extractSS(fn_pattern1, 'ER')
    binary_ER = [[i, N, round(float(M)/float(N)/rate_obs, ndigits=2), B, 'ER', j, k, rate_obs] for i, j, k in zip(binary, seeds, diffs)]
    test_res = test_res + binary_ER
    se = np.sqrt(rate_ER*(1-rate_ER)/n_total)
    print('    %.4f, %.4f'%(rate_ER, se))
    arr1 = [N, round(float(M)/float(N)/rate_obs, ndigits=2), B, rate_ER, 'ER', n_total, se, rate_obs]
    arr.append(arr1)

    # get rate_detection for PR sequences
    rate_PR, n_total, binary, seeds, diffs = extractSS(fn_pattern2, 'PR')
    binary_PR = [[i, N, round(float(M)/float(N)/rate_obs, ndigits=2), B, 'PR', j, k, rate_obs] for i, j, k in zip(binary, seeds, diffs)]
    test_res = test_res + binary_PR
    se = np.sqrt(rate_PR*(1-rate_PR)/n_total)
    print('    %.4f, %.4f'%(rate_PR, se))
    arr2 = [N, round(float(M)/float(N)/rate_obs, ndigits=2), B, rate_PR, 'PR', n_total, se, rate_obs]
    arr.append(arr2)
    return arr, test_res

# extract sufficient statistics (gcc, diff, bf) from scripts with pattern
def extractSS(fn_pattern, process):
    results = []
    for fname in glob.glob(fn_pattern):
        # print(fname)
        gcc = get_gcc_diff_bf(fname)
        results.append(gcc)
        #print(gcc)
    # total available results
    diffs = [item[1] for item in results if len(item)>1]
    diffs = np.array(diffs)
    seeds = [item[0] for item in results if len(item)>1]
    if process=='ER':
        binary = (diffs>0)*1
        rate = np.mean(diffs > 0)
        print('    Average over %d results, ER'%len(diffs))
        #print(diffs)
    else:
        binary = (diffs<0)*1
        rate = np.mean(diffs < 0)
        print('    Average over %d results, PR'%len(diffs))
        #print(diffs)
    return rate, len(diffs), binary, seeds, diffs
    
# get gcc, diff, bf from a single script
def get_gcc_diff_bf(filename):
    #gcc = []
    last_num = filename[-4::]
    last_num = int(last_num[last_num.index('.')+1:])
    seed = 13280 + last_num
    gcc = [seed]
    with open(filename) as file:
        cnt = 1
        for line in file.readlines():
            #print('Line {}: {}'.format(cnt, line.strip()))
            #if 'data gcc' in line:
            #    #print('\n gcc in this line {}'.format(cnt))
            #    gcc = [int(i) for i in line[10:].strip().split(' ')]
            #print(gcc)
            if 'The loglik difference is' in line:
                #print('\n loglik difference in this line {}'.format(cnt))
                diff = float(line[26:].strip())
                #print(diff)
                gcc.append(diff)
            if 'Bayes Factor is' in line:
                #print('\n Bayes factor in this line {}'.format(cnt))
                bf = float(line[17:].strip())
                gcc.append(bf)
            #print(bf)
            cnt = cnt+1
    if len(gcc)<=1: print('Check file:', filename, '\n')
    return gcc
    
# progress = [0.6, 1, 1.4, 1.8] # how much the percolation curve is observed
# progress = [0.40, 0.67, 0.93, 1.2]
scaled_t_obs = 0.67 # (1/1.5)
B = 500000
rate_obs_set = [0.5, 1, 1.5]
outfile = 'data10_20_30_low_r.csv'
outfile2 = 'data10_20_30_low_r_all.csv'

print('B=%d, scaled_t_obs=%.2f'%(B, scaled_t_obs))
data = [['N', 'scaled_t_obs', 'B', 'Rate_of_Detection', 'process', 'total', 'se', 'rate_obs']]
data_all = [['success', 'N', 'scaled_t_obs', 'B', 'process', 'seed', 'loglik_diff', 'rate_obs']]
for N in [10, 20, 30]:
    print ('->N=%d'%N)
    M_set = [int(round(N*1/1.5*rate_obs)) for rate_obs in rate_obs_set]
    for M, rate_obs in zip(M_set, rate_obs_set):
        print('-->M=%d, rate_obs=%.1f'%(M, rate_obs))
        arr, test_res = postAnalyze(N, M, rate_obs, B)
        data = data + arr
        data_all = data_all + test_res
    print('\n')
#ss=postAnalyze(10, 18, 1.5, 50000)

with open(outfile, "w+") as f:
    wr = csv.writer(f)
    wr.writerows(data)

with open(outfile2, "w+") as f:
    wr = csv.writer(f)
    wr.writerows(data_all)

