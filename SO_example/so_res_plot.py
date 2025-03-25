import matplotlib.pyplot as plt
from csv_so_processor import *
import numpy as np
import copy 
import json
from typing import Union

d_sc = 5
n_rounds = 2*d_sc
result_file = 'sc' + str(d_sc) + '_' + str(n_rounds) + '_SE_rds_SO.csv'

so_res = SoftOutput(result_file)
so_res.read_through_custom_counts()

num_survived_shots = so_res.shots - so_res.discards




# freq - gap plot
plt.figure()
plt.plot(so_res.gap_vals, np.array(so_res.at_gap_shots)/num_survived_shots)
plt.yscale('logit')
plt.xlabel('Soft Output (dB)')
plt.ylabel('Frequency')
plt.savefig('Freq_at_Soft_output.pdf')

# logic err rate - gap plot
plt.figure()
plt.plot(so_res.gap_vals,np.array(so_res.at_gap_hits)/(np.array(so_res.at_gap_shots)))
plt.yscale('logit')
plt.xlabel('Soft Output (dB)')
plt.ylabel('Logical Error Rate (at this gap; per shot)')
plt.savefig('LER_at_Soft_output.pdf')

# logic err rate (after post select) - gap plot

plt.figure()
plt.plot(so_res.gap_vals,np.array(so_res.geq_gap_hits)/(np.array(so_res.geq_gap_shots)))
plt.yscale('logit')
plt.grid(which='both',axis='both')
plt.xlabel('Soft Output (dB)')
plt.ylabel('Logical Error Rate (above this gap; per shot)')
plt.savefig('LER_above_Soft_output.pdf')

# survival rate (after post select) - gap plot

plt.figure()
plt.plot(so_res.gap_vals,(np.array(so_res.geq_gap_shots))/num_survived_shots)
# plt.yscale('logit')
plt.grid(which='both',axis='both')
plt.xlabel('Soft Output (dB)')
plt.ylabel('Survival Rate (above this gap; per shot)')
plt.savefig('Survival_rate_above_Soft_output.pdf')

plt.figure()
plt.plot(1-(np.array(so_res.geq_gap_shots)/num_survived_shots),
         np.array(so_res.geq_gap_hits)/(np.array(so_res.geq_gap_shots)))
plt.grid(which='both',axis='both')
plt.yscale('logit')
plt.xscale('logit')
x_right_index = 0
for i in range(len(so_res.geq_gap_hits)):
    if so_res.geq_gap_hits[i] == 0:
        x_right_index = i
        break
if x_right_index == 0:
    x_right_index = len(so_res.geq_gap_hits)-1
 
x_right = 1-(so_res.geq_gap_shots[x_right_index]/num_survived_shots)
plt.xlim(right=x_right*1.2)
plt.xlim(left=0)
plt.xlabel('Survival Rate')
plt.ylabel('Logical Error Rate (per shot)')
plt.savefig('LER_vs_Survival_Rate.pdf')

