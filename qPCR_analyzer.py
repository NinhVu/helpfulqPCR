#!/bin/python
# qPCR_analyzer.py created by Ninh Vu 2/__/16
# This program graphs your qPCR results and generate output text file. You can erase any point or points to improve standard and triangulate true concentration.
# In terminal type qPCR_analyzer.py input_file.txt results.txt

import sys, matplotlib.pyplot as plt, time, scipy, math, numpy as np
from scipy import stats

results = open(sys.argv[1], 'r')
concentration = open(sys.argv[2], 'w')

# read and discard first nine lines of results file
for i in range(9):
    header = results.readline()

# define lists use in reading file and initial table
s,t1,t2,t3,t4,t5,t6,t7,t8,t9 = [],[],[],[],[],[],[],[],[],[]

# read data and create lists
data = results.readline()
while data:
    rowItems = data.split("\t")
    for i in rowItems:
        if i=='S1' or i=='S2' or i=='S3' or i=='S4' or i=='S5' or i=='S6':
            s.append(rowItems[6])
            s = [float(i) for i in s]
        elif i == "T1_1" or i == "T1_2" or i == "T1_3" or i == "T1_4":
            t1.append(rowItems[6])
            t1 = [float(i) for i in t1]
        elif i == "T2_1" or i == "T2_2" or i == "T2_3" or i == "T2_4":
            t2.append(rowItems[6])
            t2 = [float(i) for i in t2]
        elif i == "T3_1" or i == "T3_2" or i == "T3_3" or i == "T3_4":
            t3.append(rowItems[6])
            t3 = [float(i) for i in t3]
        elif i == "T4_1" or i == "T4_2" or i == "T4_3" or i == "T4_4":
            t4.append(rowItems[6])
            t4 = [float(i) for i in t4]
        elif i == "T5_1" or i == "T5_2" or i == "T5_3" or i == "T5_4":
            t5.append(rowItems[6])
            t5 = [float(i) for i in t5]
        elif i == "T6_1" or i == "T6_2" or i == "T6_3" or i == "T6_4":
            t6.append(rowItems[6])
            t6 = [float(i) for i in t6]
        elif i == "T7_1" or i == "T7_2" or i == "T7_3" or i == "T7_4":
            t7.append(rowItems[6])
            t7 = [float(i) for i in t7]
        elif i == "T8_1" or i == "T8_2" or i == "T8_3" or i == "T8_4":
            t8.append(rowItems[6])
            t8 = [float(i) for i in t8]
        elif i == "T9_1" or i == "T9_2" or i == "T9_3" or i == "T9_4":
            t9.append(rowItems[6])
            t9 = [float(i) for i in t9]

    data = results.readline()

results.close()

# average standards as one list and average each sample as seperate list
S = np.array(s)
Savg = (S[::3] + s[1::3] + S[2::3])/3
T1 = [(a + b) / 2 for a, b in zip(t1[::2], t1[1::2])]
T2 = [(a + b) / 2 for a, b in zip(t2[::2], t2[1::2])]
T3 = [(a + b) / 2 for a, b in zip(t3[::2], t3[1::2])]
T4 = [(a + b) / 2 for a, b in zip(t4[::2], t4[1::2])]
T5 = [(a + b) / 2 for a, b in zip(t5[::2], t5[1::2])]
T6 = [(a + b) / 2 for a, b in zip(t6[::2], t6[1::2])]
T7 = [(a + b) / 2 for a, b in zip(t7[::2], t7[1::2])]
T8 = [(a + b) / 2 for a, b in zip(t8[::2], t8[1::2])]
T9 = [(a + b) / 2 for a, b in zip(t9[::2], t9[1::2])]
T1logC,T2logC,T3logC,T4logC,T5logC,T6logC,T7logC,T8logC,T9logC = [],[],[],[],[],[],[],[],[]

# calculate y = mx + b
stand = [20, 2, 0.2, 0.02, 0.002, 0.0002]
standlog = map(math.log10, stand)
slope, intercept, r_value, p_value, std_err = stats.linregress(standlog, Savg)
print('R^2: ', r_value**2)
print('slope: ', slope)
print('y-intercept: ', intercept)
print('S: ', S)
print('Savg: ', Savg)

T1logC[:] = [(x-intercept)/slope for x in T1]
T2logC[:] = [(x-intercept)/slope for x in T2]
T3logC[:] = [(x-intercept)/slope for x in T3]
T4logC[:] = [(x-intercept)/slope for x in T4]
T5logC[:] = [(x-intercept)/slope for x in T5]
T6logC[:] = [(x-intercept)/slope for x in T6]
T7logC[:] = [(x-intercept)/slope for x in T7]
T8logC[:] = [(x-intercept)/slope for x in T8]
T9logC[:] = [(x-intercept)/slope for x in T9]


plt.plot(standlog,Savg, 'bs')
plt.plot(T1logC,T1, 'ro')
plt.plot(T2logC,T2, 'ro')
plt.plot(T3logC,T3, 'ro')
plt.plot(T4logC,T4, 'ro')
plt.plot(T5logC,T5, 'ro')
plt.plot(T6logC,T6, 'ro')
plt.plot(T7logC,T7, 'ro')
plt.plot(T8logC,T8, 'ro')
plt.plot(T9logC,T9, 'ro')

fit = np.polyfit(standlog,Savg,1)
fit_fn = np.poly1d(fit)
plt.plot(standlog,Savg, 'yo', standlog, fit_fn(standlog), '-k')



plt.show()



