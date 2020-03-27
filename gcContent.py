#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 2020/3/12 21:15
# @Author  : Xu Chi
# @Email   : thisischixu@gmail.com
# @File    : gcContent.py
# @Software: PyCharm
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

######### Parsing arguments ########
desc='''plot GC contenr over sequence\n 
can be run like:\n
python3 gcContent.py -i Ht.fna -o heatmap -w 1000 -s 500 -n 0.2 -L 0.2
'''
parser=argparse.ArgumentParser(description=desc)
parser.add_argument('i',metavar='inFile', type=str)
parser.add_argument('o',metavar='outFile', type=str)
parser.add_argument('w',metavar='windowSize', type=int)
parser.add_argument('s',metavar='step', type=int)
parser.add_argument('n',metavar='numberOfSequence', type=float)
parser.add_argument('L',metavar='LengthOfSequence', type=float)
args=parser.parse_args()
######################################


def calculateGC(seq):
    return (seq.count('G') + seq.count('C'))/len(seq)

def slice(seq, wSize, step):
    seqLength = len(seq)
    GC = []
    if wSize > seqLength:
        print('Window size is longer than sequence length')
        return None
    for i in range(0, seqLength, step):
        end = i + wSize
        if end > seqLength:
            GC.append(round(100*calculateGC(seq[i:seqLength])))
            break
        GC.append(round(100*calculateGC(seq[i:end])))
    return GC

with open(args.i) as f:
    longest = 0
    seqs = {}
    for l in f:
        if l.startswith('>'):
            ID = l.rstrip()[1:10]
            seqs[ID] = seqs.get(ID, [])
        else:
            seqs[ID].append(l.rstrip().upper())
    for k in seqs:
        seqs[k] = slice(''.join(seqs[k]), args.w, args.s)
        longest = max(longest, len(seqs[k]))
data = {}
for k in seqs:
    if len(seqs[k]) > longest*args.n:
        seqs[k] += [100]*(longest-len(seqs[k]))
        data[k] = seqs[k][:int(longest*args.L)]

df = pd.DataFrame(data).T
fig = plt.figure()
sns_plot = sns.heatmap(df, cmap='rainbow', xticklabels=10, yticklabels=8)
plt.title('Window Size = '+ str(args.w) + ' Step = ' + str(args.s))
plt.savefig(args.o)
plt.show()
