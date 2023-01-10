import numpy as np
import tensorflow as tf
import sys


input_file = snakemake.input[0] #sys.argv[1] # "good.fa"
input_file2 = snakemake.input[1] #sys.argv[2]# "bad.fa"
output_file = snakemake.output[0] #"read_embedding.txt"
k = snakemake.params[0]
encoding = {"A":1,"C":2,"G":3,"T":4,"N":0}
max_sz = 1*10**6


kmers = np.zeros((max_sz,k,5),dtype='float32')
labels = np.zeros((max_sz),dtype='float32')

r = 0
f1 = open(input_file)
f2 = open(input_file2)
while (line1 := f1.readline()) and (line2 := f2.readline()):
    if line1[0] == ">":
        continue
    kmer = line1.split()[0]
    if len(kmer) < 50:
        continue
    for i,n in enumerate(kmer):
        kmers[r][i][encoding[n]] = 1
    labels[r] = 0
    r += 1
    if r >= max_sz:
        break
    kmer = line2.split()[0]
    for i,n in enumerate(kmer):
        kmers[r][i][encoding[n]] = 1
    labels[r] = 1
    r += 1
    if r >= max_sz:
        break
kmers = kmers[:r]
labels = labels[:r]

print(kmers[0])
data = tf.data.Dataset.from_tensor_slices((kmers, labels))
data.save(output_file)
