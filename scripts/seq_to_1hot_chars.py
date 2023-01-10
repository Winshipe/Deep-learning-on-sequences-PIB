import numpy as np
import tensorflow as tf
import sys
import os

if "snakemake" not in globals():
    class Snakemake:
        input = [sys.argv[1]]
        output = [sys.argv[2]]
        params = [17]
    snakemake = Snakemake() 
input_file = snakemake.input[0] #sys.argv[1] # "good.fa"
output_file = snakemake.output[0] #"read_embedding.txt"
output_dir = output_file #os.path.dirname(output_file) # strip the dataset.pb filename and get the output dir
k = snakemake.params[0]
encoding = {"A":1,"C":2,"G":3,"T":4,"N":0}
max_sz = 1*10**7


kmers = np.zeros((max_sz,k,5),dtype='float32')
labels = np.zeros((max_sz),dtype='float32')

print(len(kmers))
print(kmers[0])
r = 0
f1 = open(input_file)
for line in f1:
#while (line1 := f1.readline()) and (line2 := f2.readline()):
    if line[0] == ">":
        continue
    kmer = line.split()[0].upper()
    for i,n in enumerate(kmer):
        kmers[r][i][encoding[n]] = 1
    labels[r] = 0
    r += 1
    if r >= max_sz:
        break
    kmer = line.strip().split()[1].upper()
    for i,n in enumerate(kmer):
        kmers[r][i][encoding[n]] = 1
    labels[r] = 1
    r += 1
    if r >= max_sz:
        break
print(r)
kmers = kmers[:r]
labels = labels[:r]

print(kmers[0])
data = tf.data.Dataset.from_tensor_slices((kmers, labels))
data.save(output_dir)
