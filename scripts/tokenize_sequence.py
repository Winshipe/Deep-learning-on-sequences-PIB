import sys
import tensorflow as tf
import numpy as np
import os

max_sz = 10**7

infile_good = snakemake.input[0]#"good_overlaps.txt" #sys.argv[1]
infile_bad = snakemake.input[1]#"bad_overlaps.txt"
vocab_file = snakemake.input[2]#"5mer.vocab.txt"
k = int(snakemake.params[0])
outfile = os.path.dirname(snakemake.output[0]) #"tokenized_refined_good.tf" #sys.argv[2]
line_size = 150
n_ctx = 2 # the number of nts either side to consider, eg (AA[G]AA)
#kmers = np.zeros((max_sz, line_size),dtype='int32')
#labels = np.zeros(max_sz, dtype='int32')
"""
We need to interleave the 2 datasets s.t. the labels are 01010101 and so on otherwise
the training, val and test sets will be unbalanced and we'll just teach it to call 1 class
"""

kmers_g = np.zeros((max_sz,line_size),dtype='int32')
labels_g = np.zeros(max_sz, dtype='int32')

tokens = dict()
c = 2
#for i, path in enumerate(["5mer.txt","4mer.txt","3mer.txt"]):
with open(vocab_file) as f:
    for line in f:
        root = line.split(" ")[0]
        lr = len(root)
        for i in range(lr):
            nroot = root[:i] + "N" + root[i+1:]
            nfirst = "#"*(k-lr) + nroot
            nsecond = nroot + "#"*(k-lr)
            tokens[nfirst] = len(tokens) + 2
            tokens[nsecond] = len(tokens) + 2
with open("tokens.txt","w") as f:
    for kmer, tkn in tokens.items():
        f.write("%s\t%d\n"%(kmer, tkn))
        
print(len(tokens))
r = 0 
with open(infile_good) as f:
    for line in f:
        l = len(line.strip())
        if l < 50:
            continue
        padded = "#"*n_ctx + line.strip() + "#" * n_ctx
        for c in range(min([l,150])):
            tkn = 1 # 1 = OOV
            try:
                tkn = tokens[padded[c:(c+2*n_ctx+1)]]
            except:
                pass
            kmers_g[r][c] = tkn
            labels_g[r] = 0
        r += 1
        if max_sz <= r:
            break
kmers_g = kmers_g[:r]
labels_g = labels_g[:r]
kmers_b = np.zeros((max_sz,line_size),dtype='int32')
labels_b = np.zeros(max_sz, dtype='int32')
r = 0
with open(infile_bad) as f:
    for line in f:
        l = len(line.strip())
        if l < 50:
            continue
        padded = "#"*n_ctx + line.strip() + "#" * n_ctx
        for c in range(min([l,150])):
            tkn = 0
            try:
                tkn = tokens[padded[c:(c+2*n_ctx+1)]]
            except:
                pass
            kmers_b[r][c] = tkn
            labels_b[r] = 1
        r += 1
        if max_sz <= r:
            break
newlen = min([len(kmers_b),len(kmers_g)])
kmers = np.zeros((newlen*2,line_size),dtype='int32')
labels = np.zeros(newlen*2,dtype='int32')
kmers[0::2] = kmers_g[:newlen]
labels[0::2] = labels_g[:newlen]
kmers[1::2] = kmers_b[:newlen]
labels[1::2] = labels_b[:newlen]
del kmers_b
del kmers_g
del labels_b
del labels_g
del tokens
print(kmers.shape)
data = tf.data.Dataset.from_tensor_slices((kmers, labels))
#data = tf.data.Dataset.from_tensor_slices((kmers_g, labels_g))
#data = data.concatenate(tf.data.Dataset.from_tensor_slices((kmers_b, labels_b)))
#data = data.interleave(lambda x: tf.data.Dataset.from_tensor_slices(x),cycle_length=4, block_length=1)
data.save(outfile)
