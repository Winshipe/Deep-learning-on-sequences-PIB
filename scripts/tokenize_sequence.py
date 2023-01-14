import sys
import tensorflow as tf
import numpy as np
import os

max_sz = 10**7

infile = snakemake.input[0]#"good_overlaps.txt" #sys.argv[1]
vocab_file = snakemake.input[1]#"bad_overlaps.txt"
#vocab_file = snakemake.input[2]#"5mer.vocab.txt"
k = int(snakemake.params[0])
outfile = os.path.dirname(snakemake.output[0]) #"tokenized_refined_good.tf" #sys.argv[2]
#line_size = 150
n_ctx = int((k-1)//2) # the number of nts either side to consider, eg (AA[G]AA)
#kmers = np.zeros((max_sz, line_size),dtype='int32')
#labels = np.zeros(max_sz, dtype='int32')
"""
We need to interleave the 2 datasets s.t. the labels are 01010101 and so on otherwise
the training, val and test sets will be unbalanced and we'll just teach it to call 1 class
"""

kmers = np.zeros((max_sz,17),dtype='int32')
labels = np.zeros(max_sz, dtype='int32')

tokens = dict()
c = 2
#for i, path in enumerate(["5mer.txt","4mer.txt","3mer.txt"]):
with open(vocab_file) as f:
    for line in f:
        root = line.split(" ")[0]
        lr = len(root)
        tokens[root] = len(tokens)
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
with open(infile) as f:
    for line in f:
        l = len(line.strip())
        good, bad = line.strip().split("\t")
        #if l < 50:
        #    continue
        padded_g = "#"*n_ctx + good + "#" * n_ctx
        padded_b = "#"*n_ctx + bad + "#" * n_ctx
#        for c in range(min([l,150])):
        for c in range(len(good)):
            tkn_g = 1 # 1 = OOV
            try:
                tkn_g = tokens[padded_g[c:(c+2*n_ctx+1)]]
            except:
                pass
            kmers[r][c] = tkn_g
            labels[r] = 0
        r += 1
        if max_sz <= r:
            break
        for c in range(len(bad)):
            tkn_b = 1 # 1 = OOV
            try:
                tkn_b = tokens[padded_b[c:(c+2*n_ctx+1)]]
            except:
                pass
            labels[r] = 1
            kmers[r][c] = tkn_b
        r += 1
        if max_sz <= r:
            break

del tokens
print(kmers.shape)
data = tf.data.Dataset.from_tensor_slices((kmers, labels))
#data = tf.data.Dataset.from_tensor_slices((kmers_g, labels_g))
#data = data.concatenate(tf.data.Dataset.from_tensor_slices((kmers_b, labels_b)))
#data = data.interleave(lambda x: tf.data.Dataset.from_tensor_slices(x),cycle_length=4, block_length=1)
data.save(outfile)
