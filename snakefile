configfile: "config.yaml"

import math
vocab_k = int(config["kmer size"])
vocab_context_start = int(math.ceil(float(vocab_k/2)))

import os

print(os.getcwd())

rule all:
    input:
        expand(\
            [\
            "workspace/mates.sam",\
            "workspace/{nn_id}.overlaps.txt",\
            "workspace/{nn_id}.goodbad.txt",\
            "workspace/{nn_id}_data/dataset_spec.pb",\
            "workspace/{nn_id}_model/saved_model.pb",\
            ],\
            nn_id = config["id"]
        )

rule extract_pairs:
    input:
        "workspace/mapped.bam"
    output:
        "workspace/mates.sam"
#    conda: "configs/conda.yaml"
    shell: # get 10% of the reads using seed=42, still probably overkill, extract only mates
        """
            samtools view -f 0x2 -s 42.10 -o {output[0]} {input[0]}
        """

rule get_overlaps:
 # exactly what it says in the title.  Extract the overlapping parts of our mates
 # 1 --------------################
 # 2               ################-----------
    input:
        "workspace/mates.sam"
    output:
        "workspace/{nn_id}.overlaps.txt"
    conda: "configs/conda.yaml"
    script:
        "scripts/extract_overlaps.py"

rule separate_good_bad:
# if one reads overlapping section is perfectly matched and one has substitutions
# write to good.txt and bad.txt respectively
    input:
        "workspace/{nn_id}.overlaps.txt"
    output:
        goodbad = "workspace/{nn_id}.goodbad.txt",
        unknown = "workspace/{nn_id}.unknown.txt"
    params:
        k = config["kmer size"],
        reference = "panTro6.fa"
    threads: 4
    conda: "configs/conda.yaml"
    resources:
        runtime="08:00:00"
    shell:
        """
        ./scripts/separate_good_bad {params.reference} {input} {output.goodbad}  {output.unknown} {params.k}        
        sort {output.goodbad} | uniq > {output.goodbad}.temp
        mv {output.goodbad}.temp {output.goodbad}
        """    
    #script:
    #    "scripts/separate_good_bad_small_subs.py" #"scripts/separate_good_bad.py"

#rule generate_vocab:
# if we're doing more traditional NLP  then we need to assign tokens to each kmer
# including padding (given by #)
# eg AGGAG -> 1 , AGGAT -> 2, #AGGA -> 3 and so on...
# so here we build a list of k-mers that we've found (those outside this list will be given an out-of-vocabulary token instead)
#    input:
#        good = "workspace/{nn_id}.good.txt",
#    output:
#        "workspace/{nn_id}.vocab.txt"
#    conda: "configs/conda.yaml"
#    threads: 16
#    shell:
#        """
#        awk '{{print ">" NR "\\n" $1 "\\n" $2 "\\n"}}' {input.good} > temp.fa  
#        for (( i = {vocab_context_start}; i < {vocab_k} +1 ; i++ ))
#            do
#                jellyfish count -m $i -s 1G -t 16 -o output temp.fa
#                jellyfish dump output -c >> {output[0]}
#            done
#        rm temp.fa output
#        mv {input.good}.temp {input.good}
# """

        
rule prepare_sequences_token:
# build a tensorflow Dataset object with tokenized goods and bads interwoven
    input:
        goodbad = "workspace/{nn_id}.goodbad.txt",
        vocab = "workspace/{nn_id}.vocab.txt" # list of possible kmers 
    output:
        "workspace/{nn_id}_data/dataset_spec.pb" # really we need the whole folder but this file is always here in TensorFlow v2.  Woe unto Google if they change it!
        #directory("workspace/{nn_id}_data") #nicer but difficult to put into rule all
    params:
        k = config["kmer size"]
    conda: "configs/conda.yaml"
    resources:
        gpus=1,
        partition="gpu"
    script:
        "scripts/seq_to_1hot_chars.py" #if doing tokens use tokenize_sequence.py # can also use seq_to_1hot_chars.py if using some sort of char based network
# method 1, kmer tokenization
# AAA -> 1
# AAT -> 2 
# and so on for any k
# 
# other approach than kmers is to proceed with only characters
# so we encode these 1-hot as below
#
#  A  G  T  C
#[[1 [0 [0 [0
#  0  0  0  1
#  0  1  0  0
#  0] 0] 1] 0]]
# use  "scripts/seq_to_1hot_chars.py"

rule train_model:
    input:
        "workspace/{nn_id}_data/dataset_spec.pb"
    output:
        "workspace/{nn_id}_model/saved_model.pb"
    params:
        model_name = "workspace/{nn_id}_model"
    resources:
        gpus=1,
        partition="gpu",
        mem=4000
    conda: "configs/conda.yaml"
    threads: 8
    script:
        "scripts/kmer_sentiments.py" # can also use cnn_accuracy.py here


