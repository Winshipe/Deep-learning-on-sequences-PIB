import re
import threading
import multiprocessing
import subprocess
from operator import xor

readlock = threading.Lock()
writelock = threading.Lock()

if "snakemake" not in globals():
    import sys
    class Snakemake:
        output = ["small_sub_good.txt","small_sub_bad.txt","small_sub_unknown.txt"]
        input = ["new_overlaps.txt"]
        params = [17,4]
    snakemake = Snakemake()

good = open(snakemake.output[0],"w")
bad = open(snakemake.output[1],"w")
unknown = open(snakemake.output[2], "w")
k = snakemake.params[0]
nthreads = snakemake.params[1]
i = 0



def unroll_cigar(rolled):
    uc_rgx = re.compile("(?<=[A-Za-z])(?=[0-9])")
    chars = uc_rgx.split(rolled)
    out = []
    for mut in chars:
        mt = mut[-1] #like 150M, so get M
        n = int(mut[:-1])
        out.append(mt*n)
    return "".join(out)

def roll_cigar(c):
    if len(c) == 0:
        return ""
    ctr = 0
    prev = c[0]
    out = []
    for char in c:
        if char == prev:
            ctr += 1
        else:
            out.append(str(ctr)+prev)
            ctr = 1
        prev = char
    out.append(str(ctr)+prev)
    return "".join(out)

def extract_region(seq,  ix, sub_len, max_len):
    global k
    context = (k - sub_len) // 2
    start = ix - sub_len - context
    end = min([ix + context + 1, max_len])
    return seq[start:end]

def read_reference(ctg, start, stop):
    cmd = "samtools faidx panTro6.fa {}:{}-{} 2>/dev/null".format(ctg, start, stop)
    fasta = subprocess.check_output(cmd,shell=True).decode("ascii")
    return "".join([l.strip() for l in fasta.split("\n")[1:]])

def check_for_mismatches(i):
    global infile
    global good
    global bad
    global readlock
    global writelock
    
    while True:
        readlock.acquire()
        try:    
            line = infile.readline() ####
        finally:
            readlock.release()
        if line == '':
            break
        #print("s",i)
        sl = line.strip().split("\t")
        if len(sl) < 9:
            continue
        if "*" in sl[1] or "*" in sl[2] or "I" in sl[1] or "D" in sl[1] or "I" in sl[2] or "D" in sl[2]:
            continue
        cig1 = unroll_cigar(sl[1])
        cig2 = unroll_cigar(sl[2])
        pos = int(sl[6])
        ctg = sl[3]
        seq1 = sl[-2]
        seq2 = sl[-1]
        reference = read_reference(ctg, pos, pos + len(sl[-1])).upper()
        ctr = 0
        for ix, pair in enumerate(zip(cig1, cig2)):
            c1, c2 = pair
            if c1 == "M" and c2 == "M":
                if seq2[ix] != seq1[ix] or seq2[ix] != reference[ix]:
                    ctr += 1
                elif seq2[ix] == seq1[ix] and seq1[ix] == reference[ix]:
                    if ctr > 0:
                        ml = min([len(seq1),len(seq2), len(reference)])
                        s1 = extract_region(seq1, ix, ctr, ml)
                        s2 = extract_region(seq2, ix, ctr, ml)
                        ref = extract_region(reference, ix, ctr, ml)
                        mcig1 = extract_region(cig1, ix, ctr, ml)
                        mcig2 = extract_region(cig2, ix, ctr, ml)
                        if "S" in mcig1 or "S" in mcig2:
                            ctr = 0
                            continue
                        if len(s1) <= 10:
                            continue
                        if s1 != ref and s2 == ref:
                            writelock.acquire()
                            good.write(s2+"\n")
                            bad.write(s1+"\n")
                            writelock.release()
                        elif s1 == ref and s2 != ref:
                            writelock.acquire()
                            good.write(s1+"\n")
                            bad.write(s2+"\n")
                            writelock.release()
                        elif s1 == s2 and s1 != ref:
                            writelock.acquire()
                            unknown.write(s2 +"\t"+ref+"\n")
                            writelock.release()
                    ctr = 0    
            else:
                ctr = 0    
        #print("e",i)
    return

infile = open(snakemake.input[0])
threads = list(range(nthreads))
with multiprocessing.Pool(nthreads) as p:
    p.map(check_for_mismatches, threads)

#for i in range(nthreads):
#    threads.append(threading.Thread(target=check_for_mismatches,args=(i,)))
#for t in threads:
#    t.run()
#for t in threads:
#    t.join()

""" with open(snakemake.input[0]) as f:
    for line in f:
        sl = line.strip().split("\t")
        sl = [v for v in sl if v]
        if len(sl) < 9:
            continue
        if ("S" in sl[1] or "S" in sl[2]) and not ("*" in sl[1] or "*" in sl[2]): #initial_rgx.search(sl[1]) != None or initial_rgx.search(sl[2]) != None:
        #    print(sl)
    # ie contains a small substitution surrounded by good, ideally 40=3X50= or similar but we can narrow that down later
            cig1 = unroll_cigar(sl[1])
            cig2 = unroll_cigar(sl[2])
            ctr1 = 0
            ctr2 = 0
            for ix, pair in enumerate(zip(cig1, cig2)):
                c1, c2 = pair
                cond1 = c1 == "M" and c2 == "S"
                cond2 = c2 == "M" and c1 == "S"
                p1 = None
                p2 = None
                flag = False
                if cond1:
                    ctr1 += 1
                if cond2:
                    ctr2 += 1
                if p1 and not cond1 and c1 == "M" and c2 == "M":
                    if ctr1 <= 5:
                        flag, overlap_good, overlap_bad = extract_region(sl[-1],sl[-2], ix, ctr1)
                if p2 and not cond2 and c1 == "M" and c2 == "M":
                    if ctr2 <= 5:
                        flag, overlap_good, overlap_bad = extract_region(sl[-2],sl[-1], ix, ctr2)
                if flag:
                    good.write(overlap_good+"\n")
                    bad.write(overlap_bad+"\n")
 """
