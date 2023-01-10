import re

if "snakemake" not in globals():
    import sys
    class Snakemake:
        output = ["new_overlaps.txt"] 
        input = ["mates.sam"]
    snakemake = Snakemake()

rgx = re.compile("(?<=[A-Za-z])(?=[0-9])")
def unroll_cigar(rolled):
    chars = rgx.split(rolled)
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
olps = open(snakemake.output[0],"w")
with open(snakemake.input[0]) as f:
    i = 0
    seqs = dict()
    for line in f:
        if line[0] == "@":
            continue
        sl = line.split("\t")
        ctg = sl[2]
        cigar = sl[5]
        ix = int(sl[3])
        l = int(sl[8])
        seq = sl[9]
        name = sl[0]
        if name not in seqs:
            seqs[name] = []
        seqs[name].append((seq,ix,l,cigar))

        if len(seqs[name]) == 2:
            s1, s2 = seqs[name]
            start1 = s1[1]
            start2 = s2[1]
            end1 = start1 + abs(s1[2])
            end2 = start2 + abs(s2[2])
            o1, o2 = "",""
            cstart1 = start1
            cstart2 = start2
            cigar1 = s1[-1]
            cigar2 = s2[-1]
            if cigar1 == "*" or cigar2 == "*":
                continue
            ur_c1 = unroll_cigar(cigar1)
            ur_c2 = unroll_cigar(cigar2)
            last_s_c1 = ur_c1.rfind("S")
            first_s_c1 = ur_c1.find("S")
            last_s_c2 = ur_c2.rfind("S")
            first_s_c2 = ur_c2.find("S")

            if last_s_c1 != -1 and last_s_c1 != len(ur_c1):
                start1 = last_s_c1
            elif first_s_c1 != -1 and first_s_c1 != 0:
                end1 = first_s_c1

            if last_s_c1 != -1 and last_s_c1 != len(ur_c1):
                start2 = last_s_c1
            elif first_s_c1 != -1 and first_s_c1 != 0:
                end2 = first_s_c1    
            
            if start1 < start2:
                o1 = s1[0][start2 - start1:]
                o2 = s2[0][:len(o1)]
                cstart1 = start2
                mc1 = roll_cigar(ur_c1[start2 - start1:])
                mc2 = roll_cigar(ur_c2[:len(o1)])
            else:
                o2 = s2[0][start1 - start2:]
                o1 = s1[0][:len(o2)]
                cstart2 = start1
                mc2 = roll_cigar(ur_c2[start2 - start1:])
                mc1 = roll_cigar(ur_c1[:len(o1)])
            olps.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n"%(name, mc1,mc2,ctg, start1,start2,cstart1,o1,o2))
            del seqs[name]
olps.close()