import re
good = open(snakemake.output[0],"w")
bad = open(snakemake.output[1],"w")
i = 0
rgx = re.compile("[0-9M]")
with open(snakemake.input[0]) as f:
    for line in f:
        i += 1
        sl = line.strip().split("\t")
        sl = [v for v in sl if v]
        if len(sl) < 9:
            continue
        one = rgx.sub("",sl[1]) and not rgx.sub("",sl[2])
        two = not rgx.sub("",sl[1]) and rgx.sub("",sl[2])
        if one:
            good.write(sl[-2]+"\n")
            bad.write(sl[-1]+"\n")
        elif two:
            good.write(sl[-1]+"\n")
            bad.write(sl[-2]+"\n")
good.close()
bad.close()
