import fun
import re
import sys

def chunkstring(string):
    return list(string[0+i:3+i] for i in range(0, len(string), 3))

def getHostDNA(intervals,dna):
    new_intervals = [(intervals[0] * 3), (intervals[1] * 3)]
    mapped_dna = dna[new_intervals[0]:new_intervals[1]]
    return mapped_dna
    

def getHostIntervals(match,sequence):
    r_start = match[1]
    r_end = match[2]
    l = match[0]
    seq = match[3]
    x = re.finditer(seq, sequence)
    interval = []
    for i in x:
        interval = [i.start(), i.end(), r_start, r_end]
    return interval

# ACG,ATG,AGA,TCG

virus_dna_file = sys.argv[1]
host_dna_file = sys.argv[2]


v = open(virus_dna_file, "r")
virus = v.read()
v.close()
vheader = virus.split("\n")[0]
v_dna = "".join(virus.split("\n")[1:]).strip()
v_protein = fun.getProteinSeq(v_dna)

p = open(host_dna_file, "r")
host = p.read()
p.close()
pheader = host.split("\n")[0]
r_dna = "".join(host.split("\n")[1:]).strip()
r_protein = fun.getProteinSeq(r_dna)


matches = fun.getMEM(r_protein,v_protein)
#print match

for _match in matches:
    intervals = getHostIntervals(_match, r_protein)  # [protein,virus]
    host_sub_dna = getHostDNA([intervals[0], intervals[1]], r_dna)
    virus_mapped = [(intervals[2] * 3), (intervals[3] * 3)]
    new_dna = v_dna[:virus_mapped[0]] + host_sub_dna + v_dna[virus_mapped[1]:]


o = open("optimized_virus.fasta","w")
o.write(vheader)
o.write("\n")
o.write(new_dna)
o.close()


# print "~~~ Host-DNA:: ", host_dna
# print "Host-Protein:", r_protein
# print "virus_protein:", v_protein
# print "--" * 20
# print ("Old-DNA: %s\nOld-Protein :%s\nNew-DNA: %s\nNewProtein: %s") % (v_dna,v_protein,new_dna,fun.getProteinSeq(new_dna))

