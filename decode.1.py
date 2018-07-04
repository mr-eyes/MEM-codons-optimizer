import fun
import re


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

v_dna = "CCCACGATGAGATCGCCC"
v_protein = "PTMRSP"

r_dna = "GCAGCAACCATGCGCAGCGCA"
r_protein = "AATMRSA"

match = fun.getMEM(r_protein,v_protein)
#print match

intervals = getHostIntervals(match[-1], r_protein) #[protein,virus]

print intervals

host_dna = getHostDNA([intervals[0],intervals[1]] , r_dna)

virus_mapped = [(intervals[2] * 3), (intervals[3] * 3)]

new_dna = v_dna[:virus_mapped[0]] + host_dna + v_dna[virus_mapped[1]:]




print "~~~ Host-DNA:: ", host_dna
print "Host-Protein:", r_protein
print "virus_protein:", v_protein
print "--" * 20
print ("Old-DNA: %s\nOld-Protein :%s\nNew-DNA: %s\nNewProtein: %s") % (v_dna,v_protein,new_dna,fun.getProteinSeq(new_dna))
