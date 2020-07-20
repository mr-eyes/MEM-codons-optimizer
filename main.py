import re
import sys
from operator import itemgetter


def chunkstring(string):
    return list(string[0+i:3+i] for i in range(0, len(string), 3))


def getHostDNA(intervals, dna):
    new_intervals = [(intervals[0] * 3), (intervals[1] * 3)]
    mapped_dna = dna[new_intervals[0]:new_intervals[1]]
    return mapped_dna

def getProteinSeq(seq):
    codonUsage = {"A": ["GCG", "GCA", "GCT", "GCC"], "C": ["TGT", "TGC"], "D": ["GAT", "GAC"], "E": ["GAG", "GAA"], "F": ["TTT", "TTC"], "G": ["GGG", "GGA", "GGT", "GGC"], "H": ["CAT", "CAC"], "I": ["ATA", "ATT", "ATC"], "K": ["AAG", "AAA"], "L": ["TTG", "TTA", "CTG", "CTA", "CTT", "CTC"], "M": ["ATG"], "N": [
        "AAT", "AAC"], "P": ["CCG", "CCA", "CCT", "CCC"], "Q": ["CAG", "CAA"], "R": ["AGG", "AGA", "CGG", "CGA", "CGT", "CGC"], "S": ["AGT", "AGC", "TCG", "TCA", "TCT", "TCC"], "T": ["ACG", "ACA", "ACT", "ACC"], "V": ["GTG", "GTA", "GTT", "GTC"], "W": ["TGG"], "Y": ["TAT", "TAC"], "*": ["TGA", "TAG", "TAA"]}
    proteinSeq = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        for amino in codonUsage:
            if (codon in codonUsage[amino]):
                proteinSeq += amino
    return proteinSeq


def exists(match, selected):
    #print "match:",match
    #print "selected:", selected
    m_start = match[0]
    m_end = match[1]
    for s in selected:
        s_start = s[0]
        s_end = s[1]
        if(m_start >= s_start and m_end <= s_end):
            return True
    return False

# Windowing Function


def window(fseq, window_size=5):
    for i in range(len(fseq) - window_size + 1):
        yield fseq[i:i+window_size]


def getMEM(virusProteinSeq, hostProteinSeq):
    if ">" in virusProteinSeq:
        virusProteinSeq = "".join(virusProteinSeq.split("\n")[1:]).strip()
    if ">" in hostProteinSeq:
        hostProteinSeq = "".join(hostProteinSeq.split("\n")[1:]).strip()

    
    results = []
    #resultsFile = []
    lengthSum = 0
    windowLength = min([len(hostProteinSeq), len(virusProteinSeq)])
    #windowLength = len(virusProteinSeq)
    
    for i in range(windowLength, 3, -1):
        w = 0
        for motif in window(virusProteinSeq, i):
            results.append([i, w, w + i, motif])
            #resultsFile.append("Windows %d indexes %d-%d Motif: %s" % (i, w, w + i - 1, motif))
            w += 1


     
    # # [Optional] Writing windows to file
    # f = open("windows.txt", 'w+')
    # f.write("\n".join(resultsFile))
    # f.close()

    # Empty List of matches
    matches = []

    # Iterating through windowing results for matching
    for result in results:
        motif = result[3]  # Result->[windowSize, startIndex, EndIndex, Motif]
        # print motif
        for m in re.finditer(motif, hostProteinSeq):
            # print (motif, m.start(), m.end()-1)
            # print "--" * 30
            matches.append([len(motif), m.start(), m.end(), motif])
            #matches.append([m.start(),len(motif))

    # Sorting The Matches
    sortedMatches = sorted(matches, key=itemgetter(0), reverse=True)
    matches = sortedMatches
    return matches

    selected = []

    for i in range(len(matches)):
        length = matches[i][0]
        start = matches[i][1]
        end = matches[i][2]
        word = matches[i][3]
        if i == 0:
            selected.append([start, end])
            # Comments Here
            print (matches[i], " Selected!")
            lengthSum += matches[i][0]
            continue

        if not exists([start, end], selected):
            selected.append([start, end])
            # Comments Here
            print (matches[i], " Selected!")
            lengthSum += matches[i][0]
        else:
            pass
            # print "Pass"

        # print "oldStart: %s - oldEnd: %s - Start: %s - End: %s" % (oldStart,oldEnd,start,end)
        # print matches[i]

    #print "--" * 70
    return lengthSum


def getHostIntervals(match, sequence):
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
v_protein = getProteinSeq(v_dna)

p = open(host_dna_file, "r")
host = p.read()
p.close()
pheader = host.split("\n")[0]
r_dna = "".join(host.split("\n")[1:]).strip()
r_protein = getProteinSeq(r_dna)


matches = getMEM(r_protein, v_protein)
#print match

for _match in matches:
    intervals = getHostIntervals(_match, r_protein)  # [protein,virus]
    host_sub_dna = getHostDNA([intervals[0], intervals[1]], r_dna)
    virus_mapped = [(intervals[2] * 3), (intervals[3] * 3)]
    new_dna = v_dna[:virus_mapped[0]] + host_sub_dna + v_dna[virus_mapped[1]:]


o = open("optimized_virus.fasta", "w")
o.write(vheader + "_optimized")
o.write("\n")
o.write(new_dna)
o.close()


# print "~~~ Host-DNA:: ", host_dna
# print "Host-Protein:", r_protein
# print "virus_protein:", v_protein
# print "--" * 20
# print ("Old-DNA: %s\nOld-Protein :%s\nNew-DNA: %s\nNewProtein: %s") % (v_dna,v_protein,new_dna,getProteinSeq(new_dna))
