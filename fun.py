import re
from operator import itemgetter

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
    for i in xrange(len(fseq) - window_size + 1):
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
            print matches[i], " Selected!"
            lengthSum += matches[i][0]
            continue

        if not exists([start, end], selected):
            selected.append([start, end])
            # Comments Here
            print matches[i], " Selected!"
            lengthSum += matches[i][0]
        else:
            pass
            # print "Pass"

        # print "oldStart: %s - oldEnd: %s - Start: %s - End: %s" % (oldStart,oldEnd,start,end)
        # print matches[i]

    #print "--" * 70
    return lengthSum
