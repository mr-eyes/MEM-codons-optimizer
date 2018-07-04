#File1:
import re
import time
import fun
from operator import itemgetter
start_time = time.time()

def exists(match, selected):
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


# Reading Protein Sequences From Files
f = open("data/virus.protein")
virusProteinSeq = f.read()
f.close()
f = open("data/single_ref.protein")
hostProteinSeq = f.read()
f.close()

print fun.getMEM(virusProteinSeq,hostProteinSeq)


# Printing Execution Time
print("--- %s seconds ---" % (time.time() - start_time))
