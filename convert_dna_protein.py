import sys
import fun as fn

dnaFile = sys.argv[1]

dnaF = open(dnaFile,"r")
raw = dnaF.read()
header = raw.split("\n")[0]
dna = "".join(raw.split("\n")[1:]).strip()
dnaF.close()



pr = open(dnaFile.split(".")[0]+".protein","w")
pr.write(header + "\n")
pr.write(fn.getProteinSeq(dna))
pr.close()