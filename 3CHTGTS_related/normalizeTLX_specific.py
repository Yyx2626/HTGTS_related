import sys
from random import sample, seed

#python3 normalizeTLX.py 100 result1.tlx result2.tlx result3.tlx

def main():


    tlx_dict = {}
    labels = ""
    toNormalize = int(sys.argv[1])
    
    #parse each tlx into a dictionary of lines:file
    minNumber = 100000000000
    minTlx = ""
    for tlx in sys.argv[2:]:
        tlxread = open(tlx, "r")
        labels = tlxread.readline()
        lines = []
        for line in tlxread:
            lines.append(line)
        tlx_dict[tlx] = lines
        if len(lines) < minNumber:
            minNumber = len(lines)
            minTlx = tlx

        print(tlx + " contains " + str(len(lines)) + " junctions")



    #print("Normalizing to " + minTlx + " with " + str(minNumber) + " junctions")

    for tl in tlx_dict:
        seed(1234567)   # Yyx added just for reproducibility
        randomSubset = sample(tlx_dict[tl], toNormalize)
        tlx_dict[tl] = randomSubset

    for t in tlx_dict:
        idx = t.index(".tlx")
        normFile = t[:idx]+"."+str(toNormalize)+".tlx"
        writeNormFile = open(normFile, "w")
        writeNormFile.write(labels)
        for lines in tlx_dict[t]:
            writeNormFile.write(lines)

main()
