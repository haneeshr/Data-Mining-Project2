import numpy as np
import sys

filename = sys.argv[1]
filename = "../"+filename
genedata = open(filename, "r")
lines = genedata.readlines()

centroids = open("centroids.txt", 'w')
for i in range(2, len(sys.argv)):
	id = int(sys.argv[i])
	line = lines[id-1]
	line = line.split("\t")
	line = line[2:]
	line = "\t".join(line)
	centroids.write(line)