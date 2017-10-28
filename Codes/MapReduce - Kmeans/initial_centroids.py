import numpy as np
import sys
import random

filename = sys.argv[1]
genedata = open(filename, "r")
lines = genedata.readlines()

centroids = open("centroids.txt", 'w')

k = int(sys.argv[2])

ids = random.sample(range(0, len(lines)), k)
if(len(sys.argv) == k+3):
	ids = sys.argv[3:]

for id in ids:
	id = int(id)
	line = lines[id-1]
	line = line.split("\t")
	line = line[2:]
	line = "\t".join(line)
	centroids.write(line)