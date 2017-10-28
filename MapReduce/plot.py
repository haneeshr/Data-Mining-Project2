import random
import numpy as np
import sys

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import TruncatedSVD

def read_file(filename):
	geneData = np.genfromtxt(filename,delimiter='\t', dtype=float)
	trueClusters = geneData[:, [1]]

	geneData =  np.delete(geneData, [0,1], axis=1)

	trueClusters = open(filename, "r").readlines()
	trueClusters = [x.split("\t")[1] for x in trueClusters]
	trueClusters = np.array(trueClusters)
	return geneData, trueClusters


def readPartFile():
	filename = "output/output/part-r-00000"

	clusterIds = np.genfromtxt(filename,delimiter='\t')

	assignedClusters = [0]*len(clusterIds)


	for point in clusterIds:
		assignedClusters[int(point[0]-1)] = point[1]

	return assignedClusters

def calculateJaccardCoeff(assignedClusters, trueClusters):

	onescount = 0
	onesandzeroscount = 0
	for i in range(len(assignedClusters)):
		for j in range(len(assignedClusters)):
			tc1 = trueClusters[i]
			tc2 = trueClusters[j]
			nc1 = assignedClusters[i]
			nc2 = assignedClusters[j]

			if tc1==tc2 and nc1==nc2:
				onescount+=1
			elif tc1==tc2 or nc1==nc2:
				onesandzeroscount+=1

	return float(onescount)/float(onesandzeroscount+onescount)

def plot(reeducedDimensions, title, diseases):
    df = pd.DataFrame(dict(x=reeducedDimensions[0], y=reeducedDimensions[1], label=diseases))

    groups = df.groupby('label')
    fig, ax = plt.subplots()
    for name, group in groups:
        ax.plot(group.x, group.y, marker='o', linestyle='', label=name)
    ax.legend()
    plt.suptitle(title)
    plt.xlabel('Principle Component 1')
    plt.ylabel('Principle Component 2')


def main(filename):
	geneData, trueClusters = read_file(filename)
	assignedClusters = readPartFile()

	jaccard_coeffecient = calculateJaccardCoeff(assignedClusters, trueClusters)
	print "Jaccard Coeffecient for MapReduce K-means on " + filename.split("/")[1] + "\t" + str(jaccard_coeffecient)

	svdDim = TruncatedSVD(n_components=2).fit_transform(geneData).T

	plot(svdDim, "Kmeans on " + filename + "trueClusters", trueClusters)
	plot(svdDim, "Kmeans on " + filename + "assignedClusters", assignedClusters)
	plt.show()



main('../' + sys.argv[1])

