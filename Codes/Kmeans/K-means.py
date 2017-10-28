import random
import ast
import numpy as np
import sys

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import TruncatedSVD

class Point:

    def __init__(self, dataPoints, trueCluster):
        self.x = list(dataPoints)
        self.clusterNumber = None
        self.trueCluster = trueCluster

    def setCluster(self, cluster):
        self.clusterNumber = cluster


    def getCluster(self):
        return self.clusterNumber

    def getPointValue(self):
        return self.x


def kmeans(filename, clusters, clusterpoints, iterations):
    np.set_printoptions(suppress=True)

    geneData = np.genfromtxt(filename,delimiter='\t', dtype=float)
    trueClusters = geneData[:, [1]]

    geneData =  np.delete(geneData, [0,1], axis=1)

    trueClusters = open(filename, "r").readlines()
    trueClusters = [x.split("\t")[1] for x in trueClusters]
    trueClusters = np.array(trueClusters)


    if(clusterpoints == None or len(clusterpoints) == 0):
        sampleClusters = random.sample(range(0, geneData.shape[0]), clusters)
    else:
        sampleClusters = clusterpoints

    # mean = np.zeros((5,16))
    # counts = [0]*5

    # for i in range(0, len(trueClusters)):
    #     ind = int(trueClusters[i])-1
    #     mean[ind] += geneData[i]
    #     counts[ind] += 1

    # for i in range(0,len(counts)):
    #     mean[i] = mean[i]/counts[i]

    clusterPoints = []
    for cluster in sampleClusters:
        clusterPoints.append(geneData[cluster-1])
    #
    clusterPoints = np.array(clusterPoints)
    #

    dataPoints = []
    for i in range(0, len(geneData)):
        dataPoints.append(Point(geneData[i], trueClusters[i]))

    needsClutering = 1

    while(needsClutering and iterations != 0):
        iterations -= 1
        needsClutering = performCluster(dataPoints, clusterPoints)
        clusterPoints = reCalculateCentroids(dataPoints, np.shape(clusterPoints)[0])


    groundTruthMatrix = np.zeros(((np.shape(geneData))[0],(np.shape(geneData))[0]))
    generatedClusterMatrix = np.zeros(((np.shape(geneData))[0],(np.shape(geneData))[0]))

    for i in range(0, len(trueClusters)):
        for j in range(0, len(trueClusters)):
            if(trueClusters[i] == trueClusters[j]):
                groundTruthMatrix[i][j] = 1

    for i in range(0, len(dataPoints)):
        for j in range(0, len(dataPoints)):
            if(dataPoints[i].getCluster() == dataPoints[j].getCluster()):
                generatedClusterMatrix[i][j] = 1


    countones = 0;
    countoneandZeros = 0

    for i in range(0, len(trueClusters)):
        for j in range(0, len(trueClusters)):
            if(groundTruthMatrix[i][j] == 1 and generatedClusterMatrix[i][j] == 1):
                countones += 1
            elif(groundTruthMatrix[i][j] !=  generatedClusterMatrix[i][j]):
                countoneandZeros += 1

    print(countones/(countoneandZeros+countones))



    countones = 0
    countoneandZeros = 0
    for i in range(len(dataPoints)):
        for j in range(len(dataPoints)):
            tc1 = dataPoints[i].trueCluster
            tc2 = dataPoints[j].trueCluster
            cn1 = dataPoints[i].clusterNumber
            cn2 = dataPoints[j].clusterNumber

            if tc1 == tc2 and cn1 == cn2:
                countones+=1
            elif tc1 == tc2 or cn1 == cn2:
                countoneandZeros+=1

    # print(countones)
    # print(countoneandZeros)

    assignedClusters = [point.clusterNumber for point in dataPoints]
    svdDim = TruncatedSVD(n_components=2).fit_transform(geneData).T
    plot(svdDim, "K-Means on " + filename + "trueClusters", trueClusters.T)
    plot(svdDim, "K-Means on " + filename + "assignedClusters", assignedClusters)
    plt.show()

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


def reCalculateCentroids(dataPoints, clusterSize):

    clusterOccurence = np.zeros(clusterSize);
    pointDimension = np.shape(dataPoints[0].getPointValue())[0]
    currentClusterSum = np.zeros((clusterSize,
                                 pointDimension))

    for i in range(0, len(dataPoints)):
        clusterOccurence[dataPoints[i].getCluster()] += 1
        currentClusterSum[dataPoints[i].getCluster()] += dataPoints[i].getPointValue()

    for i in range(0, len(clusterOccurence)):
        currentClusterSum[i] = currentClusterSum[i] / clusterOccurence[i]

    return currentClusterSum



def performCluster(dataPoints, clusterPoints):
    needsClustering = 0

    for currentPoint in dataPoints:
        minValue = float("inf")
        index = 0;
        for j in range(0, len(clusterPoints)):
            currentDistance = distance(currentPoint.getPointValue(), clusterPoints[j])
            if currentDistance < minValue:
                minValue = currentDistance
                index = j

        if(currentPoint.getCluster() != index):
            needsClustering = 1
            currentPoint.setCluster(index)

    return needsClustering


def distance(point1, point2):
    dist = point1 - point2
    dist = np.square(dist)
    dist = np.sum(dist)
    return np.sqrt(dist)





filename      = sys.argv[1]
clusters      = int(sys.argv[2])

if(len(sys.argv) > 3):
    clusterPoints = ast.literal_eval(sys.argv[3])
else:
    clusterPoints = None

if(len(sys.argv) > 4):
    iterations    = int(sys.argv[4])
else:
    iterations    = 999



kmeans(filename, clusters, clusterPoints, iterations)


