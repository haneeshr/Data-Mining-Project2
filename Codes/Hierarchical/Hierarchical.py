import random
import numpy as np
import sys


import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import TruncatedSVD

class Point:

    def __init__(self, dataPoints, trueCluster):
        self.x = np.array((dataPoints))
        self.trueCluster = trueCluster
        self.clusterNumber = None

    def getTrueCluster(self):
        return self.trueCluster

    def getPointValue(self):
        return self.x

class Cluster:

    def __init__(self, point, clusterNumber):
        self.clusterNumber = clusterNumber
        self.points = [point]

    def getClusterValue(self):
        return self.clusterNumber

    def getPoints(self):
        return self.points

    def addPoints(self, newPoints):
        for point in newPoints:
            self.points.append(point)


def merge(cluster1, cluster2):
    cluster1.addPoints(cluster2.getPoints())


def hierarchical(filename, k):
    np.set_printoptions(suppress=True)

    geneData = np.genfromtxt(filename,delimiter='\t', dtype=float)
    trueClusters = geneData[:, [1]]

    geneData =  np.delete(geneData, [0,1], axis=1)

    trueClusters = open(filename, "r").readlines()
    trueClusters = [x.split("\t")[1] for x in trueClusters]
    trueClusters = np.array(trueClusters)

    sampleClusters = random.sample(range(0, geneData.shape[0]), k)

    dataPoints = []
    for i in range(0, len(geneData)):
        dataPoints.append(Point(geneData[i], trueClusters[i]))

    #Initial clustering
    clusters = []
    for i in range(0, len(geneData)):
        clusters.append(Cluster(dataPoints[i],i))


    clusteredMatrix = np.zeros((np.shape(geneData)[0], np.shape(geneData)[0]))

    for i in range(0, len(clusteredMatrix)):
        for j in range(0, len(clusteredMatrix)):
            if i == j: clusteredMatrix[i][j] = float("inf")
            else:
                #Update the distance between the points
                clusteredMatrix[i][j] = distance(dataPoints[i].getPointValue(), dataPoints[j].getPointValue())



    while(len(clusteredMatrix) != k):
        deleteIndexes = np.array(np.unravel_index(np.argmin(clusteredMatrix), clusteredMatrix.shape))
        row = deleteIndexes[0]
        col = deleteIndexes[1]
        # print(row, col)
        for i in range(0, len(clusteredMatrix)):
            if i==row:
                continue

            clusteredMatrix[row][i] = min(clusteredMatrix[row][i], clusteredMatrix[i][col])
            clusteredMatrix[i][row] = clusteredMatrix[row][i]

        clusteredMatrix = np.delete(clusteredMatrix, col, axis=0)
        clusteredMatrix = np.delete(clusteredMatrix, col, axis=1)

        merge(clusters[row], clusters[col])
        clusters.remove(clusters[col])


    print(calculateJaccardCoeff(clusters))

    assignedClusters = [point.clusterNumber for point in dataPoints]
    svdDim = TruncatedSVD(n_components=2).fit_transform(geneData).T
    # print(np.shape(trueClusters))
    plot(svdDim, "HAC clustering " + filename + "trueClusters", trueClusters.T)
    plot(svdDim, "HAC clustering " + filename + "assignedClusters", assignedClusters)
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

def calculateJaccardCoeff(clusters):

    onescount = 0
    onesandzeroscount = 0
    for cluster in clusters:
        for point in cluster.getPoints():
            point.clusterNumber = cluster.clusterNumber
            for otherCluster in clusters:
                for otherPoint in otherCluster.getPoints():
                    tc1 = point.trueCluster
                    tc2 = otherPoint.trueCluster
                    cn1 = cluster.clusterNumber
                    cn2 = otherCluster.clusterNumber
                    if tc1 == tc2 and cn1 == cn2:
                        onescount += 1
                    elif tc1 == tc2 or cn1 == cn2:
                        onesandzeroscount += 1

    return float(onescount)/float(onescount+onesandzeroscount)
                    

def distance(point1, point2):
    dist = (point1 - point2)
    dist = np.square(dist)
    dist = np.sum(dist)
    return np.sqrt(dist)


filename      = sys.argv[1]
clusters      = int(sys.argv[2])

hierarchical(filename, clusters)
