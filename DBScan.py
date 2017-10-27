import random
import numpy as np
import sys

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import TruncatedSVD

dataPoints = []

class Point:

    def __init__(self, dataPoints, trueCluster):
        self.x = np.array(dataPoints)
        self.clusterNumber = None
        self.trueCluster = trueCluster
        self.visited = False

    def setCluster(self, cluster):
        self.clusterNumber = cluster

    def getCluster(self):
        return self.clusterNumber

    def getPointValue(self):
        return self.x

    def setVisited(self):
        self.visited = True

    def getVisited(self):
        return self.visited


def dbscan(filename, eps, minpoints):
    np.set_printoptions(suppress=True)

    clusternumber = 1

    geneData = np.genfromtxt(filename,delimiter='\t', dtype=float)
    trueClusters = geneData[:, [1]]

    geneData =  np.delete(geneData, [0, 1], axis=1)

    trueClusters = open(filename, "r").readlines()
    trueClusters = [x.split("\t")[1] for x in trueClusters]
    trueClusters = np.array(trueClusters)

    global dataPoints

    for i in range(0, len(geneData)):
        dataPoints.append(Point(geneData[i], trueClusters[i]))


    for i in range(0, len(dataPoints)):

        if(dataPoints[i].getVisited() == True): continue
        dataPoints[i].setVisited()

        neighbourpoints = findNeighbours(dataPoints[i], eps)

        if(len(neighbourpoints) < minpoints):
            dataPoints[i].setCluster(-1) #Outliers
        else:
            dataPoints[i].setCluster(clusternumber)
            expandCluster(neighbourpoints, clusternumber, eps, minpoints)
            clusternumber += 1


    countones = 0
    countoneandZeros = 0
    for i in range(len(dataPoints)):
        for j in range(len(dataPoints)):
            tc1 = dataPoints[i].trueCluster
            tc2 = dataPoints[j].trueCluster
            cn1 = dataPoints[i].clusterNumber
            cn2 = dataPoints[j].clusterNumber

            if tc1 == tc2 and cn1 == cn2:
                countones += 1
            elif tc1 == tc2 or cn1 == cn2:
                countoneandZeros += 1

    print(float(countones) / float(countones + countoneandZeros))

    clusterNumberset = set()
    for i in range(0, len(dataPoints)):
        clusterNumberset.add(dataPoints[i].getCluster())

    print(clusterNumberset)

    assignedClusters = [point.clusterNumber for point in dataPoints]
    svdDim = TruncatedSVD(n_components=2).fit_transform(geneData).T

    plot(svdDim, "DBScan clustering" + filename + "trueClusters", trueClusters.T)
    plot(svdDim, "DBscan clustering" + filename + "assignedClusters", assignedClusters)
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

def findNeighbours(point, eps):



    neighbourPoints = []

    for i in range(0, len(dataPoints)):
        if distance(dataPoints[i].getPointValue(), point.getPointValue()) <= eps:
            neighbourPoints.append(dataPoints[i])


    return neighbourPoints



def expandCluster(neighbourpoints, clusternumber, eps, minpoints):

    for i in range(0, len(neighbourpoints)):

        if(neighbourpoints[i].getVisited() == False):

            neighbourpoints[i].setVisited()
            neighbourpoints[i].setCluster(clusternumber)
            currentneighbours = findNeighbours(neighbourpoints[i], eps)
            if(len(currentneighbours) >= minpoints):
                expandCluster(currentneighbours, clusternumber, eps, minpoints)

        if(neighbourpoints[i].getCluster() == None):
            neighbourpoints[i].setVisited()
            neighbourpoints[i].setCluster(clusternumber)


def distance(point1, point2):
    dist = point1 - point2
    dist = np.square(dist)
    dist = np.sum(dist)
    # return dist
    return np.sqrt(dist)


dbscan('iyer.txt', 1.03, 4)
