import random
import numpy as np
import sys

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


def kmeans(filename, clusters):
    np.set_printoptions(suppress=True)

    geneData = np.genfromtxt(filename,delimiter='\t', dtype=float)
    trueClusters = geneData[:, [0,1]]

    geneData =  np.delete(geneData, [0,1], axis=1)

    sampleClusters = random.sample(range(0, geneData.shape[0]), clusters)

    clusterPoints = []
    for cluster in sampleClusters:
        clusterPoints.append(geneData[cluster])

    clusterPoints = np.array(clusterPoints)

    dataPoints = []
    for i in range(0, len(geneData)):
        dataPoints.append(Point(geneData[i], trueClusters[i]))

    needsClutering = 1;

    while(needsClutering):
        needsClutering = performCluster(dataPoints, clusterPoints)
        clusterPoints = reCalculateCentroids(dataPoints, np.shape(clusterPoints)[0])

    # for point in dataPoints:
    #     print(point.getCluster())


def reCalculateCentroids(dataPoints, clusterSize):

    clusterOccurence = np.zeros(clusterSize);
    # print(np.shape(dataPoints[0].getPointValue())[0])
    pointDimension = np.shape(dataPoints[0].getPointValue())[0]
    currentClusterSum = np.zeros((clusterSize,
                                 pointDimension))

    for i in range(0,len(dataPoints)):
        clusterOccurence[dataPoints[i].getCluster()] += 1
        currentClusterSum[dataPoints[i].getCluster()] += dataPoints[i].getPointValue()

    for i in range(0,len(clusterOccurence)):
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
    dist = dist*dist
    dist =  sum(dist)
    # dist = dist^0.5
    return dist



# filename = sys.argv[1]
# clusters = int(sys.argv[2])
kmeans('cho.txt', 5)


