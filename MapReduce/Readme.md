K-means Clustering using Hadoop MapReduce

1. Place the input data file into the Hadoop Distributed File System.
2. Make sure this file is in the folder "kmeans/data/".
3. Execute the command "./run.sh arg1 arg2 arg3 arg4".


arg1 - The filename(dataset) on which we need to run K-means Clustering.
arg2 - K, which specifies the number of clusters the dataset need to be clustered into
arg3 - K space separated ids as a string, which specifies the initial cluster centroids (or) an empty string to randomly initialize cluster centroids.
arg4 - N, number of iterations for which clustering is to be performed.