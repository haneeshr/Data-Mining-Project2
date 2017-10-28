#! /bin/sh

filename=$1
k=$2
cluster_ids=$3
num_iteration=$4

python initial_centroids.py $filename $cluster_ids

inputfile=kmeans/data/$filename
rm -R output/
mkdir output
hadoop com.sun.tools.javac.Main Kmeans.java
jar cf kmeans.jar Kmeans*.class
hadoop jar kmeans.jar Kmeans $inputfile kmeans/output/ $k $num_iteration
rm Kmeans*.class
rm kmeans.jar
hdfs dfs -get kmeans/output/ output/
hdfs dfs -rm -r kmeans/output
python plot.py $filename