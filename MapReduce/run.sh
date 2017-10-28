#! /bin/sh

filename=$1
k=$2
cluster_ids=$3
num_iteration=$4

python initial_centroids.py $filename $k $cluster_ids
inputfile=kmeans/data/$filename
rm -R output/
mkdir output
hadoop com.sun.tools.javac.Main Kmeans.java
jar cf kmeans.jar Kmeans*.class
start_time=`date +%s`
hadoop jar kmeans.jar Kmeans $inputfile kmeans/output/ $k $num_iteration
hdfs dfs -get kmeans/output/ output/
hdfs dfs -rm -r kmeans/output
python plot.py $filename
end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
rm Kmeans*.class
rm kmeans.jar