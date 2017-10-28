Running K-Means clustering:

Note: Make sure the necessary input files are placed in this folder before running the program.

1. Invoke the program as follows - python3 K-means.py <filename> <Number of clusters> [Initial centroid list] [Iterations]
   Note that the arguments [Initial centroid list] [Iterations] are optinal and by default random cluster values are taken along with max iterations as 999.

2. The program will output the jaccard and clusters are plotted using PCA, which can be visualized for checking the performance of the algorithm.


   Examples:
   python3 K-means.py cho.txt 5                           - Missing optional arguments
   python3 K-means.py cho.txt 5 [378,55,51,35,236] 20     - With correct optional arguments
