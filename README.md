# Data Splitting Methods Based on Data Distribution Consistency  
### Systematic: 
Systematic sampling is another deterministic approach in which every kth observation is sampled. If the data are ordered in some way, this implicitly generates a stratified sample, with stratification on the ordinal variable. One
approach is to sort the data along the output variable dimension to obtain a representative sample of the output variable distribution (Baxter, Stanley, Zhang, & Smith, 2000). This approach is easy to implement, as it assumes that the output variable can be mapped to a unique input state. However, this assumption may not hold in multivariate datasets where multiple input states might give rise to the same output, where the method cannot ensure that representative input–output combinations will be sampled, since only the output variable is considered. 
### DUPLEX: 
CADEX, or Kennard–Stone sampling (Kennard & Stone, 1969), is one of the earliest algorithms designed for data splitting. The approach iteratively draws samples based on distance, selecting points farthest away from those already included in the sample, and ensures maximum coverage of the data. An improved version called DUPLEX was proposed by Snee (1977); it is used widely in the field of chemometrics, including several ANN applications (Despagne & Massart, 1998; Sprevak, Azuaje, & Wang, 2004). However, the computational complexity of this algorithm may prohibit its use on large datasets.
### SBSS: 
SOM-based stratified sampling The SBSS approach involves two steps (Bowden et al., 2002). In the first step, the data are partitioned into K strata (clusters) using a selforganizing map (SOM; Kohonen, 1990), which considers the distances between data points. In the second step, data for the calibration and evaluation subsets are obtained by sampling from each of these strata. For SBSS-P, the sampling is done in proportion to the number of samples in each stratum. For SBSS-N, the sample allocation is increased for strata that contain a larger number of data points, or where the data points within a stratum have a larger variance (May et al., 2010).  
### SOMPLEX: 
Effectively combine the strengths of SOM and DUPLEX  
### MDUPLEX: 
Improved and enhanced version of the DUPLEX method 


---
`How to use`: You can specify only the most basic parameter: inputFileName; Whether "I" is included in the partition calculation; splitting method  
example:  
`Win:DataSplittingAlgorithm.exe data.txt false SOMPLEX`  
`Linux: ./DataSplittingAlgorithm data.txt false SOMPLEX`  
Or you can specify more parameters, as shown below:  
`seed` = argv[4]: Provide seeds for some random methods, default = 1000  
`trainFraction` = argv[5]: Set the fraction of the training set, default = 0.6  
`testFraction` = argv[6]: Set the fraction of the test set, default = 0.2  
`outputTrain` = argv[7]: Set the file name of the output training set, default = train.txt  
`outputTest` = argv[8]: Set the file name of the output test set, default = test.txt  
`outputValid` = argv[9]: Set the file name of the output validation set, default = valid.txt  
example:  
`Win:DataSplittingAlgorithm.exe data.txt false SOMPLEX 1001 0.5 0.3 Tr.txt Ts.txt Vd.txt`  
`Linux: ./DataSplittingAlgorithm data.txt false SOMPLEX 1001 0.5 0.3 Tr.txt Ts.txt Vd.txt`  
