# Data Splitting Methods Based on Data Distribution Consistency  
  
Several data partitioning algorithms:  
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
I provide a typical data file ("data.txt"), the program will read this file, and the data in the file can be divided into three subsets with consistent distribution characteristics (training set, test set, validation set). This is the rainfall and runoff data of a certain catchment for 10 years (precision is daily).  
  
The "data.txt" file gives a typical data format for this program. The first line of your file should contain "Idex","I",... ,"I","O"  
`"Idex"`: This column is a subscript column that simply identifies each row of data. It is convenient to use a continuous sequence, but this is not necessary  
  
`"I"..."I"`: This is a list of input vectors. In the field of hydrological models, there are rainfall (P), evapotranspiration (E), rainfall with a lag of N days (Pn), and so on. There can be multiple input vectors separated by a Tab '\t'(i.e., the Tab key)  
  
Note that you can choose not to include the input vectors("I"..."I") in the program operation, so that the program will only subset the output vector("O"); Otherwise, the program will compute the Euclidean distance between the rows to divide the subsets  
  
`"O"`: This is the output vector, which corresponds to the runoff Q in the hydrological domain  

`How to run`: (Default) You can only specify the most basic parameter: 1.inputFileName; 2.Whether "I"(i.e., the input vectors) is included in the partition calculation; 3.splitting method name  
example:  
`Windows:DataSplittingAlgorithm.exe data.txt false SOMPLEX`  
`Linux: ./DataSplittingAlgorithm data.txt false SOMPLEX`  
  
Or you can specify more parameters, as shown below:  
`seed` = argv[4]: Provide seeds for some random methods, default = 1000  
`trainFraction` = argv[5]: Set the fraction of the training set, default = 0.6  
`testFraction` = argv[6]: Set the fraction of the test set, default = 0.2  
`outputTrain` = argv[7]: Set the file name of the output training set, default = train.txt  
`outputTest` = argv[8]: Set the file name of the output test set, default = test.txt  
`outputValid` = argv[9]: Set the file name of the output validation set, default = valid.txt  
example:  
`Windows:  DataSplittingAlgorithm.exe data.txt false SOMPLEX 1001 0.5 0.3 Tr.txt Ts.txt Vd.txt`  
`Linux:  ./DataSplittingAlgorithm data.txt false SOMPLEX 1001 0.5 0.3 Tr.txt Ts.txt Vd.txt`  
  
My suggestion is to use my new two proposed methods to partition your data, i.e., SOMPLEX and MDUPLEX. Because the data distribution characteristics between the various subsets obtained by these two methods are more consistent, which is confirmed in a large number of hydrological historical data.
  
These codes are written in Visual Studio 2019

