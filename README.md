Data splitting methods - source codes
Part of my graduate student project including several main data splitting methods:

Systematic: Systematic sampling is another deterministic approach in which every kth observation is sampled. If the data are ordered in some way, this implicitly generates a stratified sample, with stratification on the ordinal variable. One
approach is to sort the data along the output variable dimension to obtain a representative sample of the output variable distribution (Baxter, Stanley, Zhang, & Smith, 2000). This approach is easy to implement, as it assumes that the output variable can be mapped to a unique input state. However, this assumption may not hold in multivariate datasets where multiple input states might give rise to the same output, where the method cannot ensure that representative input–output combinations will be sampled, since only the output variable is considered.

DUPLEX: CADEX, or Kennard–Stone sampling (Kennard & Stone, 1969), is one of the earliest algorithms designed for data splitting. The approach iteratively draws samples based on distance, selecting points farthest away from those already included in the sample, and ensures maximum coverage of the data. An improved version called DUPLEX was proposed by Snee (1977); it is used widely in the field of chemometrics, including several ANN applications (Despagne & Massart, 1998; Sprevak, Azuaje, & Wang, 2004). However, the computational complexity of this algorithm may prohibit its use on
large datasets.

SBSS: SOM-based stratified sampling The SBSS approach involves two steps (Bowden et al., 2002). In the first step, the data are partitioned into K strata (clusters) using a selforganizing map (SOM; Kohonen, 1990), which considers the distances between data points. In the second step, data for the calibration and evaluation subsets are obtained by sampling from each of these strata. For SBSS-P, the sampling is done in proportion to the number of samples in each stratum. For SBSS-N, the sample allocation is increased for strata that contain a larger number of data points, or where the data points within a stratum have a larger variance (May et al., 2010).

SOMPLEX: Effectively combine the strengths of SOM and DUPLEX

MDUPLEX: Improved and enhanced version of the DUPLEX method
