# Codes for "Extending Data Fission to Unsupervised Learning: Implications for Post-Clustering Inference"
An example of simulated data generated using the single-cell data simulator, scDesign3. We configure the simulation with a cell count of 1000 and vary the total number of genes at 1000, 2000, 3000, 4000, and 5000, while fixing the number of true differentially expressed (DE) genes at 40. We introduce gene-gene correlations among 100 non-DE genes. Differential expression analysis is performed using a naive method directly based on the `Seurat` pipeline, which leverages PCA for dimensionality reduction, the Louvain method for clustering via the FindClusters function, and the Wilcoxon test for differential expression testing via the FindMarkers function. Data fission is implemented based on the negative binomial distribution. Specifically, suppose `X \sim \text{NB}(r, \theta)`. We draw ` Y \sim \text{Bin}(X, p) ` with `p = 0.5`. Then, we define `f(X) = Y`  and `g(X) = X - Y`. Following data fission, differential expression analysis is likewise conducted using the  `Seurat` pipeline.
