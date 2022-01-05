## Changes in version 1.9.3 (2021-12-16)
+ Add LINCS2 database with file name of `lincs2020.h5`

## Changes in version 1.9.2 (2021-12-06)
+ Add `dest_path` parameter to `getCmapCEL` function

## Changes in version 1.2.1 (2020-08-14)
+ Edited vignette containing 'DEGs and Cutoffs Definition' subsections to
document how to defining query DEGs and gene set reference database by setting
LFC scores and FDR cutoffs.

## Changes in version 1.1.0 (2019-10-23)
+ Initial version 

## Changes in version 0.99.10 (2019-10-22)
+ Stored data in ExperimentHub
+ Saved LINCS gctx file to hdf5 file in batches
+ Saved cmap databases (cmap, cmap_rank, cmap_expr) to hdf5 files
+ Loaded hdf5 file into R as SummarizedExperiment object
  - hdf5 file includes matrix, rownames and colnames

## Changes in version 0.99.0 (2019-04-02)
+ Needed 50Gb memory to load matrix in LINCS gctx file and save as HDF5 backed SE object
