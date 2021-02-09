# Pyclone pipeline
Pipeline to run Pyclone, from sequenza, mutect and mpileup results.

1. given a fixed number *n* of mutations to input to PyClone, and that we have total number of *t* mutations, in which *m* are coding:
  * if *t<=n*, use all *t* mutations
  * if *m>n*, then randomly select *n* coding mutations  
  * if *m==n*, then use all *n* coding mutations  
  * if *m<n*, then randomly select *(n-m)* from *(t-m)* non-coding mutations to top up  
2. can choose if you want to use mutations in copy number neutral only by using flag variable *CN_neutral* (default *TRUE*) and numeric variable *CN* (*CN* is the number of major and minor allele copy number in copy-neutral regions, *CN* should be in .csv masterfile input to this pipeline)
  * *CN_neutral==FALSE*, randomly select mutations as indicated in **1**
  * if *CN_neutral==TRUE*, only mutations in intersected copy-neutral regions (i.e. regions with *A==B==CN* in sequenza) of all samples of the patient will be used. If we have *k* mutations in copy-neutral regions, treat *k* as *t* in **1**: 
    + If *k==n*, all and only these *k* mutations will be used  
    + If *k>n*, *n* mutations are to be randomly selected from *k*, as in **1**.  
    + If there are not enough mutations in the copy-neutral regions (i.e. *100<=k<n*), then run PyClone twice:
      + with only these *k* mutations in the copy-neutral regions, 
      + with randomly chosen mutation (i.e. as if *CN_neutral==FALSE*)  
    + If there are too few mutations in the intersected copy-neutral regions (*k<100*), then run as if *CN_neutral==FALSE*  
    Output number of mutations in intersected copy-neutral regions and whether a second run is conducted (i.e. if there are too few mutations in the copy-neutral regions)


Built on top of Hechuan's Pyclone pipeline.
