library(sesame)

## As sesame and sesameData are under active development, this documentation is
## specific to the following version of R, sesame, sesameData and ExperimentHub:
sesame_checkVersion()

## CRITICAL: After a new installation, one must cache the associated annotation data using the following command. This needs to be done only once per SeSAMe installation/update. Caching data to local guarantees proper data retrieval and saves internet traffic.
sesameDataCache()
## You can find the location of the cached annotation data on your local computer using:
tools::R_user_dir("ExperimentHub", which="cache")

#### The openSesame Pipeline ####

## The following use the idat_dir in system.file,
## Replace it with your actual IDAT path
idat_dir = system.file("extdata/", package = "sesameData")
betas = openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(2))

##  The above openSesame call is equivalent to:
betas = do.call(cbind, BiocParallel::bplapply(
    searchIDATprefixes(idat_dir), function(pfx) {
        getBetas(prepSesame(readIDATpair(pfx), "QCDPB"))
}, BPPARAM = BiocParallel::MulticoreParam(2)))

## or even more explicitly (if one needs to control argument passed
## to a specific preprocessing function)
betas = do.call(cbind, BiocParallel::bplapply(
    searchIDATprefixes(idat_dir), function(pfx) {
        getBetas(noob(pOOBAH(dyeBiasNL(inferInfiniumIChannel(qualityMask(
            readIDATpair(pfx)))))))
}, BPPARAM = BiocParallel::MulticoreParam(2)))

betas = openSesame(idat_dir, func = getBetas) # getBetas is the default
sdfs = openSesame(idat_dir, func = NULL) # return SigDF list
allele_freqs = openSesame(idat_dir, func = getAFs) # SNP allele frequencies
sdfs = openSesame(sdfs, prep = "Q", func = NULL)   # take and return SigDFs

pvals = openSesame(idat_dir, func = pOOBAH, return.pval=TRUE)

#### Data Preprocessing ####

sdf = sesameDataGet('EPIC.1.SigDF')
sdf_preped = openSesame(sdf, prep="DB", func=NULL)

#### Preprocessing Function Code ####

prepSesameList()

