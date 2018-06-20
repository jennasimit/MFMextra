# MTFMextra

This R package provides the simulation functions used to assess the joint fine-mapping methods of the MTFM R package.

## Simulation of null data that forms the basis for simulating a specific region for two diseases with shared controls

Below is an example of how to use hapgen2 (http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) to simulate 100,000 individuals with no associated SNPs in a gene (IL2RA) on chromosome 10. 
This requires map, legend and hap files, which could be obtained from a reference panel (e.g. CEU of 1000 Genomes: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) vcf file by using vcftools (https://vcftools.github.io/index.html) with the --IMPUTE option.
To output only a subset of SNPs, we specify the SNP positions, line-by-line, in the file keep-snps.txt.

./hapgen2 \
-m ./genetic_map_chr10_combined_b37.txt #map file \
-l ./IL2RA.impute.legend # legend file \
-h ./IL2RA.impute.hap # hap file \
-n 100000 0 # 100,000 controls, 0 cases \
-no_haps_output -no_gens_output \
-t ./keep-snps.txt  # output for subset of SNPs that are listed in this file \
-Ne 11418 # effective population size; 11418 recommended for CEU \
-o ./null_100k  # prefix of output files

