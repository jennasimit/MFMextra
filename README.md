# MTFMextra

This R package provides the simulation functions used to assess the joint fine-mapping methods of the MTFM R package.

## Simulation of case-control data for one disease

Below is an example of how to use hapgen2 (http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) to simulate N0 cases and N1 controls in a gene (IL2RA) on chromosome 10, where we set two causal variants.
The positions of the causal variants are $c1 and $c2, which have $ORhet1 and $ORhet2 as the respective odds ratios relating the odds of disease in heterozygote carriers of the 
non-reference allele (this is specified by the '1' following each of $c1 and $c2; '0' indicates with respect to the reference allele) compared to the homozygote reference allele. 
We assume a multiplicative model so that OR for homozygote disease risk is OR^2. For each simulation the causal variants were randomly selected from a certain SNP group.

This requires map, legend and hap files, which could be obtained from a reference panel (e.g. CEU of 1000 Genomes: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) vcf file by 
using vcftools (https://vcftools.github.io/index.html) with the --IMPUTE option. 
To output only a subset of SNPs, we specify the SNP positions, line-by-line, in the file keep-snps.txt.

ORhom1=$(bc -l <<< "$ORhet1 ^2")
ORhom2=$(bc -l <<< "$ORhet2 ^2")

./hapgen2 \
-m ./genetic_map_chr10_combined_b37.txt #map file \
-l ./IL2RA.impute.legend # legend file \
-h ./IL2RA.impute.hap # hap file \
-n $N0 $N1 # N0 controls, N1 cases \
-dl $c1 1 $ORhet1 $ORhom1 $c2 1 $ORhet2 $ORhom2 \
-no_haps_output -no_gens_output \
-t ./keep-snps.txt  # output for subset of SNPs that are listed in this file \
-Ne 11418 # effective population size; 11418 recommended for CEU \
-o ./CCdata  # prefix of output files




## Simulation of null data that forms the basis for simulating a specific region for two diseases with shared controls

Below is an example of how to use hapgen2 (http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) to simulate 100,000 individuals with no associated SNPs in a gene (IL2RA) on chromosome 10. 
This requires map, legend and hap files, which could be obtained from a reference panel (e.g. CEU of 1000 Genomes: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) vcf file by using vcftools (https://vcftools.github.io/index.html) with the --IMPUTE option.
By not specifying the causal variants, the default is to simulate null data.
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

