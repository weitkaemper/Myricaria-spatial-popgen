# Spatial population genetics pipeline for target capture sequence data 


## Preparing `CAPTUS` output for analyses: `bash` shell, Box A (Figure 2)

### Make BAMS (Binary Alignment Maps)

Use the program `Strobealign` to match Illumina Sequencing reads to a reference genome. In our case, this is specimen 14, sequenced from fresh plant material.

Input material is `CAPTUS` output in the form `/01_clean_reads` with default folder settings
From the folder `/01_clean_reads`:

```sh
../strobealign ../ref_clean.fa{1}_R1.fastq.gz {1}_R2.fastq.gz
samtools sort -o {1}.bam
```

### Perform variant calling

To perform variant calling:

```sh
bcftools mpileup -f ../ref_clean.fa -b bam_list.txt -a FORMAT /DP, AD -q20 -Q20 -Ou /bcftools call -mv -Oz -o
```
Name of output file is `all.raw.vcf.gz`. To extract biallelic SNPs, input this to:

```sh
bcftools view -m2 -M2 -v snps -Oz -o all.biallelic.vcf.gz all.raw.vcf.gz
```
Now, `all.biallelic.vcf.gz` is the output file. 

### Check data quality

You can now check quality of data with:

```sh
plotvcfstats
```

If there are too many transversions, apply the following quality-control step on file `all.biallelic.vcf.gz`:

```sh
bcftools view -i QUAL >=500 -Oz -o all.bifiltered.vcf.gz all.biallelic.vcf.gz
```
The output file is `all.bifiltered.vcf.gz`. 

### Prune data

In order to prune for linkage disequilibrium on `all.bifiltered.vcf.gz`:

```sh
bcftools +prune -m 0.2 -w 50000 all.pruned.vcf.gz all.bifiltered.vcf.gz
```
The output file is `all.pruned.vcf.gz`. 

Now we can remove sites with too much missing data, where less than 90% of individuals are represented:

```sh
vcftools --gzvcf all.pruned.vcf.gz --recode -INFO -all --max-missing 0.9 --out all.pruned.miss10 --recode
```

Output file is `all.pruned.miss10`. 

Conversely, we now remove individuals with too much missing data, i.e. individuals with more than 15% sites missing. 

```sh
vcftools --gzvcf all.pruned.miss10 --missing -indiv --out missing-indiv 
```

From this information, manually create a list of such individuals and list them in a file called `missing_indiv.imiss`. Using this file, run:

```sh
bcftools view -S ^missing_indiv.imiss --force -samples -Oz -o all.pruned.nomiss_samples.vcf.gz all.pruned.miss10.recode.vcf
```
Output file is `all.pruned.nomiss_samples.vcf.gz`. 

### Convert to genotype matrix

Convert this into genotype matrix with:

```sh
vcftools --gzvcf all.pruned.nomiss_sample.vcf.gz --012 --out output_geno
```
This outputs three files, namely ".012", ".0l2.indv", and ".012.pos". 

## Preparing genotype and coordinate matrices: `R`, Box B (Figure 2)



### Convert `output_geno` into `allele_frequency_matrix.tsv` 

Now format the three output_geno files into one matrix in 0,0.5,1 style, reordered so that specimen numbers are numerically ascending:

```r
output_geno <- as.matrix(read.table("output_geno.012",header=FALSE,sep="\t"))
output_geno <- output_geno[,2:2184]
output_geno[output_geno == -1] <- NA
output_geno <- output_geno/2
indv <- read.table('output_geno.012.indv')[,1]
indv <- stringr::str_split_i(indv,stringr::fixed("."),1)
rownames(output_geno) <- indv 
pos <- read.table('output_geno.012.pos')
pos <- paste(pos[,1], pos[,2], sep="_")
colnames(output_geno) <- pos
ordered <- order(as.integer(rownames(output_geno)))
allele_frequency_matrix <- output_geno[ordered,] 

```
### Load and format coordinate data

Load your coordinate data into `R` in the correct format, as a two column matrix, with X coordinates as the first column and Y coordinates as the second column. Row names and row order should be the same as the allele frequency matrix. (See `conStruct` vignette for formatting data.)

### Create distance matrix

Load package geodist:

```r
install.packages('geodist')
library.packages('geodist')
```

Using geodist, make a distances file of the location of each individual compared with all other individuals, which includes the row and column names. 

```r
Dist <- geodist(Coords,measure="geodesic")
rownames(Dist) <- rownames(Coords)
colnames(Dist) <- rownames(Coords)
```

### Subset datasets for analysis

If necessary, take only coordinate data in your dataset which also has genetic data. Do this step if some of your sequences are missing, but you still included the geographical data in your matrix.

```r
Coords_sub <- Coords[rownames(Coords) %in% rownames(allele_frequency_matrix),]
```

## `conStruct` analyses: `R`, Box C (Figure 2)


Install and load `conStruct`:

```r
install.packages('conStruct')
library('conStruct')
```

Now we should focus on a geographical area to test population structure. In this example, we focus on Siberia, here defined by individuals in our dataset whose X coordinate > 76 and Y coordinate > 43.

### `ConStruct` analysis for geographically demarkated specimens

```r
Coords_siberia <- Coords_sub[Coords_sub[,"X"]>76,]
Coords_siberia <- Coords_siberia[Coords_siberia[,"Y"]>43,]
```
The allele frequency needs to contain the same individuals as the coordinate dataset for Siberia.

```r
freqs_siberia <- allele_frequency_matrix[rownames(allele_frequency_matrix) %in% rownames(Coords_siberia),]
```

The distance matrix also needs to contain the same individuals as the coordinate dataset for Siberia.
```r
Dist_siberia <- Dist[rownames(Dist) %in% rownames(Coords_siberia),colnames(Dist) %in% rownames(Coords_siberia)]
```
Run spatial model:

```r
run_siberia_SPATIAL <- conStruct(spatial = TRUE, 
                    K = 2, 
                    freqs = freqs_siberia,
                    geoDist = Dist_siberia, 
                    coords = Coords_siberia,
                    prefix = "siberia_SPATIAL",
		    n.chains = 20,
		    n.iter = 2000)
```

Check for chain with largest maximal log-posterior density by comparing 
```r
run_siberia_SPATIAL$chain_N$MAP$lpd
```
where N varies from 1 to 20. In this case, the maximum is attained by chain 2. 

Run non-spatial model:

```r
run_siberia_NONSPATIAL <- conStruct(spatial = FALSE, 
                    K = 2, 
                    freqs = freqs_siberia,
                    geoDist = Dist_siberia, 
                    coords = Coords_siberia,
                    prefix = "siberia_NONSPATIAL",
		    n.chains = 20,
		    n.iter = 2000)
```

Check for chain with largest maximal log-posterior density by comparing 
```r
run_siberia_NONSPATIAL$chain_N$MAP$lpd
```
where N varies from 1 to 20. In this case, the maximum is attained by chain 9. 

### Subset data for `EEMS` analysis

Subset the data according to the `conStruct` predominant layer in a `conStruct` run: 

First, create table of admixture proportions

```r
admix_proportion_siberia <- run_siberia_SPATIAL$chain_2$MAP$admix.proportions
rownames(admix_proportion_siberia) <- rownames(Coords_siberia)
```
Next, subset the data by the gene proportions in each layer.

```r
siberia_layer1 <- admix_proportion_siberia[admix_proportion_siberia[,1]>0.5,]
siberia_layer2 <- admix_proportion_siberia[admix_proportion_siberia[,2]>0.5,]
```
Now subset the other input objects to `EEMS`, such as the coordinates.

```r
freqs_siberia_layer1 <- allele_frequency_matrix[rownames(allele_frequency_matrix) %in% rownames(siberia_layer1),]
freqs_siberia_layer2 <- allele_frequency_matrix[rownames(allele_frequency_matrix) %in% rownames(siberia_layer2),]

Dist_siberia_layer1 <- Dist[rownames(Dist) %in% rownames(siberia_layer1),colnames(Dist) %in% rownames(siberia_layer1)]
Dist_siberia_layer2 <- Dist[rownames(Dist) %in% rownames(siberia_layer2),colnames(Dist) %in% rownames(siberia_layer2)]

Coords_siberia_layer1 <- Coords[rownames(Coords) %in% rownames(siberia_layer1),]
Coords_siberia_layer2 <- Coords[rownames(Coords) %in% rownames(siberia_layer2),]

```
### `conStruct` for specimens defined by species

We read a `csv` file containing all specimens with their putative species label and subset the conStruct input files for those of the putative species, in this case `M. germanica`:

```r
species_df <- read.table('SpecimenID_speciesName.csv', sep=';')
germanicas <- species_df[species_df[,2] == "M.germanica", 1]
Coords_germanica <- Coords_sub[strtoi(rownames(Coords_sub)) %in% germanicas,]

freqs_germanica <- allele_frequency_matrix[rownames(allele_frequency_matrix) %in% rownames(Coords_germanica),]

Dist_germanica <- Dist[rownames(Dist) %in% rownames(Coords_germanica),colnames(Dist) %in% rownames(Coords_germanica)]

```
Now run the spatial `conStruct` analysis on those:
```r
run_germanica_SPATIAL <- conStruct(spatial = TRUE, 
                    K = 2, 
                    freqs = freqs_germanica,
                    geoDist = Dist_germanica, 
                    coords = Coords_germanica,
                    prefix = "germanica_SPATIAL",
		    n.chains = 20,
		    n.iter = 2000)
```
Check for chain with largest maximal log-posterior density by comparing 
```r
run_germanica_SPATIAL$chain_N$MAP$lpd
```
where N varies from 1 to 20. In this case, the maximum is attained by chain 13. 

We now run the non-spatial `conStruct` analysis:
```r
run_germanica_NONSPATIAL <- conStruct(spatial = FALSE, 
                    K = 2, 
                    freqs = freqs_germanica,
                    geoDist = Dist_germanica, 
                    coords = Coords_germanica,
                    prefix = "germanica_NONSPATIAL",
		    n.chains = 20,
		    n.iter = 2000)
```
Check for chain with largest maximal log-posterior density by comparing 
```r
run_germanica_NONSPATIAL$chain_N$MAP$lpd
```
where N varies from 1 to 20. In this case, the maximum is attained by chain 4. 

Format conStruct plots so that the colours are consistent for presentation.
```r
load("germanica_NONSPATIAL_LONG_data.block.Robj")
nonspatial_block <- data.block

make.plots.matched(
	conStruct.results = run_germanica_NONSPATIAL_LONG, 
	prefix = "germanica_nonspatial_matched", 
	data.block = nonspatial_block, 
	sample.names = rownames(nonspatial_block$coords))

load("germanica_SPATIAL_LONG_data.block.Robj")
spatial_block <- data.block

make.plots.matched(
	conStruct.results = run_germanica_SPATIAL_LONG, 
	prefix = "germanica_spatial_matched", 
	data.block = spatial_block, 
	sample.names = rownames(spatial_block$coords))

load("siberia_NONSPATIAL_LONG_data.block.Robj")
siberia_nonspatial_block <- data.block

make.plots.matched(
	conStruct.results = run_siberia_NONSPATIAL_LONG, 
	prefix = "siberia_nonspatial_matched", 
	data.block = siberia_nonspatial_block, 
	sample.names = rownames(siberia_nonspatial_block$coords))

load("siberia_SPATIAL_LONG_data.block.Robj")
siberia_spatial_block <- data.block

make.plots.matched(
	conStruct.results = run_siberia_SPATIAL_LONG, 
	prefix = "siberia_spatial_matched", 
	data.block = siberia_spatial_block, 
	sample.names = rownames(siberia_spatial_block$coords))
```
## `EEMS` analyses: `R`, Box D (Figure 2)

Install and load `reems` alongside the recommended packages `rworldmap` and `rworldxtra` used to overlay our `EEMS` output on a political map of the study area. 

```r
install.packages(c('rworldmap', 'rworldxtra','reems'))
library('reems')
```

We create a custom maps for the area outline for `EEMS`
```r
layer1_outline <- matrix(c(110.5, 88.8, 88.8, 110.5, 110.5, 54.6, 52, 49.9, 49.9, 54.6), ncol = 2)
layer2_outline <- matrix(c(76, 77, 80, 83, 90, 90, 76, 49, 42.5, 43.5, 46.5, 50, 52.5, 49), ncol = 2)
```
We run an `EEMS` analysis for our Siberian data as follows. 

```r
eems_siberia_layer1 <- eems(
  freqs = 2*freqs_siberia_layer1,
  coords = Coords_siberia_layer1,
  mcmcpath = file.path(getwd(),"eems_siberia_layer1"),
  nChains = 3,
  outer = layer1_outline,
  nDemes = 400,
  numMCMCIter = 2000000,
  numBurnIter = 1000000,
  numThinIter = 9999,
  parallel = FALSE
)
 
eems_siberia_layer2 <- eems(
  freqs = 2*freqs_siberia_layer2,
  coords = Coords_siberia_layer2,
  mcmcpath = file.path(getwd(),"eems_siberia_layer2"),
  nChains = 3,
  outer = layer2_outline,
  nDemes = 400,
  numMCMCIter = 2000000,
  numBurnIter = 1000000,
  numThinIter = 9999,
  parallel = FALSE
)
```
Finally, we  plot the `EEMS` output.
```r
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(
  mcmcpath = eems_siberia_layer1,
  plotpath = "eems_siberia_layer1_grid_deme",
  longlat = TRUE,
  out.png = FALSE,
  xpd = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  col.demes = "black",
  projection.in = projection_none,
  projection.out = projection_mercator,
  add.outline = TRUE,
  col.outline = "gray90",
  add.map = TRUE,
)


eems.plots(
  mcmcpath = eems_siberia_layer2,
  plotpath = "eems_siberia_layer2_grid_deme",
  longlat = TRUE,
  out.png = FALSE,
  xpd = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  col.demes = "black",
  projection.in = projection_none,
  projection.out = projection_mercator,
  add.outline = TRUE,
  col.outline = "gray90",
  add.map = TRUE,
)

```
