# Spatial population genetics pipeline for target capture sequence data

## Prerequisites

Make sure the required command-line tools and the `R` system are available in a recent version:

* [`CAPTUS`](https://github.com/edgardomortiz/captus)
* [`strobealign`](https://github.com/ksahlin/StrobeAlign)
* [`samtools`](https://github.com/samtools/samtools)
* [`bcftools`](https://github.com/samtools/bcftools)
* [`vcftools`](https://github.com/vcftools/vcftools)
* [`R`](https://www.r-project.org/)

## Cleaning data using `CAPTUS`

Clean raw Illumina reads using `CAPTUS` default settings:

```sh
captus\_assembly clean -r
```

## Preparing `CAPTUS` output for analyses: `bash` shell, Box A (Figure 2)

### Make BAMS (Binary Alignment Maps)

Use the program `Strobealign` to match Illumina Sequencing reads to a reference genome. In our case, this is specimen 14, sequenced from fresh plant material.

Input material is `CAPTUS` output in the form `/01\_clean\_reads` with default folder settings
From the folder `/01\_clean\_reads`:

```sh
../strobealign ../ref\_clean.fa{1}\_R1.fastq.gz {1}\_R2.fastq.gz
samtools sort -o {1}.bam
```

### Perform variant calling

To perform variant calling:

```sh
bcftools mpileup -f ../ref\_clean.fa -b bam\_list.txt -a FORMAT /DP, AD -q20 -Q20 -Ou /bcftools call -mv -Oz -o
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

From this information, manually create a list of such individuals and list them in a file called `missing\_indiv.imiss`. Using this file, run:

```sh
bcftools view -S ^missing\_indiv.imiss --force -samples -Oz -o all.pruned.nomiss\_samples.vcf.gz all.pruned.miss10.recode.vcf
```

Output file is `all.pruned.nomiss\_samples.vcf.gz`.

### Convert to genotype matrix

Convert this into genotype matrix with:

```sh
vcftools --gzvcf all.pruned.nomiss\_sample.vcf.gz --012 --out output\_geno
```

This outputs three files, namely ".012", ".0l2.indv", and ".012.pos".

## Preparing genotype and coordinate matrices: `R`, Box B (Figure 2)



### Convert `output\_geno` into `allele\_frequency\_matrix.tsv`

Now format the three output\_geno files into one matrix in 0,0.5,1 style, reordered so that specimen numbers are numerically ascending:

```r
output\_geno <- as.matrix(read.table("output\_geno.012",header=FALSE,sep="\\t"))
output\_geno <- output\_geno\[,2:2184]
output\_geno\[output\_geno == -1] <- NA
output\_geno <- output\_geno/2
indv <- read.table('output\_geno.012.indv')\[,1]
indv <- stringr::str\_split\_i(indv,stringr::fixed("."),1)
rownames(output\_geno) <- indv 
pos <- read.table('output\_geno.012.pos')
pos <- paste(pos\[,1], pos\[,2], sep="\_")
colnames(output\_geno) <- pos
ordered <- order(as.integer(rownames(output\_geno)))
allele\_frequency\_matrix <- output\_geno\[ordered,] 

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
Coords\_sub <- Coords\[rownames(Coords) %in% rownames(allele\_frequency\_matrix),]
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
Coords\_siberia <- Coords\_sub\[Coords\_sub\[,"X"]>76,]
Coords\_siberia <- Coords\_siberia\[Coords\_siberia\[,"Y"]>43,]
```

The allele frequency needs to contain the same individuals as the coordinate dataset for Siberia.

```r
freqs\_siberia <- allele\_frequency\_matrix\[rownames(allele\_frequency\_matrix) %in% rownames(Coords\_siberia),]
```

The distance matrix also needs to contain the same individuals as the coordinate dataset for Siberia.

```r
Dist\_siberia <- Dist\[rownames(Dist) %in% rownames(Coords\_siberia),colnames(Dist) %in% rownames(Coords\_siberia)]
```

Run spatial model:

```r
run\_siberia\_SPATIAL <- vector("list", 5)

for (K in 1:5){
run\_siberia\_SPATIAL[[K]] <- conStruct(spatial = TRUE, 
                    K = K, 
                    freqs = freqs\_siberia,
                    geoDist = Dist\_siberia, 
                    coords = Coords\_siberia,
                    prefix = "siberia\_SPATIAL",
		    n.chains = 10,
		    n.iter = 2000)
}
```

Check for chain with largest maximal log-posterior density by comparing

```r
run\_siberia\_SPATIAL[[K]]$chain\_N$MAP$lpd
```

where N varies from 1 to 10. 

Run non-spatial model:

```r
run\_siberia\_NONSPATIAL <- vector("list", 5)

for (K in 1:5){
run\_siberia\_NONSPATIAL[[K]] <- conStruct(spatial = FALSE, 
                    K = K, 
                    freqs = freqs\_siberia,
                    geoDist = Dist\_siberia, 
                    coords = Coords\_siberia,
                    prefix = "siberia\_NONSPATIAL",
		    n.chains = 10,
		    n.iter = 2000)
}
```

Check for chain with largest maximal log-posterior density by comparing

```r
run\_siberia\_NONSPATIAL[[K]]$chain\_N$MAP$lpd
```

where N varies from 1 to 10. 

### Subset data for `EEMS` analysis

Subset the data according to the `conStruct` predominant layer in a `conStruct` run:

First, create table of admixture proportions

```r
admix\_proportion\_siberia <- run\_siberia\_SPATIAL$chain\_2$MAP$admix.proportions
rownames(admix\_proportion\_siberia) <- rownames(Coords\_siberia)
```

Next, subset the data by the gene proportions in each layer.

```r
siberia\_layer1 <- admix\_proportion\_siberia\[admix\_proportion\_siberia\[,1]>0.5,]
siberia\_layer2 <- admix\_proportion\_siberia\[admix\_proportion\_siberia\[,2]>0.5,]
```

Now subset the other input objects to `EEMS`, such as the coordinates.

```r
freqs\_siberia\_layer1 <- allele\_frequency\_matrix\[rownames(allele\_frequency\_matrix) %in% rownames(siberia\_layer1),]
freqs\_siberia\_layer2 <- allele\_frequency\_matrix\[rownames(allele\_frequency\_matrix) %in% rownames(siberia\_layer2),]

Dist\_siberia\_layer1 <- Dist\[rownames(Dist) %in% rownames(siberia\_layer1),colnames(Dist) %in% rownames(siberia\_layer1)]
Dist\_siberia\_layer2 <- Dist\[rownames(Dist) %in% rownames(siberia\_layer2),colnames(Dist) %in% rownames(siberia\_layer2)]

Coords\_siberia\_layer1 <- Coords\[rownames(Coords) %in% rownames(siberia\_layer1),]
Coords\_siberia\_layer2 <- Coords\[rownames(Coords) %in% rownames(siberia\_layer2),]

```

### `conStruct` for specimens defined by species

We read a `csv` file containing all specimens with their putative species label and subset the conStruct input files for those of the putative species, in this case `M. germanica`:

```r
species\_df <- read.table('SpecimenID\_speciesName.csv', sep=';')
germanicas <- species\_df\[species\_df\[,2] == "M.germanica", 1]
Coords\_germanica <- Coords\_sub\[strtoi(rownames(Coords\_sub)) %in% germanicas,]

freqs\_germanica <- allele\_frequency\_matrix\[rownames(allele\_frequency\_matrix) %in% rownames(Coords\_germanica),]

Dist\_germanica <- Dist\[rownames(Dist) %in% rownames(Coords\_germanica),colnames(Dist) %in% rownames(Coords\_germanica)]

```

Now run the spatial `conStruct` analysis on those:

```r
run\_germanica\_SPATIAL <- vector("list", 5)

for (K in 1:5){
run\_germanica\_SPATIAL[[K]] <- conStruct(spatial = TRUE, 
                    K = K, 
                    freqs = freqs\_germanica,
                    geoDist = Dist\_germanica, 
                    coords = Coords\_germanica,
                    prefix = "germanica\_SPATIAL",
		    n.chains = 10,
		    n.iter = 2000)
}
```

Check for chain with largest maximal log-posterior density by comparing

```r
run\_germanica\_SPATIAL[[K]]$chain\_N$MAP$lpd
```

where N varies from 1 to 10.

We now run the non-spatial `conStruct` analysis:

```r
run\_germanica\_NONSPATIAL <- vector("list", 5)

for (K in 1:5){
run\_germanica\_NONSPATIAL[[K]] <- conStruct(spatial = FALSE, 
                    K = K, 
                    freqs = freqs\_germanica,
                    geoDist = Dist\_germanica, 
                    coords = Coords\_germanica,
                    prefix = "germanica\_NONSPATIAL",
		    n.chains = 10,
		    n.iter = 2000)
}
```
Check for chain with largest maximal log-posterior density by comparing

```r
run\_germanica\_NONSPATIAL[[K]]$chain\_N$MAP$lpd
```

where N varies from 1 to 10.

### Cross-validation for conStruct runs

Run cross validation for the M. germanica subset for all of Eurasia

ger.xvals <- x.validation(train.prop = 0.9,
n.reps = 8,
K = 1:5,
freqs = freqs\_germanica,
data.partitions = NULL,
geoDist = Dist\_germanica,
coords = Coords\_germanica,
prefix = "ger\_CV",
n.iter = 1e3,
make.figs = FALSE,
save.files = FALSE,
parallel = TRUE,
n.nodes = 4)

Run cross validation for the M. germanica and M. longifolia for the Mountains of Southern Siberia and Kazakhstan subset

siberia.xvals <- x.validation(train.prop = 0.9,
n.reps = 8,
K = 1:5,
freqs = freqs_siberia,
data.partitions = NULL,
geoDist = Dist_siberia,
coords = Coords_siberia,
prefix = "siberia_CV",
n.iter = 1e3,
make.figs = FALSE,
save.files = FALSE,
parallel = TRUE,
n.nodes = 4)

### Plotting conStruct results consistently

We define a function make.plots.matched. This function ensures consistency in coloring the layers to prevent label switching.

```r
make.plots.matched <- function(
    conStruct.results,
    prefix,
    data.block,
    sample.names = rownames(data.block$coords),
    layer.colors = NULL,
    pdf.width = 10,
    pdf.height = 5,
    pie.width = 6,
    pie.height = 6
) {

    ## ============================================================
    ## Checks
    ## ============================================================

    if (!is.list(conStruct.results)) {
        stop(
            "\n\"conStruct.results\" must be the list containing ",
            "the results from all chains\n\n"
        )
    }

    if (is.null(names(conStruct.results)) ||
        !any(grepl("chain", names(conStruct.results), ignore.case = TRUE))) {

        stop(
            "\nyou must specify conStruct results across all chains\n",
            "i.e. from \"conStruct.results\" rather than ",
            "\"conStruct.results[[1]]\"\n\n"
        )
    }

    if (is.null(data.block$K)) {
        stop("\n\"data.block\" does not contain \"K\"\n\n")
    }

    if (is.null(data.block$N)) {
        stop("\n\"data.block\" does not contain \"N\"\n\n")
    }

    K <- data.block$K
    N <- data.block$N

    n.chains <- length(conStruct.results)
    chain.names <- names(conStruct.results)


    ## ============================================================
    ## Sample names
    ## ============================================================

    if (is.null(sample.names)) {
        sample.names <- rownames(data.block$coords)
    }

    if (is.null(sample.names)) {
        stop(
            "\n\"sample.names\" was not supplied and ",
            "rownames(data.block$coords) are NULL\n\n"
        )
    }

    if (length(sample.names) != N) {
        stop(
            "\nlength(sample.names) must equal data.block$N\n",
            "length(sample.names) = ", length(sample.names), "\n",
            "data.block$N = ", N, "\n\n"
        )
    }


    ## ============================================================
    ## Layer colours
    ## ============================================================

    if (is.null(layer.colors)) {

        layer.colors <- c(
            "blue",
            "red",
            "goldenrod1",
            "forestgreen",
            "darkorchid1",
            "deepskyblue",
            "darkorange1",
            "seagreen2",
            "yellow1",
            "black"
        )

        if (K > length(layer.colors)) {
            stop(
                "\nyou have specified more layers than there are ",
                "default colors.\n",
                "you must specify your own \"layer.colors\"\n\n"
            )
        }

        layer.colors <- layer.colors[seq_len(K)]

    } else {

        if (length(layer.colors) != K) {
            stop(
                "\nyou must specify one color per layer\n\n"
            )
        }
    }


    ## ============================================================
    ## Extract MAP admixture proportions from each chain
    ## ============================================================

    admix <- lapply(
        seq_len(n.chains),
        function(i) {

            x <- conStruct.results[[i]]

            if (is.null(x$MAP$admix.proportions)) {
                stop(
                    "\nchain ", i,
                    " does not contain \"MAP$admix.proportions\"\n\n"
                )
            }

            a <- as.matrix(x$MAP$admix.proportions)

            if (nrow(a) != N) {
                stop(
                    "\nchain ", i,
                    " contains ", nrow(a),
                    " samples, but data.block$N = ", N, "\n\n"
                )
            }

            if (ncol(a) != K) {
                stop(
                    "\nchain ", i,
                    " contains ", ncol(a),
                    " layers, but data.block$K = ", K, "\n\n"
                )
            }

            a
        }
    )


    ## ============================================================
    ## Match layer labels across chains
    ##
    ## Chain 1 is the reference.
    ##
    ## match.layers.x.runs() returns the permutation needed to
    ## reorder the layers in chain i so that they correspond to
    ## the layers in chain 1.
    ## ============================================================

    layer.orders <- vector("list", n.chains)

    layer.orders[[1]] <- seq_len(K)

    if (n.chains > 1) {

        for (i in 2:n.chains) {

            layer.orders[[i]] <- match.layers.x.runs(
                admix[[1]],
                admix[[i]]
            )
        }
    }

    names(layer.orders) <- chain.names


    ## ============================================================
    ## Apply matching
    ## ============================================================

    matched.admix <- lapply(
        seq_len(n.chains),
        function(i) {

            a <- admix[[i]][,
                layer.orders[[i]],
                drop = FALSE
            ]

            colnames(a) <- paste0(
                "Layer_",
                seq_len(K)
            )

            rownames(a) <- sample.names

            a
        }
    )

    names(matched.admix) <- chain.names


    ## ============================================================
    ## STRUCTURE / bar plots
    ##
    ## One page per chain.
    ##
    ## sample.names are passed directly to make.structure.plot()
    ## so that the names occur below the bars.
    ## ============================================================

    grDevices::pdf(
        file = paste0(
            prefix,
            "_matched_structure_plots.pdf"
        ),
        width = pdf.width,
        height = pdf.height
    )

    for (i in seq_len(n.chains)) {

        make.structure.plot(
            matched.admix[[i]],
            mar = c(7, 4, 2, 2),
            sample.order = NULL,
            layer.order = NULL,
            sample.names = sample.names,
            sort.by = NULL,
            layer.colors = layer.colors
        )

        graphics::title(
            main = paste0(
                "Chain ",
                i
            ),
            line = 0.5
        )
    }

    grDevices::dev.off()


    ## ============================================================
    ## PIE MAPS
    ##
    ## The SAME matched.admix matrices are used here, so the layer
    ## colours correspond exactly to those in the STRUCTURE plots.
    ## ============================================================

    if (K > 1) {

        grDevices::pdf(
            file = paste0(
                prefix,
                "_matched_pie_maps.pdf"
            ),
            width = pie.width,
            height = pie.height
        )

        for (i in seq_len(n.chains)) {

            make.admix.pie.plot(
                matched.admix[[i]],
                data.block$coords,
                layer.colors,
                radii = 2.7,
                add = FALSE,
                x.lim = NULL,
                y.lim = NULL
            )

            graphics::title(
                main = paste0(
                    "Chain ",
                    i
                )
            )
        }

        grDevices::dev.off()
    }


    ## ============================================================
    ## Save matched results and permutations
    ## ============================================================

    save(
        matched.admix,
        layer.orders,
        file = paste0(
            prefix,
            "_matched.RData"
        )
    )


    ## ============================================================
    ## Return useful objects invisibly
    ## ============================================================

    invisible(
        list(
            matched.admix = matched.admix,
            layer.orders = layer.orders,
            sample.names = sample.names,
            layer.colors = layer.colors
        )
    )
}

```

Format conStruct plots so that the colours are consistent for presentation.

```r
load("germanica\_NONSPATIAL\_LONG\_data.block.Robj")
nonspatial\_block <- data.block

make.plots.matched(
	conStruct.results = run\_germanica\_NONSPATIAL\_LONG, 
	prefix = "germanica\_nonspatial\_matched", 
	data.block = nonspatial\_block, 
	sample.names = rownames(nonspatial\_block$coords))

load("germanica\_SPATIAL\_LONG\_data.block.Robj")
spatial\_block <- data.block

make.plots.matched(
	conStruct.results = run\_germanica\_SPATIAL\_LONG, 
	prefix = "germanica\_spatial\_matched", 
	data.block = spatial\_block, 
	sample.names = rownames(spatial\_block$coords))

load("siberia\_NONSPATIAL\_LONG\_data.block.Robj")
siberia\_nonspatial\_block <- data.block

make.plots.matched(
	conStruct.results = run\_siberia\_NONSPATIAL\_LONG, 
	prefix = "siberia\_nonspatial\_matched", 
	data.block = siberia\_nonspatial\_block, 
	sample.names = rownames(siberia\_nonspatial\_block$coords))

load("siberia\_SPATIAL\_LONG\_data.block.Robj")
siberia\_spatial\_block <- data.block

make.plots.matched(
	conStruct.results = run\_siberia\_SPATIAL\_LONG, 
	prefix = "siberia\_spatial\_matched", 
	data.block = siberia\_spatial\_block, 
	sample.names = rownames(siberia\_spatial\_block$coords))
```

## `EEMS` analyses: `R`, Box D (Figure 2)

Install and load `reems` alongside the recommended packages `rworldmap` and `rworldxtra` used to overlay our `EEMS` output on a political map of the study area.

```r
install.packages(c('rworldmap', 'rworldxtra','reems'))
library('reems')
```

We create a custom maps for the area outline for `EEMS`

```r
layer1\_outline <- matrix(c(110.5, 88.8, 88.8, 110.5, 110.5, 54.6, 52, 49.9, 49.9, 54.6), ncol = 2)
layer2\_outline <- matrix(c(76, 77, 80, 83, 90, 90, 76, 49, 42.5, 43.5, 46.5, 50, 52.5, 49), ncol = 2)
```

We run an `EEMS` analysis for our Siberian data as follows.

```r
eems\_siberia\_layer1 <- eems(
  freqs = 2\*freqs\_siberia\_layer1,
  coords = Coords\_siberia\_layer1,
  mcmcpath = file.path(getwd(),"eems\_siberia\_layer1"),
  nChains = 3,
  outer = layer1\_outline,
  nDemes = 400,
  numMCMCIter = 2000000,
  numBurnIter = 1000000,
  numThinIter = 9999,
  parallel = FALSE
)
 
eems\_siberia\_layer2 <- eems(
  freqs = 2\*freqs\_siberia\_layer2,
  coords = Coords\_siberia\_layer2,
  mcmcpath = file.path(getwd(),"eems\_siberia\_layer2"),
  nChains = 3,
  outer = layer2\_outline,
  nDemes = 400,
  numMCMCIter = 2000000,
  numBurnIter = 1000000,
  numThinIter = 9999,
  parallel = FALSE
)
```

Finally, we  plot the `EEMS` output.

```r
projection\_none <- "+proj=longlat +datum=WGS84"
projection\_mercator <- "+proj=merc +datum=WGS84"

eems.plots(
  mcmcpath = eems\_siberia\_layer1,
  plotpath = "eems\_siberia\_layer1\_grid\_deme",
  longlat = TRUE,
  out.png = FALSE,
  xpd = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  col.demes = "black",
  projection.in = projection\_none,
  projection.out = projection\_mercator,
  add.outline = TRUE,
  col.outline = "gray90",
  add.map = TRUE,
)


eems.plots(
  mcmcpath = eems\_siberia\_layer2,
  plotpath = "eems\_siberia\_layer2\_grid\_deme",
  longlat = TRUE,
  out.png = FALSE,
  xpd = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  col.demes = "black",
  projection.in = projection\_none,
  projection.out = projection\_mercator,
  add.outline = TRUE,
  col.outline = "gray90",
  add.map = TRUE,
)

```

