# RNAsnv
## Section 1 - Overview

RNAsnv is a pipeline for calling genomic SNVs and RNA edit sites from RNA-seq data alone while
maintaining both high sensitivity and specificity. 

RNAsnv is called using the following format:

```
/path/to/RNAsnv.pl -c <configuration file> [options] [BAM1 BAM2 ...]

    Options:
            -help           this help message
            -man            displays RNAsnv man page
            -c              configuration file
            -train          train model using parameters in config file
            -predict        predict variant site origin using parameters in config file
            -wgs            WGS BAM file file file file
```

Note that RNAsnv expects a provided WGS bam file to be sorted in coordinate format. This can be done using picard tools, e.g.

```
java -Xmx10g -Djava.io.tmpdir=tmp -jar SortSam.jar \ 
	I=WGS.bam O=WGS.sorted.bam SO=coordinate 
```

Note: The tmpdir option establishes a local directory for temporary files used during sorting, which prevents small partitions
such as /tmp from filling too quickly.

## Section 1.1 - Compiling RNAsnv

Before compiling RNAsnv, make sure that bamtools is installed (see Section 2) and modify Makefile to point at bamtools /lib and /include directories
Also make sure to set LD_LIBRARY_PATH to include bamtools /lib (e.g. export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/src/bamtools/lib) when running RNAsnv, 
otherwise you'll get an error message with something like: 

```
"./remove_3prime_5prime_variants: error while loading shared libraries: libbamtools.so.2.3.0: 
cannot open shared object file: No such file or directory"
```

Additionally, a configuration file template is provided with the RNAsnv distribution. A trained model for
predicting variant site origin can be downloaded from X.



## Section 1.2 - Citing RNAsnv

When citing RNAsnv, please use the following reference:





#####################################################################################################################
#
## Section 2 - RNAsnv dependencies
 
RNAsnv requires the following programs to be installed correctly. See section X for step-by-step examples of how 
to install the programs on a linux system.

Program | Tested version | Download link
------------- | -------------| -------------
Bamtools | 2.3.0 		| 		https://github.com/pezmaster31/bamtools
samtools |  1.1 and 1.2	| 		http://www.htslib.org/download/
bcftools |  1.1 and 1.2	| 		http://www.htslib.org/download/
picard |  1.08			| 		http://broadinstitute.github.io/picard/
bedtools |  2.23		| 			https://github.com/arq5x/bedtools2

## Section 2.1 - RNAsnv dependencies step-by-step installation

# samtools
```
git clone https://github.com/samtools/samtools.git
cd samtools
make
make prefix=/usr/users/bdownie install
ls -l samtools
```

# bcftools
```
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
make prefix=/usr/users/bdownie install
ls -l bcftools
```

# Bamtools (requires cmake)
```
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools/
mkdir build
cd build/
cmake ..
make
ls -l bin/bamtools
```

# Bedtools 
```
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
make
ls -l bin/bedtools
```

#####################################################################################################################

#####################################################################################################################
#
## Section 3 - Generating RNAsnv BED files

Generating bed files appropriate for use in RNAsnv, using hg19 as an example
Note that UCSC stores human chromosomes names as "chr1", chr2, etc, while
Ensembl chromosome names are stored by chromosome number. For all programs to work correctly (RNAsnv, bedtools, etc)
chromosome names must be consistent across all files. For our purposes, we remove "chr" from all annotation files.

#
## Section 3.1 - Repeat masker track. Change as appropriate for genome of interest.

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome
gunzip -c rmsk.txt.gz  | awk '{print $6 "\t" $7 "\t" $8 }' | bedtools complement -i - -g hg19.genome | awk '{print $1 "\t" $2 "\t" $3 "\tNotRepeat\tNotRepeat"}'   > hg19.not_repeats.bed &
gunzip -c rmsk.txt.gz | awk '{print $6 "\t" $7 "\t" $8 "\t" $13 "\t" $12}' |   sed 's/\?//g' | cat - hg19.not_repeats.bed | bedtools sort | sed 's/^chr//' >  hg19.rmsk.txt.bed &
~/src/RNAsnv/parse_repeats.pl hg19.rmsk.txt.bed > hg19.repeat_intervals.bed
```

All intermediary files can be deleted after hg19.repeat_intervals.bed is generated


## Section 3.2 - Intron track  from UCSC

Download tracks from UCSC table browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
Step-by-step instructions available from: https://www.biostars.org/p/13290/

```
gunzip -c ucsc.hg19.introns.bed.gz | sed 's/^chr//' > ucsc.hg19.introns.bed
```

## Section 3.3 - RNA edit sites from RADAR

Download edit sites from RADAR (http://rnaedit.com/download/)
Note that sites are only available for hg19, mm9, and dm3. If later versions of these genomes are used for alignment,
use UCSC liftover tool (https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert coordinates to the appropriate reference release version

```
wget http://www.stanford.edu/~gokulr/database/Human_AG_all_hg19_v2.txt
tail -n +2 Human_AG_all_hg19_v2.txt   | sed 's/^chr//' | awk '{print $1 "\t" $2 "\t" $2 "\t" $3}' > Human_AG_all_hg19_v2.bed
```
