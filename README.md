# nftide-atacseq
Nextflow pipeline for bulk ATAC-seq

## Introduction ##
The pipeline uses cutadapt to remove nextera adaptors, aligns fastqs with bowtie2, remove PCR duplicates with picard, call peaks with macs2, generate fragment bed files with bedtools and genome track files (bw) with bedGraphToBigWig.

## Software dependencies ##
Dependencies  | Version
------------- | -------------
nextflow | 25.10
openjdk | 17.0.8
Perl | 5.32.1
python | 3.10.12
cutadapt | 4.5
samtools | 1.18
bowtie2 | 2.5.1
bedtools | 2.31.0
picard | 3.4.0
macs2 | 2.2.9.1
bedGraphToBigWig | 2.10
matplotlib | 3.10.8
pandas | 3.0.2
numpy | 2.4.3

If genome coverage files (.bw) are desired, please install bedGraphToBigWig manually (or with conda), and get the chromsize file from UCSC (download or using fetchChromSizes).

## Installation ##
(1) Create a conda environment with 
```
mamba create -n nftide-atacseq nextflow=25.10.2 python cutadapt samtools bowtie2 bedtools picard macs2 pandas numpy matplotlib ucsc-bedgraphtobigwig
``` 
and activate the environment with  
```
mamba activate nftide-atacseq
```

(2) Clone the repository with `git clone`, and execute
```
cd nftide-atacseq
```

(3) You also need to prepare bowtie2 index files for your genome. Please refer to their manuals to generate indexed genome files.

## Usage ##
(1) Prepare the `samplesheet.csv`. The csv file __must__ contain 3 columns with defined column names:  
`sample`: Name of the sequenced library. For example, `demo-1`. It will be the prefix of the output. Note: Different fastqs with same sample name will be merged before processing.  
`fastq_1`: Path to read 1.  
`fastq_2`: Path to read 2.  

(2) Run nextflow pipeline.
```
nextflow run atacseq_pe.nf \
  -output-dir outdir \
  --genome hg38 \
  --genomeDir bowtie2_index \
  --input_csv samplesheet.csv \
  --run_macs2 true \
  --gsize hs \
  --make_bw true \
  --chromsize hg38.chrom.sizes \
  -with-report outdir/nf_atac_report.html \
  -with-timeline outdir/nf_atac_timeline.html
```
where:
__--genome__ gives the prefix of the outputs.  
__--samplesheet.csv__ must have 3 columns named "sample", "fastq_1" and "fastq_2". Refer to samplesheet.csv in the repository.  
__--genomeDir__ is the bowtie2 index.  
__-output-dir__ is the output directory.  
If __--run_macs2__ is set to true, __--gsize__ must be set as input for the __-g__ parameter of `macs2 callpeak`.  
If __--make_bw__ is set to true, __--chromsize__ must be set as input for `bedtools genomecov` and `bedGraphToBigWig`. Refer to hg38.chrom.sizes in the repository.  

`-output-dir`: Path to the output directory.  
`--input_csv`: Path to samplesheet.csv as described in **step (1)**.  
`--genome`: the genome prefix of the outputs.  
`--genomeDir`: path to bowtie2 index.  
If `--run_macs2` is set to true, `--gsize` must be set as input for the `-g` parameter of `macs2 callpeak`.  
If `--make_bw` is set to true, `--chromsize` must be set as input for `bedtools genomecov` and `bedGraphToBigWig`. Refer to hg38.chrom.sizes in the repository.  

By default, the pipeline allows 2 samples to be processed in parallel. To change this behavior, modify _maxForks_ in __nextflow.config__.

## Expected output ##
The pipeline creates subfolders (named by samples in the samplesheet) in -output-dir.  
__Note__: the start and end coordinates of *_fragments.tsv files have been shifted by +4/-5 bp.



