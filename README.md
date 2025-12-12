# nftide-atacseq
Nextflow pipeline for bulk ATAC-seq

## Introduction ##
The pipeline uses cutadapt to remove nextera adaptors, aligns fastqs with bowtie2, remove PCR duplicates with picard, call peaks with macs2, generate fragment bed files with bedtools and genome track files (bw) with bedGraphToBigWig.

## Software dependencies ##
Dependencies  | Version
------------- | -------------
nextflow | 25.10
Perl | 5.32.1
python | 3.10.12
cutadapt | 4.5
samtools | 1.18
bowtie2 | 2.5.1
bedtools | 2.31.0
picard | 3.4.0
macs2 | 2.2.9.1
bedGraphToBigWig | 2.10

If genome coverage files (.bw) are desired, please install bedGraphToBigWig manually, and get the chromsize file from UCSC (download or using fetchChromSizes).

## Usage ##
(1) Create a conda environment with 
```
conda create -n nftide-atacseq nextflow=25.10.2 python cutadapt samtools bowtie2 bedtools picard macs2
``` 
and activate the environment with  
```
conda activate nftide-atacseq
```

(2) Clone the repository with `git clone`, and execute
```
cd nftide-atacseq
```

(3) Run nextflow pipeline.
```
nextflow atacseq_pe.nf \
  -output-dir outdir \
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
__--samplesheet.csv__ must have 3 columns named "sample", "fastq_1" and "fastq_2".  
__--genomeDir__ is the bowtie2 index.  
__-output-dir__ is the output directory.  
If __--run_macs2__ is set to true, __--gsize__ must be set as input for the __-g__ parameter of `macs2 callpeak`.  
If __--make_bw__ is set to true, __--chromsize__ must be set as input for `bedtools genomecov` and `bedGraphToBigWig`.

By default, the pipeline allows 2 samples to be processed in parallel. To change this behavior, modify _maxForks_ in __nextflow.config__.

## Expected output ##
The pipeline creates subfolders (named by samples in the samplesheet) in -output-dir.



