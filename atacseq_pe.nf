#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow


// nextflow atacseq_pe.nf -with-report nf_atac_report.html -with-timeline nf_atac_timeline.html
// mandatory field: genomeDir, input_csv, gsize, -output-dir

params.genome = "hg38"
params.genomeBaseDir = "/home/zeemeeuw/data/ref"
params.genomeDir = "${params.genomeBaseDir}/${params.genome}/${params.genome}-bowtie2/${params.genome}"
params.input_csv = '/home/zeemeeuw/data/nextflow/data/atacseq/samplesheet.csv'
params.run_masc2 = true
params.gsize = 'hs'
params.make_bw = true
params.chromsize = '/home/zeemeeuw/data/nextflow/nf_pipe/nftide-atacseq/hg38.chrom.sizes'


process CUTADAPT {
    tag "cutadapt on ${id}"

    input:
    tuple val(id), path(read1), path(read2)

    output:
    tuple val(id), path("*_cutadapt_R1.fq.gz"), path("*_cutadapt_R2.fq.gz"), emit: trimmed_reads
    tuple val(id), path("*_cutadapt.log"), emit: cutadapt_log

    script:
    """
    cutadapt \
        -j ${task.cpus} -m 2:2 \
        -a "CTGTCTCTTATACACATCT" \
        -A "CTGTCTCTTATACACATCT" \
        --pair-filter=any \
        -o ${id}_cutadapt_R1.fq.gz \
        -p ${id}_cutadapt_R2.fq.gz \
        ${read1} \
        ${read2} > "${id}_cutadapt.log"

    """
}


process BOWTIE2 {
    tag "bowtie2 on ${id}"

    input:
    val genome
    tuple val(id), path(read1), path(read2)

    output:
    tuple val(id), path("*${genome}_aligned.sam"), emit: aligned_sam
    tuple val(id), path("*${genome}_bowtie2.log"), emit: bowtie2_log

    script:

    """
    echo "=== Starting bowtie2 alignment against ${genome} ===" > "${id}_${genome}_bowtie2.log"
    bowtie2 \
        -p ${task.cpus} -t -X2000 \
        --no-mixed --no-discordant \
        -x ${params.genomeDir} \
        -1 ${read1} \
        -2 ${read2} > ${id}_${genome}_aligned.sam 2>> "${id}_${genome}_bowtie2.log"
    """

}

process ADDREPLACERG {
    tag "AddreplaceRG to ${sam}"

    input:
    val genome
    tuple val(id), path(sam)

    output:
    tuple val(id), path("*${genome}_aligned.bam"), emit: rgtagged_bam
    
    script:

    """
    samtools addreplacerg -r "@RG\tID:RG1\tSM:${id}_${genome}" ${sam} -o ${id}_${genome}_aligned.bam
    """

}

process SORTBAM {
    tag "Sorting ${bam}"

    input:
    val genome
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*${genome}_sorted.bam"), emit: sorted_bam
    
    script:

    """
    samtools sort -@ ${task.cpus} ${bam} > ${id}_${genome}_sorted.bam
    """

}

process MARKDUPLICATE {
    tag "picard MarkDuplicates on ${bam}"

    input:
    val genome
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*${genome}_markdup.bam"), emit: markdup_bam
    
    script:

    """
    picard MarkDuplicates \
    --REMOVE_DUPLICATES false \
    -I ${bam} \
    -O ${id}_${genome}_markdup.bam \
    -M ${id}_${genome}_dups.txt
    """

}

process FILTERBAM {
    tag "Filtering PCR duplicates of ${bam}"

    input:
    val genome
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*${genome}_filtered.bam"), emit: filtered_bam
    
    script:

    """
    samtools view -b -f 2 -F 1804 -q 30 ${bam} -o ${id}_${genome}_filtered.bam
    """

}


process MACS2 {
    tag "Peak calling on ${bam} with macs2"

    input:
    val genome
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*${genome}_peaks*"), emit: called_peaks
    
    script:

    """
    macs2 callpeak \
        -f BAM \
        -g ${params.gsize} \
        -q 0.05 \
        --nomodel --shift 100 --extsize 200 --keep-dup all \
        -t ${bam} \
        -n ${id}_${genome}
    """

}

process MAKEFRAGMENT {
    tag "Making fragments from ${bam} with bedtools"

    input:
    val genome
    tuple val(id), path(bam)

    output:
    // tuple val(id), path("*${genome}_rmdup_sortbyname.bam"), emit: sortname_bam
    tuple val(id), path("*${genome}_rmchrM_fragments.tsv"), emit: filtfrags
    tuple val(id), path("*${genome}_fragments.tsv"), emit: nofiltfrags
    
    script:

    """
    samtools sort -n -@ ${task.cpus} ${bam} > ${id}_${genome}_rmdup_sortbyname.bam
    bedtools bamtobed -i ${id}_${genome}_rmdup_sortbyname.bam -bedpe | \
        perl -alne 'BEGIN{\$,="\t";} if(\$F[8] eq "+"){print(\$F[0],\$F[1]+4,\$F[5]-5);}elsif(\$F[8] eq "-"){print(\$F[0],\$F[4]+4,\$F[2]-5);}' | \
        sort -k1,1 -k2,2n > ${id}_${genome}_fragments.tsv
    grep -v -P "^chrM|^GL|^GH|^KI|^JH|^chrY|^MT" ${id}_${genome}_fragments.tsv > ${id}_${genome}_rmchrM_fragments.tsv
    """

}

process MAKEBW {
    tag "Making bw from ${bed}"

    input:
    val genome
    tuple val(id), path(bed)

    output:
    tuple val(id), path("*${genome}_ATAC.bw"), emit: bw

    script:
    """    
    bedtools genomecov -bg -i ${bed} -g ${params.chromsize} > ${id}_${genome}_ATAC.bedgraph
    bedGraphToBigWig ${id}_${genome}_ATAC.bedgraph ${params.chromsize} ${id}_${genome}_ATAC.bw
    """
}


workflow {
    main:
    def genome_basename
    if(!params.genome) {
        genome_basename = new File(params.genomeDir).getName()
    } else {
        genome_basename = params.genome
    }

    read_pairs_ch = channel.fromPath(params.input_csv)
    .splitCsv(header:true)
    .map { row -> [row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)] }

    CUTADAPT(read_pairs_ch)
    BOWTIE2(genome_basename, CUTADAPT.out.trimmed_reads)
    ADDREPLACERG(genome_basename, BOWTIE2.out.aligned_sam)
    SORTBAM(genome_basename, ADDREPLACERG.out.rgtagged_bam)
    MARKDUPLICATE(genome_basename, SORTBAM.out.sorted_bam)
    FILTERBAM(genome_basename, MARKDUPLICATE.out.markdup_bam)
    if(params.run_masc2){
        called_peaks = MACS2(genome_basename, FILTERBAM.out.filtered_bam)
    }else{
        called_peaks = channel.empty()
    }
    MAKEFRAGMENT(genome_basename, FILTERBAM.out.filtered_bam)

    if(params.run_masc2){
        ch_bw = MAKEBW(genome_basename, MAKEFRAGMENT.out.nofiltfrags)
    }else{
        ch_bw = channel.empty()
    }
    

    publish:
    cutadapt_fastqs = CUTADAPT.out.trimmed_reads
    cutadapt_logs = CUTADAPT.out.cutadapt_log
    bowtie2_logs = BOWTIE2.out.bowtie2_log
    markdup_bams = MARKDUPLICATE.out.markdup_bam
    filtdup_bams = FILTERBAM.out.filtered_bam
    out_peaks = called_peaks
    // sortname_bams = MAKEFRAGMENT.out.sortname_bam
    filtfrags = MAKEFRAGMENT.out.filtfrags
    nofiltfrags = MAKEFRAGMENT.out.nofiltfrags
    bws = ch_bw

}

output {
    cutadapt_fastqs {
        path { sample, _f1, _f2 -> "${sample}/cutadapt" }
    }
    cutadapt_logs {
        path { sample, _f1 -> "${sample}/cutadapt" }
    }
    bowtie2_logs {
        path { sample, _f1 -> "${sample}/bowtie2" }
    }
    markdup_bams {
        path { sample, _f1 -> "${sample}/bowtie2" }
    }
    filtdup_bams {
        path { sample, _f1 -> "${sample}/bowtie2" }
    }
    out_peaks {
        path { sample, _f1 -> "${sample}/macs2" }
    }
    // sortname_bams {
    //     path { sample, _f1 -> "${sample}/bowtie2" }
    // }
    filtfrags {
        path { sample, _f1 -> "${sample}/fragments" }
    }
    nofiltfrags {
        path { sample, _f1 -> "${sample}/fragments" }
    }
    bws {
        path { sample, _f1 -> "${sample}/fragments" }
    }
}