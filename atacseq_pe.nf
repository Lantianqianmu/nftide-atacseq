#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow


// nextflow run atacseq_pe.nf -with-report nf_atac_report.html -with-timeline nf_atac_timeline.html -bg
// mandatory field: genomeDir, input_csv, gsize, -output-dir

// nextflow run atacseq_pe.nf -with-report nf_atac_report.html -with-timeline nf_atac_timeline.html -resume

params.genome = "hg38"
params.genomeBaseDir = "/data/xrz/ref"
params.genomeDir = "${params.genomeBaseDir}/${params.genome}/${params.genome}-bowtie2/${params.genome}"
params.input_csv = '/data/xrz/LXYdata/nextflow/samplesheet.csv'
params.run_masc2 = false
params.gsize = 'hs'
// params.gsize = 3182
params.make_bw = false 
params.chromsize = '/data/xrz/charseq/nftide-atacseq/hg38.chrom.sizes'
params.calc_tss = true 
params.tss_file = "${params.genomeBaseDir}/${params.genome}/${params.genome}_ArchR_TSS.tsv"

process MERGE_FQ {
    tag "Merging fastq files of ${meta.id}..."
    
    input:
    tuple val(meta), path(r1s) , path(r2s)
    
    output:
    tuple val(meta), path("*_merged_R1.fq.gz"), path("*_merged_R2.fq.gz"), emit: merged_fq
    
    script:
    """
    cat ${r1s.join(' ')} > ${meta.id}_merged_R1.fq.gz
    cat ${r2s.join(' ')} > ${meta.id}_merged_R2.fq.gz
    """
}

process CUTADAPT {
    tag "cutadapt on ${meta.id}"

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("*_cutadapt_R1.fq.gz"), path("*_cutadapt_R2.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("*_cutadapt.log"), emit: cutadapt_log

    script:
    """
    cutadapt \
        -j ${task.cpus} -m 2:2 \
        -a "CTGTCTCTTATACACATCT" \
        -A "CTGTCTCTTATACACATCT" \
        --pair-filter=any \
        --overlap 1 \
        -o ${meta.id}_cutadapt_R1.fq.gz \
        -p ${meta.id}_cutadapt_R2.fq.gz \
        ${read1} \
        ${read2} > "${meta.id}_cutadapt.log"

    """
}


process BOWTIE2 {
    tag "bowtie2 on ${meta.id}"

    input:
    val genome
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("*${genome}_aligned.sam"), emit: aligned_sam
    tuple val(meta), path("*${genome}_bowtie2.log"), emit: bowtie2_log

    script:

    """
    echo "=== Starting bowtie2 alignment against ${genome} ===" > "${meta.id}_${genome}_bowtie2.log"
    bowtie2 \
        -p ${task.cpus} -t -X2000 \
        --no-mixed --no-discordant \
        -x ${params.genomeDir} \
        -1 ${read1} \
        -2 ${read2} > ${meta.id}_${genome}_aligned.sam 2>> "${meta.id}_${genome}_bowtie2.log"
    """

}

process ADDREPLACERG {
    tag "AddreplaceRG to ${sam}"

    input:
    val genome
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*${genome}_aligned.bam"), emit: rgtagged_bam
    
    script:

    """
    samtools addreplacerg -r "@RG\tID:RG1\tSM:${meta.id}_${genome}" ${sam} -o ${meta.id}_${genome}_aligned.bam
    """

}

process SORTBAM {
    tag "Sorting ${bam}"

    input:
    val genome
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*${genome}_sorted.bam"), emit: sorted_bam
    
    script:

    """
    samtools sort -@ ${task.cpus} ${bam} > ${meta.id}_${genome}_sorted.bam
    """

}

process MARKDUPLICATE {
    tag "picard MarkDuplicates on ${bam}"

    input:
    val genome
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*${genome}_markdup.bam"), emit: markdup_bam
    
    script:

    """
    picard MarkDuplicates \
    --REMOVE_DUPLICATES false \
    -I ${bam} \
    -O ${meta.id}_${genome}_markdup.bam \
    -M ${meta.id}_${genome}_dups.txt
    """

}

process FILTERBAM {
    tag "Filtering PCR duplicates of ${bam}"

    input:
    val genome
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*${genome}_filtered.bam"), emit: filtered_bam
    
    script:

    """
    samtools view -b -f 2 -F 1804 -q 10 ${bam} -o ${meta.id}_${genome}_filtered.bam
    """

}


process MACS2 {
    tag "Peak calling on ${bam} with macs2"

    input:
    val genome
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*${genome}_peaks*"), emit: called_peaks
    
    script:

    """
    macs2 callpeak \
        -f BAM \
        -g ${params.gsize} \
        -q 0.05 \
        --nomodel --shift 100 --extsize 200 --keep-dup all \
        -t ${bam} \
        -n ${meta.id}_${genome}
    """

}

process MAKEFRAGMENT {
    tag "Making fragments from ${bam} with bedtools"
    errorStrategy 'ignore'

    input:
    val genome
    tuple val(meta), path(bam)

    output:
    // tuple val(meta), path("*${genome}_rmdup_sortbyname.bam"), emit: sortname_bam
    tuple val(meta), path("*${genome}_rmchrM_fragments.tsv"), emit: filtfrags
    tuple val(meta), path("*${genome}_fragments.tsv"), emit: nofiltfrags
    
    script:

    """
    samtools sort -n -@ ${task.cpus} ${bam} > ${meta.id}_${genome}_rmdup_sortbyname.bam
    bedtools bamtobed -i ${meta.id}_${genome}_rmdup_sortbyname.bam -bedpe | \
        perl -alne 'BEGIN{\$,="\t";} if(\$F[8] eq "+"){print(\$F[0],\$F[1]+4,\$F[5]-5);}elsif(\$F[8] eq "-"){print(\$F[0],\$F[4]+4,\$F[2]-5);}' | \
        sort -k1,1 -k2,2n > ${meta.id}_${genome}_fragments.tsv
    grep -v -P "^chrM|^GL|^GH|^KI|^JH|^chrY|^MT" ${meta.id}_${genome}_fragments.tsv > ${meta.id}_${genome}_rmchrM_fragments.tsv
    """

}


process PLOTSIZE {
    tag "Plotting fragment length distribution from ${tsv}"
    errorStrategy 'ignore'

    input:
    val genome
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*png"), emit: fragsize_plots
    tuple val(meta), path("${meta.id}_${genome}_length_counts.txt"), emit: length_counts

    script:

    """
    python3 << 'EOF'
    import pandas as pd
    import matplotlib.pyplot as plt
    
    df = pd.read_csv("${tsv}", sep='\\t', header=None, names=['chr', 'start', 'end'])
    df['length'] = df['end'] - df['start']
    
    length_counts = df['length'].value_counts().sort_index().reset_index()
    length_counts.columns = ['length', 'count'] 
    length_counts.to_csv("${meta.id}_${genome}_length_counts.txt", sep='\\t', index=False)

    plt.figure(figsize=(8, 6))
    plt.plot(length_counts['length'], length_counts['count'], color='#8D1616', linewidth=1)
    
    plt.xlabel('Fragment Length (bp)')
    plt.ylabel('Count')
    plt.title(f'Fragment Length Distribution: ${meta.id}')
    plt.grid(True, alpha=0.3, linestyle=':')
    
    max_len = min(1000, length_counts['length'].max())
    plt.xlim(0, max_len)
    
    plt.tight_layout()
    plt.savefig("${meta.id}_${genome}_fragment_size.png", dpi=300)
    plt.close()
    EOF
    """ 
}


process MAKEBW {
    tag "Making bw from ${bed}"

    input:
    val genome
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*${genome}_ATAC.bw"), emit: bw

    script:
    """    
    bedtools genomecov -bg -i ${bed} -g ${params.chromsize} > ${meta.id}_${genome}_ATAC.bedgraph
    bedGraphToBigWig ${meta.id}_${genome}_ATAC.bedgraph ${params.chromsize} ${meta.id}_${genome}_ATAC.bw
    """
}

process PLOTTSS {
    tag "Calculating TSS enrichment from ${tsv}"

    input:
    val genome
    tuple val(meta), path(tsv)
    path pyscript_tss_enrichment

    output:
    tuple val(meta), path("*_tss_enrichment.png"), path("*_tss_enrichment.txt"), emit: tss_enrichment

    script:
    """    
    python ${pyscript_tss_enrichment} -i ${tsv} -t ${params.tss_file} -o ${meta.id}_${genome} -p ${task.cpus}
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

    ch_read_pairs = channel.fromPath(params.input_csv)
    .splitCsv(header:true)
    .map { row -> 
        [
            row.sample,
            row
        ]
    }
    .groupTuple()
    .map { _sample, rows -> 
        rows.withIndex().collect { row, index ->
            row + [rep: index + 1]
        }
    }
    .flatMap { item -> item }
   .map { row -> 

        [
            [
                id: row.sample,
                rep: row.rep,

            ], 
            [
                file(row.fastq_1, checkIfExists: true), 
                file(row.fastq_2, checkIfExists: true)
            ]
        ]
    }
    .map{meta, files -> [meta.subMap(['id']), files]}
    .groupTuple()
    .map { meta, filePairs ->
        [ meta, filePairs.collect { pair -> pair[0] }, filePairs.collect { pair -> pair[1] }]
    }

    MERGE_FQ(ch_read_pairs)
    CUTADAPT(MERGE_FQ.out.merged_fq)
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
    PLOTSIZE(genome_basename, MAKEFRAGMENT.out.filtfrags)

    if(params.make_bw){
        ch_bw = MAKEBW(genome_basename, MAKEFRAGMENT.out.nofiltfrags)
    }else{
        ch_bw = channel.empty()
    }
    if(params.calc_tss){
        ch_tss = PLOTTSS(genome_basename, MAKEFRAGMENT.out.nofiltfrags, "${projectDir}/../utils/tss_enrichment.py")
    }else{
        ch_tss = channel.empty()
    }

    

    publish:
    merged_fastqs = MERGE_FQ.out.merged_fq
    cutadapt_fastqs = CUTADAPT.out.trimmed_reads
    cutadapt_logs = CUTADAPT.out.cutadapt_log
    bowtie2_logs = BOWTIE2.out.bowtie2_log
    markdup_bams = MARKDUPLICATE.out.markdup_bam
    filtdup_bams = FILTERBAM.out.filtered_bam
    out_peaks = called_peaks
    filtfrags = MAKEFRAGMENT.out.filtfrags
    nofiltfrags = MAKEFRAGMENT.out.nofiltfrags
    bws = ch_bw
    fragsize_plots = PLOTSIZE.out.fragsize_plots
    fragsize_txt = PLOTSIZE.out.length_counts
    tss_out = ch_tss

}

output {
    merged_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    cutadapt_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/cutadapt" }
    }
    cutadapt_logs {
        path { meta, _f1 -> "${meta.id}/cutadapt" }
    }
    bowtie2_logs {
        path { meta, _f1 -> "${meta.id}/bowtie2" }
    }
    markdup_bams {
        path { meta, _f1 -> "${meta.id}/bowtie2" }
    }
    filtdup_bams {
        path { meta, _f1 -> "${meta.id}/bowtie2" }
    }
    out_peaks {
        path { meta, _f1 -> "${meta.id}/macs2" }
    }
    filtfrags {
        path { meta, _f1 -> "${meta.id}/fragments" }
    }
    nofiltfrags {
        path { meta, _f1 -> "${meta.id}/fragments" }
    }
    fragsize_plots {
        path { meta, _f1 -> "${meta.id}/fragments" }
    }
    fragsize_txt {
        path { meta, _f1 -> "${meta.id}/fragments" }
    }
    bws {
        path { meta, _f1 -> "${meta.id}/bws" }
    }
    tss_out {
        path { meta, _f1, _f2 -> "${meta.id}/fragments" }
    }

}
