#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * pipeline input parameters --can be overiiden at command line
 */

params.reads = "s3://alm-dx-test/fastq/*_R{1,2}.fastq.gz"
params.genome_index = "s3://alm-dx-test/Homo_sapiens.GRCh38.dna.primary_assembly.ERCC.Homo_sapiens.GRCh38.91.tar.gz"
params.rrna_index = "s3://alm-dx-test/human_rrna"
params.gtf = "s3://alm-dx-test/Homo_sapiens.GRCh38.91.gtf.gz"
params.hk_genes = "/data/jyoung/rnaseq_nf/refs/hk_gene_ids.txt"
params.outdir = "s3://alm-pl-tmp"
params.version = "3.0.0b"

log.info """
    A L M A C    R N A - S E Q    P I P E L I N E
    =============================================
    version   : ${params.version}
    reads     : ${params.reads}
    ref genome: ${params.genome_index}
    ref GTF   : ${params.gtf}
    outdir    : ${params.outdir}
    """
    .stripIndent()


/* read_pairs.into{fastqc_input; rrna_input; star_input}*/

process FASTQC {
        publishDir "$params.outdir/fastqc", mode:'copy'

        input:
        tuple val(pair_id), path(reads)

        output:
        path "fastqc_${pair_id}_logs"

        script:
        """
        mkdir fastqc_${pair_id}_logs
        fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
        """
}

process RRNA_QUANT {

        publishDir "$params.outdir/bwamem", mode:'copy'
        container "791141132915.dkr.ecr.us-east-2.amazonaws.com/rrna_quant"

        input:
        tuple val(pair_id), path(reads)
        path rrna_index

        output:
        path "*.txt"

        script:
        """
        bash /bin/rrna_quant.sh "$rrna_index" ${reads[0]} ${reads[1]} ${pair_id}_rrna_stats.txt

        """
}

process STAR {

        publishDir "$params.outdir/STAR", mode:'copy'

        input:
        tuple val(pair_id), path(reads)
        path index

        output:
        tuple val(pair_id), path("*.bam"), path("*_Log.final.out")

        script:
        """
        tar -xzf $index
        folname="\$(basename $index .tar.gz)"
        /data/jyoung/rnaseq_nf/bin/STAR-2.6.1a/source/./STAR --runThreadN 8 --genomeDir \
        "\$folname" \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 1 \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMismatchNmax 3 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --genomeLoad NoSharedMemory \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix ${pair_id}_
        """
}

process MULTIQC {

        publishDir "$params.outdir/multiqc", mode:'copy'

        input:
        tuple val(pair_id), path(bam)
        path '*_logs'

        output:
        path 'multiqc_report.html'

        script:
        """
        multiqc .
        """
}

workflow {

        Channel
                .fromFilePairs( params.reads, checkIfExists: true)
                .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
                .set { read_pairs_ch }

        fastqc_ch = FASTQC(read_pairs_ch)
        bwa_ch = RRNA_QUANT(read_pairs_ch, params.rrna_index)
        star_ch = STAR_ALIGN(read_pairs_ch, params.genome_index)
        MULTIQC(bam_ch.collect(), fastqc_ch.collect())
    }


