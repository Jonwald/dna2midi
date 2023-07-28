#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * pipeline input parameters --can be overiiden at command line
 */

params.reads = "/home/ec2-user/fastq/*_R{1,2}.fastq.gz"
params.genome_index = "/data/jyoung/rnaseq_nf/refs/Homo_sapiens.GRCh38.dna.primary_assembly.ERCC.Homo_sapiens.GRCh38.91"
params.rrna_index = "/home/ec2-user/refs/human_rrna"
params.gtf = "/data/jyoung/rnaseq_nf/refs/Homo_sapiens.GRCh38.91.gtf"
params.hk_genes = "/data/jyoung/rnaseq_nf/refs/hk_gene_ids.txt"
params.outdir = "results"
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

	container "791141132915.dkr.ecr.us-east-2.amazonaws.com/fastqc"

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

process BWAMEM {

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

workflow {

        Channel
                .fromFilePairs( params.reads, checkIfExists: true)
                .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
                .set { read_pairs_ch }

        fastqc_ch = FASTQC(read_pairs_ch)
        bwa_ch = BWAMEM(read_pairs_ch, params.rrna_index)

    }


