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

process FASTQ_TRANSFER {
        publishDir "$params.outdir/fastq_data", mode:'copy'

        input:
        tuple val(pair_id), path(reads)
        path bs_config   /* ## add this*/

        output:
        path "*.fastq.gz"

        script:
        """
        #bs  list project -c $bs_config | cut -f2
        # list runs | grep project id;
        # list datasets |grep "project id" |grep "run" | > list
        # for dataset in list; bs download fastq
        fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
        """
}

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

process STAR_ALIGN {

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
        tuple val(pair_id), path(bam), 
        path '*_logs'

        output:
        path 'multiqc_report.html'

        script:
        """
        multiqc .
        """
}


process MARK_DUPS {

        publishDir "$params.outdir/mark_dups", mode:'copy'

        input:
        tuple val(pair_id), path(bam), path(logs)

        output:
        tuple val(pair_id), path("*.bam"), path("*_marked_dup_metrics.txt")

        script:
        """
        java -jar picard.jar MarkDuplicates I=${bam[0]} \
        O=${pair_id}_Aligned.sortedByCoord.out.marked.bam \
        M=${pair_id}_marked_dup_metrics.txt \
        REMOVE_DUPLICATES=False \

        """

}

process REM_DUPS {

        publishDir "$params.outdir/rem_dups", mode:'copy'

        input:
        tuple val(pair_id), path(bam), path(logs)

        output:
        tuple val(pair_id), path("*.bam"), path("*_rem_dup_metrics.txt")

        script:
        """
        java -jar picard.jar MarkDuplicates I=${bam[0]} \
        O=${pair_id}_Aligned.sortedByCoord.out.removed.bam \
        M=${pair_id}_rem_dup_metrics.txt \
        REMOVE_DUPLICATES=True \

        """
}

process STRINGTIE {

        publishDir "$params.outdir/stringtie", mode:'copy'

        input:
        tuple val(pair_id), path(bam), path(logs)
        path gtf

        output:
        tuple val(pair_id), path("*.gtf"), path("*.tab"), path("*_t_data.ctab")

        script:
        """
        stringtie -e -B -p 4 --rf -c 0.001 -G $gtf -A ${pair_id}_nodedup_gene_abund.tab -o ${pair_id}_stringtie.gtf ${bam[0]}
        mv t_data.ctab ${pair_id}_nodedup_t_data.ctab

        """
}

process STRINGTIE_DEDUP {

        publishDir "$params.outdir/stringtie_dedup", mode:'copy'

        input:
        tuple val(pair_id), path(bam)
        path gtf

        output:
        path "*.tab"

        script:
        """
        /stringtie -e -B -p 4 --rf -c 0.001 -G $gtf -A ${pair_id}_dedup_gene_abund.tab -o ${pair_id}_stringtie.gtf ${bam[0]}
        mv t_data.ctab ${pair_id}_t_data.ctab

        """
}

process POST_PROCESSING {

        publishDir "$params.outdir/QC", mode:'copy'

        input:

        path "*"
        path hk_genes

        output:
        path "*_qc_metrics.txt"
        path "*.tsv"

        script:
        """
        ## compile gene matirces
        python /data/jyoung/rnaseq_nf/bin/merge_tsv.py -p \$PWD -e "_dedup_gene_abund.tab" -k "Gene ID" -c "Coverage" -o test_hkout.COV
        python /data/jyoung/rnaseq_nf/bin/merge_tsv.py -p \$PWD -e "_nodedup_gene_abund.tab" -k "Gene ID" -c "TPM" -o test_out_tmp
        python /data/jyoung/rnaseq_nf/bin/merge_tsv.py -p \$PWD -e "_nodedup_gene_abund.tab" -k "Gene ID" -c "FPKM" -o test_out_fpkm
        python /data/jyoung/rnaseq_nf/bin/merge_tsv.py -p \$PWD -e "_nodedup_t_data.ctab" -k "t_name" -c "FPKM" -o test_out_tr_fpkm

        ## get sample names
        for file in *.gtf;
         do sname=`basename \$file .gtf`;
         echo \$sname \$file;
         done > snames

        ## compile counts
        python /data/jyoung/rnaseq_nf/bin/prepDE.py -i snames

        ## unzip fastqc results
        mv **/*.zip .
        for file in *.zip; do unzip \$file; done

        ## compile qc file
        python /data/jyoung/rnaseq_nf/bin/parse_qc.py -f \$PWD

        """
}


workflow {
        
        Channel
                .fromFilePairs( params.reads, checkIfExists: true)
                .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
                .set { read_pairs_ch }

        fastqc_ch = FASTQC(read_pairs_ch)
        bwa_ch = RRNA_QUANT(read_pairs_ch, params.rrna_index)
        bam_ch = STAR_ALIGN(read_pairs_ch, params.genome_index)
        multiqc_ch = MULTIQC(bam_ch.collect(), fastqc_ch.collect())
        }

        markduplicates_ch = MARK_DUPS(bam_ch)
        removeduplicates_ch = REM_DUPS(bam_ch)
        stringtie_all_ch = STRINGTIE(markduplicates_ch, params.gtf)
        stringtie_dedup_ch = STRINGTIE_DEDUP(removeduplicates_ch, params.gtf)

        pp_1 = fastqc_ch | collect
        pp_2 = bwa_ch | collect
        pp_3 = markduplicates_ch | map { it[2] } | collect
        pp_4 = bam_ch | map { it[2] } | collect
        pp_5 = stringtie_dedup_ch | collect
        pp_6 = stringtie_all_ch | map { it[1,2,3] } | collect
        pp_in = pp_1 | mix (pp_2, pp_3, pp_4, pp_5, pp_6) | collect

        POST_PROCESSING(pp_in, params.hk_genes)
    }


