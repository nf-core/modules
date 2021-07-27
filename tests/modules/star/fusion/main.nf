#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE               } from '../../../../modules/star/genomegenerate/main.nf' addParams( options: [args: '--genomeSAindexNbases 9'] )
include { STAR_ALIGN as STAR_FOR_STARFUSION } from '../../../../modules/star/align/main.nf'          addParams( options: [args: '--readFilesCommand zcat --outSAMtype BAM Unsorted --outReadsUnmapped None --twopassMode Basic --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30'] )
include { STAR_FUSION                       } from '../../../../modules/star/fusion/main.nf' addParams( options: [:] )

workflow test_star_alignment_paired_end_for_starfusion {
    input = [ [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    genome_dir  = Channel.value(file(params.starfusion_genome_dir))

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_FOR_STARFUSION ( input, STAR_GENOMEGENERATE.out.index, gtf )
    STAR_FUSION (STAR_FOR_STARFUSION.junction, genome_dir)
}
