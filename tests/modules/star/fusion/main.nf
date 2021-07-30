#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE               } from '../../../../modules/star/genomegenerate/main.nf' addParams( options: [args: '--genomeSAindexNbases 9'] )
include { STAR_ALIGN                        } from '../../../../modules/star/align/main.nf'          addParams( options: [args: '--readFilesCommand zcat --outSAMtype BAM Unsorted --outReadsUnmapped None --twopassMode Basic --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30'] )
include { STAR_FUSION                       } from '../../../../modules/star/fusion/main.nf'         addParams( options: [:] )

workflow test_star_alignment_paired_end_for_starfusion {

    params.starfusion_genome_dir = "${launchDir}/ctat_genome_lib_build_dir/"

    def input = [ [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true) ]
            ]
    def fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    def gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    def genome_dir  = file(params.starfusion_genome_dir, checkIfExists: false) // Tested locally

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf )
    STAR_FUSION (STAR_ALIGN.out.junction, genome_dir)
}
