include { RTGTOOLS_VCFEVAL      } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { HAPPY_HAPPY           } from '../../../modules/nf-core/happy/happy/main'

workflow VCF_VALIDATE_SNPS {

    take:
    ch_vcf // channel: [ meta, vcf, tbi, truth_vcf, truth_tbi, bed ]
    ch_fasta // channel: []
    tools  // A comma-delimited list of the tools to use for validation (happy,vcfeval)

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    list_tools = tools.tokenize(",")

    if("happy" in list_tools){
        happy_input = ch_vcf
            .map { meta, vcf, tbi, truth_vcf, truth_tbi, bed ->
                [ meta, vcf, truth_vcf, bed, [] ] // TODO write logic for regions and target bed
            }
        HAPPY_HAPPY (

        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    }


    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

