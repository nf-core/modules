include { RTGTOOLS_VCFEVAL      } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { HAPPY_HAPPY           } from '../../../modules/nf-core/happy/happy/main'

workflow VCF_VALIDATE_SNPS {

    take:
    ch_vcf                          // [mandatory] channel: [ meta, vcf, tbi, truth_vcf, truth_tbi, bed ]
    ch_fasta                        // [happy only] channel: [ meta, fasta ]
    ch_fasta_fai                    // [happy only] channel: [ meta, fasta_fai ]
    ch_vcfeval_sdf                  // [vcfeval only] channel: [ meta, sdf ]
    ch_happy_false_positive_regions // [optional] channel: [ meta, false_positives_bed ]
    ch_happy_stratification_tsv     // [optional] channel: [ meta, stratification_tsv ]
    ch_happy_stratification_beds    // [optional] channel: [ meta, [stratification_beds] ]
    tools                           // [mandatory] value: A comma-delimited list of the tools to use for validation (happy,vcfeval)

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    ch_happy_vcfs = Channel.empty()
    ch_happy

    list_tools = tools.tokenize(",")

    if("happy" in list_tools){
        happy_input = ch_vcf
            .map { meta, vcf, tbi, truth_vcf, truth_tbi, bed ->
                [ meta, vcf, truth_vcf, bed, [] ] // TODO write logic for regions and target bed
            }
        HAPPY_HAPPY (
            happy_input,
            ch_fasta,
            ch_fasta_fai,
            ch_happy_false_positive_regions,
            ch_happy_stratification_tsv,
            ch_happy_stratification_beds
        )
        ch_versions = ch_versions.mix(HAPPY_HAPPY.out.versions.first())
        ch_happy_vcfs = HAPPY_HAPPY.out.vcf.join(HAPPY_HAPPY.out.tbi)
    }

    if("vcfeval" in list_tools){
        RTGTOOLS_VCFEVAL(
            ch_vcf,
            ch_vcfeval_sdf
        )
        ch_versions = ch_versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

        // TODO add RTGTOOLS_ROCPLOT to generate images for roc TSVs
    }

    emit:
    happy_vcfs      = ch_happy_vcfs            // channel: [ meta, vcf, tbi ]


    versions = ch_versions                     // channel: [ versions.yml ]
}

