// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls

include { GATK4_CALCULATECONTAMINATION } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS      } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES     } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_MUTECT2                } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input // channel: [ val(meta), [ input ], [ input_index ], [] ]
    ch_fasta // channel: /path/to/reference/fasta
    ch_fai // channel: /path/to/reference/fasta/index
    ch_dict // channel: /path/to/reference/fasta/dictionary
    ch_alleles // channel: /path/to/alleles
    ch_alleles_tbi // channel: /path/to/alleles/index
    ch_germline_resource // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_panel_of_normals // channel: /path/to/panel/of/normals
    ch_panel_of_normals_tbi // channel: /path/to/panel/of/normals/index
    ch_interval_file // channel: /path/to/interval/file

    main:
    //Perform variant calling using mutect2 module in tumor single mode.
    GATK4_MUTECT2(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_alleles,
        ch_alleles_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
    )

    //Generate pileup summary table using getpileupsummaries.

    ch_pileup_input = ch_input
        .combine(ch_interval_file)
        .map { meta, input_file, input_index, _which_norm, intervals ->
            [meta, input_file, input_index, intervals]
        }

    GATK4_GETPILEUPSUMMARIES(
        ch_pileup_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
    )

    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //[] is a placeholder for the optional input where the matched normal sample would be passed in for tumor-normal samples, which is not necessary for this workflow.
    ch_pileup = GATK4_GETPILEUPSUMMARIES.out.table.collect().map { meta, table -> [meta, table, []] }


    GATK4_CALCULATECONTAMINATION(ch_pileup)

    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.

    ch_vcf = GATK4_MUTECT2.out.vcf.collect()
    ch_tbi = GATK4_MUTECT2.out.tbi.collect()
    ch_stats = GATK4_MUTECT2.out.stats.collect()
    ch_segment = GATK4_CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = GATK4_CALCULATECONTAMINATION.out.contamination.collect()

    ch_filtermutect_in = ch_vcf
        .join(ch_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_stats, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_segment, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_contamination, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi, stats, segment, contamination -> [meta, vcf, tbi, stats, [], segment, contamination, []] }

    GATK4_FILTERMUTECTCALLS(ch_filtermutect_in, ch_fasta, ch_fai, ch_dict)

    emit:
    contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination.collect() // channel: [ val(meta), [ contamination ] ]
    filtered_stats      = GATK4_FILTERMUTECTCALLS.out.stats.collect() // channel: [ val(meta), [ stats ] ]
    filtered_tbi        = GATK4_FILTERMUTECTCALLS.out.tbi.collect() // channel: [ val(meta), [ tbi ] ]
    filtered_vcf        = GATK4_FILTERMUTECTCALLS.out.vcf.collect() // channel: [ val(meta), [ vcf ] ]
    mutect2_tbi         = GATK4_MUTECT2.out.tbi.collect() // channel: [ val(meta), [ tbi ] ]
    mutect2_stats       = GATK4_MUTECT2.out.stats.collect() // channel: [ val(meta), [ stats ] ]
    mutect2_vcf         = GATK4_MUTECT2.out.vcf.collect() // channel: [ val(meta), [ vcf ] ]
    pileup_table        = GATK4_GETPILEUPSUMMARIES.out.table.collect() // channel: [ val(meta), [ table ] ]
    segmentation_table  = GATK4_CALCULATECONTAMINATION.out.segmentation.collect() // channel: [ val(meta), [ segmentation ] ]
}
