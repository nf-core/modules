//
// Run GATK haplotypecaller for all input samples, merge them with genomicsdbimport, perform joint genotyping with genotypeGVCFS and recalibrate with variantrecalibrator & applyvqsr.
//

include { GATK4_HAPLOTYPECALLER } from '../../../modules/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../modules/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTRECALIBRATOR } from '../../../modules/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR } from '../../../modules/gatk4/applyvqsr/main'

workflow GATK_JOINT_GERMLINE_VARIANT_CALLING {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], [] ]
    run_haplotc      // channel: true/false run haplotypecaller portion of workflow or skip to genomicsdbimport when false
    run_vqsr         // channel: true/false run vqsr portion of subworkflow
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    sites            // channel: /path/to/known/sites/file
    sites_index      // channel: /path/to/known/sites/index
    joint_id         // channel: joint id to replace individual sample ids with
    joint_intervals  // channel: joint intervals to be used for genomicsdbimport and genotypegvcfs
    allelespecific   // channel: true/false run allelespecific mode of vqsr modules
    resources        // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    annotation       // channel: [annotations, to, use, for, variantrecal, filtering]
    mode             // channel: which mode to run variantrecal: SNP/INDEL/BOTH
    create_rscript   // channel: true/false whether to generate rscript plots in variantrecal
    truthsensitivity // channel: 0-100.0 truthsensitivity cutoff for applyvqsr

    main:
    ch_versions = Channel.empty()

    // haplotypecaller can be skipped if input samples are already in gvcf format, essentially making the subworkflow joint genotyping.
    if (run_haplotc) {
        haplotc_input = channel.from(input)
        //
        //Perform variant calling using haplotypecaller module. Additional argument "-ERC GVCF" used to run in gvcf mode.
        //
        GATK4_HAPLOTYPECALLER ( haplotc_input, fasta, fai, dict, sites, sites_index)

        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())
        ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.collect{it[1]}.toList()
        ch_index = GATK4_HAPLOTYPECALLER.out.tbi.collect{it[1]}.toList()

    } else {
        // if haplotypecaller is skipped, this channels the input to genomicsdbimport instead of the output vcfs and tbis that normally come from haplotypecaller
        direct_input = channel.from(input)
        ch_vcf = direct_input.collect{it[1]}.toList()
        ch_index = direct_input.collect{it[2]}.toList()
    }

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    gendb_input = Channel.of([[ id:joint_id ]]).combine(ch_vcf).combine(ch_index).combine([joint_intervals]).combine(['']).combine([dict])
    gendb_input.view()
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )

    ch_versions       = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    //
    //Joint genotyping performed using GenotypeGVCFs
    //
    ch_genotype_in = GATK4_GENOMICSDBIMPORT.out.genomicsdb.collect()
    //[] is a placeholder for the input where the vcf tbi file would be passed in for non-genomicsdb workspace runs, which is not necessary for this workflow as it uses genomicsdb workspaces.
    ch_genotype_in.add([])
    genotype_input = ch_genotype_in.combine([joint_intervals])
    GATK4_GENOTYPEGVCFS ( genotype_input, fasta, fai, dict, sites, sites_index)

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    // setting run_vqsr to false skips the VQSR process, for if user does not wish to perform VQSR,
    // or want to run the hard filtering recommended by gatk best practices for runs with a low number of samples instead.
    if (run_vqsr) {
        //
        //Perform first step in VQSR using VariantRecalibrator
        //
        ch_gvcf       = GATK4_GENOTYPEGVCFS.out.vcf.collect()
        ch_gtbi       = GATK4_GENOTYPEGVCFS.out.tbi.collect()
        ch_vrecal_in  = ch_gvcf.combine(ch_gtbi, by: 0)

        GATK4_VARIANTRECALIBRATOR ( ch_vrecal_in, fasta, fai, dict, allelespecific, resources, annotation, mode, create_rscript )

        ch_versions   = ch_versions.mix(GATK4_VARIANTRECALIBRATOR.out.versions)

        //
        //Perform second step in VQSR using ApplyVQSR
        //
        ch_recal      = GATK4_VARIANTRECALIBRATOR.out.recal.collect()
        ch_idx        = GATK4_VARIANTRECALIBRATOR.out.idx.collect()
        ch_tranches   = GATK4_VARIANTRECALIBRATOR.out.tranches.collect()
        ch_vqsr_in    = ch_vrecal_in.combine(ch_recal, by: 0).combine(ch_idx, by: 0).combine(ch_tranches, by: 0)

        GATK4_APPLYVQSR ( ch_vqsr_in, fasta, fai, dict, allelespecific, truthsensitivity, mode )

        ch_versions   = ch_versions.mix(GATK4_APPLYVQSR.out.versions)

    }

    emit:
    versions       = ch_versions                                     // channel: [ versions.yml ]
    genomicsdb     = GATK4_GENOMICSDBIMPORT.out.genomicsdb.collect() // channel: [ val(meta), [ genomicsdb ] ]
    genotype_vcf   = GATK4_GENOTYPEGVCFS.out.vcf.collect()           // channel: [ val(meta), [ vcf ] ]
    genotype_index = GATK4_GENOTYPEGVCFS.out.vcf.collect()           // channel: [ val(meta), [ tbi ] ]
    haplotc_vcf    = run_haplotc ? GATK4_HAPLOTYPECALLER.out.vcf.collect()       : [] // channel: [ val(meta), [ vcf ] ]
    haplotc_index  = run_haplotc ? GATK4_HAPLOTYPECALLER.out.tbi.collect()       : [] // channel: [ val(meta), [ tbi ] ]
    recal_file     = run_vqsr ? GATK4_VARIANTRECALIBRATOR.out.recal.collect()    : [] // channel: [ val(meta), [ recal ] ]
    recal_index    = run_vqsr ? GATK4_VARIANTRECALIBRATOR.out.idx.collect()      : [] // channel: [ val(meta), [ idx ] ]
    recal_tranches = run_vqsr ? GATK4_VARIANTRECALIBRATOR.out.tranches.collect() : [] // channel: [ val(meta), [ tranches ] ]
    vqsr_vcf       = run_vqsr ? GATK4_APPLYVQSR.out.vcf.collect()                : [] // channel: [ val(meta), [ vcf ] ]
    vqsr_index     = run_vqsr ? GATK4_APPLYVQSR.out.tbi.collect()                : [] // channel: [ val(meta), [ tbi ] ]
}
