//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    ch_vcf              // channel: [ val(meta), path(vcf), [path(custom_file1), path(custom_file2)... (optional)]]
    ch_fasta            // channel: [ val(meta2), path(fasta) ] (optional)
    val_genome          //   value: genome to use
    val_species         //   value: species to use
    val_cache_version   //   value: cache version to use
    ch_cache            // channel: [ val(meta3), path(cache) ] (optional)
    ch_extra_files      // channel: [ path(file1), path(file2)... ] (optional)

    main:
    ENSEMBLVEP_VEP(
        ch_vcf,
        val_genome,
        val_species,
        val_cache_version,
        ch_cache,
        ch_fasta,
        ch_extra_files,
    )

    ch_vcf_tbi = ENSEMBLVEP_VEP.out.vcf.join(ENSEMBLVEP_VEP.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    emit:
    vcf_tbi = ch_vcf_tbi                // channel: [ val(meta), path(vcf), path(tbi) ]
    json    = ENSEMBLVEP_VEP.out.json   // channel: [ val(meta), path(json) ]
    tab     = ENSEMBLVEP_VEP.out.tab    // channel: [ val(meta), path(tab) ]
    reports = ENSEMBLVEP_VEP.out.report // channel: [ path(html) ]
}
