//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    ch_vcf               // channel: [ val(meta), path(vcf) ]
    ch_fasta             // channel: [ path(fasta)]      fasta to use (optionnal)
    ch_vep_genome        // channel: [ val(genome)]      genome to use
    ch_vep_species       // channel: [val(species)]      species to use
    ch_vep_cache_version // channel: [val(species)]      cache version to use
    ch_vep_cache         // channel: [ path(cache)]      path: (/path/to/vep/cache (optionnal)
    ch_vep_extra_files   // channel: [ path(file1), path(file2)......] (optionnal additional files)

    main:
    ch_versions = Channel.empty()

    ENSEMBLVEP_VEP(ch_vcf, ch_vep_genome, ch_vep_species, ch_vep_cache_version, ch_vep_cache, ch_fasta, ch_vep_extra_files)
    TABIX_TABIX(ENSEMBLVEP_VEP.out.vcf)

    ch_vcf_tbi = ENSEMBLVEP_VEP.out.vcf.join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), path(vcf), path(tbi) ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), path(json) ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), path(tab) ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ path(html) ]
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
