#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_DOWNLOAD } from '../../../../../modules/nf-core/ensemblvep/download/main'

vep_cache_version = "110"
vep_genome = "WBcel235"
vep_species = "caenorhabditis_elegans"
vep_cache_input = Channel.of([[id:"${vep_cache_version}_${vep_genome}"], vep_genome, vep_species, vep_cache_version])

workflow test_ensemblvep_download {
    ENSEMBLVEP_DOWNLOAD(vep_cache_input)
}
