#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_PHASE_SHAPEIT5 } from '../../../../subworkflows/nf-core/vcf_phase_shapeit5/main.nf'
include { BCFTOOLS_VIEW      } from '../../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX     } from '../../../../modules/nf-core/bcftools/index/main.nf'

workflow test_vcf_phase_shapeit5_panel {
    
    ref_panel = Channel.of([
        [ id:'1000GP',region:"chr21:16600115-16799989"],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true),
        [],
        "chr21:16600115-16799989"
    ])

    VCF_PHASE_SHAPEIT5 ( ref_panel,
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[]]).collect()
    )
}

workflow test_vcf_phase_shapeit5_ind {
    input_vcf = Channel.of([
        [ id:'NA12878', region:"chr21:16600115-16799989"], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi",
            checkIfExists: true)
    ])
    ref_panel = Channel.of([
        [ panel:'1000GP'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf",
            checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi",
            checkIfExists: true)        
    ])
    ch_map = Channel.of([
        [ map:'chr21.b38'],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/chr21.b38.gmap.gz", checkIfExists: true)
    ])
    BCFTOOLS_VIEW(input_vcf, [], [], [])
    BCFTOOLS_INDEX(BCFTOOLS_VIEW.out.vcf)

    input_vcf = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_INDEX.out.csi)
        .combine(Channel.of([[],"chr21:16600115-16799989"]))

    VCF_PHASE_SHAPEIT5 ( input_vcf,
        ref_panel.collect(),
        Channel.of([[],[],[]]).collect(),
        ch_map.collect()
    )

}
