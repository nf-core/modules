#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SHAPEIT5_LIGATE                   } from '../../../../../modules/nf-core/shapeit5/ligate/main.nf'
include { SHAPEIT5_PHASECOMMON              } from '../../../../../modules/nf-core/shapeit5/phasecommon/main.nf'
include { BCFTOOLS_VIEW                     } from '../../../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX                    } from '../../../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX2 } from '../../../../../modules/nf-core/bcftools/index/main.nf'

workflow test_shapeit5_ligate {
    input_vcf = [
        [ id:'NA12878_1X', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/NA12878.chr21.s.1x.vcf.gz.csi", checkIfExists: true),
    ]
    
    ref_panel = Channel.of([
        [ id:'REF_1000GP', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/glimpse/1000GP.chr21.noNA12878.s.bcf.csi", checkIfExists: true)
    ]).collect()

    scaffold  = Channel.of([[],[],[]]).collect()
    map       = Channel.of([[],[]]).collect()

    BCFTOOLS_VIEW ( input_vcf, [], [], [] )
    BCFTOOLS_INDEX ( BCFTOOLS_VIEW.out.vcf )

    region = Channel.of("chr21:16600000-16750000", "chr21:16650000-16800000")
    sample = Channel.of([[]])

    phase_input = Channel.of([[ id:'NA12878_1X', single_end:false ]])
                        .combine(BCFTOOLS_VIEW.out.vcf.collect().map{it[1]})
                        .combine(BCFTOOLS_INDEX.out.csi.collect().map{it[1]})
                        .combine(region)
                        .combine(sample)

    SHAPEIT5_PHASECOMMON ( phase_input, ref_panel, scaffold, map )
    
    BCFTOOLS_INDEX2 ( SHAPEIT5_PHASECOMMON.output.phased_variant )

    ligate_input = SHAPEIT5_PHASECOMMON.output.phased_variant.groupTuple()
                                        .join(BCFTOOLS_INDEX2.out.csi.groupTuple())

    ligate_input.view()
    SHAPEIT5_LIGATE ( ligate_input )
}
