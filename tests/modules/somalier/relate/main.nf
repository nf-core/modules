#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_RELATE } from '../../../../modules/somalier/relate/main.nf'

workflow test_somalier_relate {

    input = [ [ id:'cohort', single_end:false ], // meta map
        [file("./tests/modules/somalier/relate/normal.somalier", checkIfExists: true),
        file("./tests/modules/somalier/relate/tumour.somalier", checkIfExists: true)]
    ]

    SOMALIER_RELATE (input,[],[])
}


workflow test_somalier_relate_ped_groups {

    input = [ [ id:'cohort', single_end:false ], // meta map
        [file("./tests/modules/somalier/relate/normal.somalier", checkIfExists: true),
        file("./tests/modules/somalier/relate/tumour.somalier", checkIfExists: true)]
    ]
    groups = file("./tests/modules/somalier/relate/groups.txt", checkIfExists: true)
    ped = file("./tests/modules/somalier/relate/families.ped", checkIfExists: true)

    SOMALIER_RELATE (input,groups,ped)
}
