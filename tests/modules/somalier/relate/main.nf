#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_RELATE } from '../../../../modules/somalier/relate/main.nf'


workflow test_somalier_relate {

    input = [ [ id:'cohort', single_end:false ], // meta map
        [ file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/normal.somalier", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/tumour.somalier", checkIfExists: true) ]
    ]

    SOMALIER_RELATE (input,[],[])
}


workflow test_somalier_relate_ped_groups {

input = [ [ id:'cohort', single_end:false ], // meta map
    [ file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/normal.somalier", checkIfExists: true),
    file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/tumour.somalier", checkIfExists: true) ]
]

    groups = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/groups.txt", checkIfExists: true)
    ped = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/family.ped", checkIfExists: true)

    SOMALIER_RELATE (input,groups,ped)
}
