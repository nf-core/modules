

process MTMALIGN_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtm-align:20220104--h4ac6f70_0':
        'biocontainers/mtm-align:20220104--h4ac6f70_0' }"

    input:
    tuple val(meta), path('*.pdb', arity: '2..*')

    output:
    tuple val(meta), path("./mTM_result/${prefix}.aln")  , emit: alignment
    tuple val(meta), path("./mTM_result/${prefix}.pdb")    , emit: structure
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls *.pdb | sed s/\\ /\\n/ > input_list.txt
    mtm-align -i input_list.txt -o ${prefix}.pdb
    # -o does not affect the fasta naming, so move it to the new name
    mv ./mTM_result/result.fasta ./mTM_result/${prefix}.aln

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\\(Version [[:digit:]]*\\)" | grep -oe "[[:digit:]]*" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir mTM_result
    touch mTM_result/${prefix}.aln
    touch mTM_result/${prefix}.pdb

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\\(Version [[:digit:]]*\\)" | grep -oe "[[:digit:]]*" )
    END_VERSIONS
    """
}
