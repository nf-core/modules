

process MTMALIGN_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mtm-align=20220104"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtm-align:20220104--h4ac6f70_0':
        'biocontainers/mtm-align:20220104--h4ac6f70_0' }"

    input:
    tuple val(meta), path(pdbs, name='*.pdb', arity: '1..*')

    output:
    tuple val(meta), path("./mTM_result/result.fasta")  , emit: alignment
    tuple val(meta), path("./mTM_result/result.pdb")    , emit: structure
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // unused, as mtmalign does not have any optional arguments
    """
    ls *.pdb | sed s/\ /\n/ > input_list.txt
    mtm-align -i input_list.txt

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\(Version [[:digit:]]*\)" | grep -oe "[[:digit:]]*" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir mTM_result
    touch mTM_result/result.fasta

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\(Version [[:digit:]]*\)" | grep -oe "[[:digit:]]*" )
    END_VERSIONS
    """
}
