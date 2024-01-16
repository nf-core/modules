

process MTMALIGN_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5bcf71dc66dac33d8e003c5e78043b80f5c7f269':
        'biocontainers/mulled-v2-5bcf71dc66dac33d8e003c5e78043b80f5c7f269' }"

    input:
    tuple val(meta), path('*.pdb', arity: '2..*')

    output:
    tuple val(meta), path("./mTM_result/${prefix}.aln.gz")  , emit: alignment
    tuple val(meta), path("./mTM_result/${prefix}.pdb.gz")    , emit: structure
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

    # compress both output files
    pigz -p ${task.cpus} ./mTM_result/${prefix}.aln ./mTM_result/${prefix}.pdb

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
    touch mTM_result/${prefix}.aln.gz
    touch mTM_result/${prefix}.pdb.gz

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\\(Version [[:digit:]]*\\)" | grep -oe "[[:digit:]]*" )
    END_VERSIONS
    """
}
