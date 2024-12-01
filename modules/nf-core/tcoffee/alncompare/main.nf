process TCOFFEE_ALNCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
        'biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

    input:
    tuple val(meta), path(msa), path(ref_msa)

    output:
    tuple val(meta), path("*.scores"), emit: scores
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ? task.ext.prefix : "${msa.baseName}"
    def args = task.ext.args ? task.ext.args.contains('compare_mode') ? task.ext.args  : (task.ext.args +'-compare_mode tc ' ) : '-compare_mode tc '
    def metric_name = args.split('compare_mode ')[1].split(' ')[0]
    def header = meta.keySet().join(",")
    def values = meta.values().join(",")

    """
    # check whether it is compressed
    if [[ "${msa}" == *.gz ]]; then
        unpigz -c ${msa} > uncompressed_msa.fa
    else
        ln ${msa} uncompressed_msa.fa
    fi

    export TEMP='./'
    t_coffee -other_pg aln_compare \
        -al1 ${ref_msa} \
        -al2 uncompressed_msa.fa \
        ${args} \
        | grep -v "seq1" | grep -v '*' | \
        awk '{ print \$4}' ORS="\t" \
        >> "scores.txt"

    # Add metadata info to output file
    echo "${header},${metric_name}" > "${prefix}.scores"

    # Add values
    scores=\$(awk '{sub(/[[:space:]]+\$/, "")} 1' scores.txt | tr -s '[:blank:]' ',')
    echo "${values},\$scores" >> "${prefix}.scores"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch "${prefix}.scores"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
