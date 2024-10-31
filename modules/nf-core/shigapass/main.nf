process SHIGAPASS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigapass:1.5.0--hdfd78af_0':
        'biocontainers/shigapass:1.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: report
    tuple val(meta), path("*_ShigaPass_Flex_summary.tsv"), optional: true, emit: flex_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # ShigaPass does not accept compressed FASTA files
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Convert our genome path to a file with a path in it
    ls $fasta_name > ${fasta_name}_tmp.txt

    # Run ShigaPass
    ShigaPass.sh \\
        -l ${fasta_name}_tmp.txt \\
        $args \\
        -p /usr/local/share/shigapass-1.5.0/db \\
        -t $task.cpus \\
        -o ${prefix}

    # Remove the temporary file from above
    rm ${fasta_name}_tmp.txt

    # Convert to tab delimited and move to the pwd
    sed 's/;/\t/g' ${prefix}/ShigaPass_summary.csv > ${prefix}.tsv

    # Convert to tab delimited and move to the pwd
    [ ! -f ${prefix}/ShigaPass_Flex_summary.csv ] || sed 's/;/\t/g' ${prefix}/ShigaPass_Flex_summary.csv > ${prefix}_Flex_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh --version 2>&1) | sed 's/^.*ShigaPass version //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh --version 2>&1) | sed 's/^.*ShigaPass version //' )
    END_VERSIONS
    """
}
