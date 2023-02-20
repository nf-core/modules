process TIARA_TIARA {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "conda-forge::tiara=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiara:1.0.3' :
        'quay.io/biocontainers/tiara:1.0.3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("classification_*.out.{txt,txt.gz}")      , emit: classifications
    tuple val(meta), path("log_*.out.{txt,txt.gz}")                 , emit: log
    tuple val(meta), path("*.{fasta,fasta.gz}")                     , emit: fasta, optional: true
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tiara -i ${fasta} \
        -o ${prefix}.out.txt \
        --threads ${task.cpus} \
        ${args}

    ## Rename classification files
    if [ -f ${prefix}.out.txt ]; then
        mv ${prefix}.out.txt classification_${prefix}.out.txt
    else
        mv ${prefix}.out.txt.gz classification_${prefix}.out.txt.gz
    fi

    ## fix gzip flag weirdness and ensure consistent .fasta filename output
    ## check if fasta files are being output
    if echo "${args}" | grep -q "tf\|to-fasta"; then
        ## check if we've asked for gzip output, then rename files consistently
        if echo "${args}" | grep -q "gz"; then
            find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}.fasta.gz"' \\;
        else
            find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}.fasta"' \\;
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """
}
