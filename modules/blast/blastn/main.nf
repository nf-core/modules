process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0'
    } else {
        container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'
    }

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.blastn.txt'), emit: txt
    path "versions.yml"                  , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/.ndb//'`
    blastn \\
        -num_threads $task.cpus \\
        -db \$DB \\
        -query $fasta \\
        $args \\
        -out ${prefix}.blastn.txt
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
