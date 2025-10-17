process LINKS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/links:2.0.1--h4ac6f70_5':
        'biocontainers/links:2.0.1--h4ac6f70_5' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(reads)

    output:
    tuple val(meta), path("*.log"),                         emit: log
    tuple val(meta), path("*.pairing_distribution.csv"),    emit: pairing_distribution,  optional: true
    tuple val(meta), path("*.pairing_issues"),              emit: pairing_issues
    tuple val(meta), path("*.scaffolds"),                   emit: scaffolds_csv
    tuple val(meta), path("*.scaffolds.fa"),                emit: scaffolds_fasta
    tuple val(meta), path("*.bloom"),                       emit: bloom
    tuple val(meta), path("*.gv"),                          emit: scaffolds_graph
    tuple val(meta), path("*.assembly_correspondence.tsv"), emit: assembly_correspondence
    tuple val(meta), path("*.simplepair_checkpoint.tsv"),   emit: simplepair_checkpoint, optional: true
    tuple val(meta), path("*.tigpair_checkpoint.tsv"),      emit: tigpair_checkpoint
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Currently LINKS does not support more than 4 threads
    def nthreads = "${task.cpus}" < 4 ? "${task.cpus}" : 4
    def args = task.ext.args ?: ""
    """
    if [[ ${assembly} == *.gz ]];
    then
        gzip -dc ${assembly} > assembly.fa
    else
        ln -s ${assembly} assembly.fa
    fi

    for read_file in ${reads};
        do
            if [[ \$read_file == *.gz ]];
            then 
                gzip -dc \$read_file > \$(basename \$read_file .gz)
                echo \$(basename \$read_file .gz) >> readfile.fof
            else
                echo \$read_file >> readfile.fof
            fi
        done

    LINKS -f assembly.fa \\
        -s readfile.fof \\
        -j $nthreads \\
        -b ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LINKS: \$(echo \$(LINKS | grep -o 'LINKS v.*' | sed 's/LINKS v//'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.log
    touch ${prefix}.pairing_distribution.csv
    touch ${prefix}.pairing_issues
    touch ${prefix}.scaffolds
    touch ${prefix}.scaffolds.fa
    touch ${prefix}.bloom
    touch ${prefix}.gv
    touch ${prefix}.assembly_correspondence.tsv
    touch ${prefix}.simplepair_checkpoint.tsv
    touch ${prefix}.tigpair_checkpoint.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LINKS: \$(echo \$(LINKS | grep -o 'LINKS v.*' | sed 's/LINKS v//'))
    END_VERSIONS
    """
    }
