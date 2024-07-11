process UNTAR {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:e1e4ff8dd129544f' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:10736d98b4d693d8' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$prefix"), emit: untar
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix    = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.baseName.toString().replaceFirst(/\.tar$/, ""))

    """
    mkdir $prefix

    ## Ensures --strip-components only applied when top level of tar contents is a directory
    ## If just files or multiple directories, place all in prefix
    if [[ \$(tar -tavf ${archive} | grep -o -P "/\$" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -C $prefix --strip-components 1 \\
            -xavf \\
            $args \\
            $archive \\
            $args2
    else
        tar \\
            -C $prefix \\
            -xavf \\
            $args \\
            $archive \\
            $args2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.toString().replaceFirst(/\.[^\.]+(.gz)?$/, ""))
    """
    mkdir ${prefix}
    ## Dry-run untaring the archive to get the files and place all in prefix
    if [[ \$(tar -tavf ${archive} | grep -o -P "/\$" | uniq | wc -l) -eq 1 ]]; then
        for i in `tar -tf ${archive}`;
        do
            if [[ \$(echo "\${i}" | grep -E "/\$") == "" ]];
            then
                touch \${i}
            else
                mkdir -p \${i}
            fi
        done
    else
        for i in `tar -tf ${archive}`;
        do
            if [[ \$(echo "\${i}" | grep -E "/\$") == "" ]];
            then
                touch ${prefix}/\${i}
            else
                mkdir -p ${prefix}/\${i}
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
