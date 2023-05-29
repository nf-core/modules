process SHAPEIT5_SWITCH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::shapeit5=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:1.0.0--h0c8ee15_0':
        'biocontainers/shapeit5:1.0.0--h0c8ee15_0'}"

    input:
        tuple val(meta) , path(estimate), path(estimate_index), val(region), path(pedigree)
        tuple val(meta2), path(truth)   , path(truth_index)
        tuple val(meta3), path(freq)    , path(freq_index)

    output:
        tuple val(meta), path("*.txt.gz"), emit: errors
        path "versions.yml"              , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args         = task.ext.args   ?: ''
        def prefix       = task.ext.prefix ?: "${meta.id}"
        def freq_cmd     = freq            ? "--frequency ${freq}"   : ""
        def pedigree_cmd = pedigree        ? "--pedigree ${pedigree}": ""

        """
        SHAPEIT5_switch \\
            $args \\
            --estimation $estimate \\
            --region $region \\
            --validation $truth \\
            $freq_cmd \\
            $pedigree_cmd \\
            --thread $task.cpus \\
            --output ${prefix}

        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                shapeit5: "\$(SHAPEIT5_switch | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
        END_VERSIONS
        """
}
