process IPHOP_PREDICT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iphop:1.3.2--pyhdfd78af_0':
        'biocontainers/iphop:1.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path iphop_db

    output:
    tuple val(meta), path("*.Host_prediction_to_genus_m*.csv")  , emit: iphop_genus
    tuple val(meta), path("*.Host_prediction_to_genome_m*.csv") , emit: iphop_genome
    tuple val(meta), path("*.Detailed_output_by_tool.csv")      , emit: iphop_detailed_output
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export PERL5LIB=/usr/local/lib/perl5/site_perl/5.22.0

    INPUT=$fasta
    if [[ $fasta == *.gz ]]
    then
        gunzip -c $fasta > ${prefix}.fasta
        INPUT=${prefix}.fasta
    fi

    iphop \\
        predict \\
        --fa_file \$INPUT \\
        --out_dir iphop_results \\
        --db_dir $iphop_db \\
        --num_threads $task.cpus \\
        $args

    for FILE in iphop_results/Host_prediction_to_*.csv; do
        mv \${FILE} ./${prefix}.\${FILE##*/}
    done

    mv iphop_results/Detailed_output_by_tool.csv ./${prefix}.Detailed_output_by_tool.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_score = args.contains('--min_score') ? args.split('--min_score ')[1] : '90'
    """
    mkdir -p iphop_results
    touch ${prefix}.Host_prediction_to_genus_m${min_score}.csv
    touch ${prefix}.Host_prediction_to_genome_m${min_score}.csv
    touch ${prefix}.Detailed_output_by_tool.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """
}
