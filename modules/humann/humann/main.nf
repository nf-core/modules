// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HUMANN_HUMANN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::humann=3.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/humann:3.0.0--pyh5e36f6f_1"
    } else {
        container "quay.io/biocontainers/humann:3.0.0--pyh5e36f6f_1"
    }

    input:
    tuple val(meta), path(input), path(profile)
    path nucleotide_db
    path protein_db

    output:
    tuple val(meta), path("*_genefamilies.tsv.gz") , emit: genefamilies
    tuple val(meta), path("*_pathabundance.tsv.gz"), emit: pathabundance
    tuple val(meta), path("*_pathcoverage.tsv.gz") , emit: pathcoverage
    tuple val(meta), path("*.log")                 , emit: log
    path "versions.yml"                            , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    if [[ $nucleotide_db == *.tar.gz ]]; then
        mkdir nuc_database
        tar -xzvf $nucleotide_db -C nuc_database
    else
        mv $nucleotide_db nuc_database
    fi
    if [[ $protein_db == *.tar.gz ]]; then
        mkdir prot_database
        tar -xzvf $protein_db -C prot_database
    else
        mv $protein_db prot_database
    fi

    humann \\
        --threads $task.cpus \\
        --input $input \\
        $options.args \\
        --taxonomic-profile $profile \\
        --nucleotide-database nuc_database \\
        --protein-database prot_database \\
        --output out \\
        --output-basename $prefix \\
        --o-log ${prefix}.log

    gzip -n out/*.tsv
    cp out/${prefix}_genefamilies.tsv.gz ${prefix}_genefamilies.tsv.gz
    cp out/${prefix}_pathabundance.tsv.gz ${prefix}_pathabundance.tsv.gz
    cp out/${prefix}_pathcoverage.tsv.gz ${prefix}_pathcoverage.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( humann --version 2>&1 | sed 's/humann v//' )
    END_VERSIONS
    """
}
