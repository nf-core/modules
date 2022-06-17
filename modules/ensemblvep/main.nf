process ENSEMBLVEP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ensembl-vep=105.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:105.0--pl5321h4a94de4_1' :
        'quay.io/biocontainers/ensembl-vep:105.0--pl5321h4a94de4_1' }"

    input:
    tuple val(meta), path(vcf)
    val   genome
    val   species
    val   cache_version
    path  cache
    path  fasta
    path  extra_files

    output:
    tuple val(meta), path("*.ann.vcf")  , optional:true, emit: vcf
    tuple val(meta), path("*.ann.tab")  , optional:true, emit: tab
    tuple val(meta), path("*.ann.json") , optional:true, emit: json
    path "*.summary.html"               , emit: report
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""

    """
    vep \\
        -i $vcf \\
        -o ${prefix}.ann.${file_extension} \\
        $args \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
