process EXOMISER {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::exomiser-cli=13.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/exomiser-cli:13.2.1--hdfd78af_0':
        'biocontainers/exomiser-cli:13.2.1--hdfd78af_0' }"

    input:
    tuple val(meta) , path(vcf), val(assembly), path(analysis_yml), path(job_yml) , path(output_yml), path(ped), path(sample_phenopacket_yml)
    tuple val(meta2), path(exomiser_cache, stageAs: "data/*")
    tuple val(meta2), val(exomiser_cache_version)

    output:
    tuple val(meta), path("*.vcf.gz")       , optional:true, emit: vcf
    tuple val(meta), path("*.html")         , optional:true, emit: html
    tuple val(meta), path("*.json")         , optional:true, emit: json
    tuple val(meta), path("*genes.tsv")     , optional:true, emit: genes_tsv
    tuple val(meta), path("*variants.tsv")  , optional:true, emit: variants_tsv

    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def EXOMISER_VERSION = 13.2.1
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def analysis_yml_cmd            = analysis_yml              ? "--analysis ${analysis_yml}"          : ""
    def assembly_cmd                = assembly                  ? "--assembly ${assembly}"              : ""
    def job_yml_cmd                 = job_yml                   ? "--job ${job_yml}"                    : ""
    def output_yml_cmd              = output_yml                ? "--output ${output_yml}"              : ""
    def ped_cmd                     = ped                       ? "--ped ${ped}"                        : ""
    def sample_phenopacket_yml_cmd  = sample_phenopacket_yml    ? "--sample ${sample_phenopacket_yml}"  : ""

    """
    ## generate application.properties to indicate where the reference data is and which version to use
    cat <<-END_PROPERTIES > application.properties
    exomiser.data-directory=./data
    exomiser.hg19.data-version=${exomiser_cache_version}
    exomoser.hg38.data-version=${exomiser_cache_version}
    exomiser.phenotype.data-version=${exomiser_cache_version}
    END_PROPERTIES

    exomiser-cli \\
        ${analysis_yml_cmd} \\
        ${assembly_cmd} \\
        ${job_yml_cmd } \\
        ${output_yml_cmd} \\
        --output-directory . \\
        --output-filename ${prefix} \\
        ${ped_cmd} \\
        ${sample_phenopacket_yml_cmd} \\
        --vcf ${vcf} \\
        --spring.config.location=./application.properties
        $args

    ## $args can contain
    ## --output-format
    ## --output-prefix (deprecated)
    ## --preset (mandatory)


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : ${EXOMISER_VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ## generate application.properties to indicate where the reference data is and which version to use
    cat <<-END_PROPERTIES > application.properties
    exomiser.data-directory=./data
    exomiser.hg19.data-version=${exomiser_cache_version}
    exomoser.hg38.data-version=${exomiser_cache_version}
    exomiser.phenotype.data-version=${exomiser_cache_version}
    END_PROPERTIES

    touch ${prefix}.vcf
    touch ${prefix}.html
    touch ${prefix}.json
    touch ${prefix}.genes.tsv
    touch ${prefix}.variants.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : ${EXOMISER_VERSION}
    END_VERSIONS
    """
}
