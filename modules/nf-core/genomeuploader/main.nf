process GENOMEUPLOADER {
    tag "$meta.id"
    label 'process_single'

    secret secrets.ENA_WEBIN ? "ENA_WEBIN" : ""
    secret secrets.ENA_WEBIN_PASSWORD ? "ENA_WEBIN_PASSWORD" : ""

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genome-uploader:2.5.0--pyhdfd78af_0':
        'biocontainers/genome-uploader:2.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(metadata_tsv)
    path(fastas, stageAs: "genomes/*")

    output:
    tuple val(meta), path("upload_output/*")                                        , emit: upload_output_dir
    tuple val(meta), path("upload_output/MAG_upload/submission.xml")                , emit: submission
    tuple val(meta), path("upload_output/MAG_upload/registered_MAGs_${prefix}.tsv") , emit: registered_mags
    tuple val(meta), path("upload_output/MAG_upload/genome_samples.xml")            , emit: genome_samples
    tuple val(meta), path("upload_output/MAG_upload/manifests_${prefix}/*.manifest"), emit: manifests
    path "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    genome_upload \\
        $args \\
        --upload_study "${meta.study_accession}" \\
        --centre_name "${meta.center_name}" \\
        --genome_info ${metadata_tsv} \\
        --out upload_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome_upload: \$( genome_upload --version | sed 's/genome_uploader //' )
        ena-webin-cli: \$( ena-webin-cli -version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p upload_output/MAG_upload/manifests_${prefix}
    touch upload_output/MAG_upload/submission.xml
    touch upload_output/MAG_upload/registered_MAGs_${prefix}.tsv
    touch upload_output/MAG_upload/genome_samples.xml
    touch upload_output/MAG_upload/manifests_${prefix}/test_mag.manifest

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome_upload: \$( genome_upload --version | sed 's/genome_uploader //' )
        ena-webin-cli: \$( ena-webin-cli -version )
    END_VERSIONS
    """
}
