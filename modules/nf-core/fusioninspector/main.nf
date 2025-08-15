process FUSIONINSPECTOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/139b94a1f10c6e23a8c27eaed1e5a689db978a513d0ee155e74d35f0970814fe/data' :
        'community.wave.seqera.io/library/fusion-inspector_igv-reports_perl-json-xs_pysam_pruned:c6147971d107ab11'}"

    input:
    tuple val(meta), path(reads), path(fusion_list)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*FusionInspector.fusions.tsv"), emit: tsv         , optional:true
    tuple val(meta), path("fi_workdir/*.gtf")            , emit: out_gtf     , optional:true
    tuple val(meta), path("*FusionInspector.log")        , emit: log         , optional:true
    tuple val(meta), path("*html")                       , emit: html        , optional:true
    tuple val(meta), path("*abridged.tsv")               , emit: abridged_tsv, optional:true
    tuple val(meta), path("IGV_inputs")                  , emit: igv_inputs  , optional:true
    tuple val(meta), path("fi_workdir")                  , emit: fi_workdir  , optional:true
    tuple val(meta), path("chckpts_dir")                 , emit: chckpts_dir , optional:true
    path  "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta  = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    """
    FusionInspector \\
        --fusions $fusion_list \\
        --genome_lib ${reference} \\
        $fasta \\
        --CPU ${task.cpus} \\
        -O . \\
        --out_prefix $prefix \\
        --vis $args $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FusionInspector: \$(FusionInspector --version 2>&1 | grep -i 'version' | sed -e 's/FusionInspector version: //' -e 's/[[:space:]]//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch FusionInspector.log
    touch ${prefix}.FusionInspector.fusions.abridged.tsv
    touch ${prefix}.FusionInspector.fusions.tsv
    touch ${prefix}.fusion_inspector_web.html
    mkdir -p chckpts_dir
    mkdir -p fi_workdir
    touch fi_workdir/${prefix}.gtf
    mkdir -p IGV_inputs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FusionInspector: \$(FusionInspector --version 2>&1 | grep -i 'version' | sed -e 's/FusionInspector version: //' -e 's/[[:space:]]//g')
    END_VERSIONS
    """
}
