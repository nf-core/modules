process FCS_FCSADAPTOR {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.5.0/fcs-adaptor.sif':
        'docker.io/ncbi/fcs-adaptor:0.5.0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.cleaned_sequences.fa.gz"), emit: cleaned_assembly, optional: true
    tuple val(meta), path("*.fcs_adaptor_report.txt") , emit: adaptor_report
    tuple val(meta), path("*.fcs_adaptor.log")        , emit: log
    tuple val(meta), path("*.pipeline_args.yaml")     , emit: pipeline_args
    tuple val(meta), path("*.skipped_trims.jsonl")    , emit: skipped_trims
    path "versions.yml"                               , emit: versions

    // Downstream handling of optional cleaned_assembly
    //
    // ch_cleaned_assembly  = FCS_FCSADAPTOR ( ch_input_assembly ).cleaned_assembly

    // ch_input_assembly
    // | join (ch_cleaned_assembly, by:0, remainder:true )
    // | map { meta, input, cleaned ->
    //     [ meta, cleaned ?: input ]
    // }
    // | set { ch_downstream_channel }

    // FCS_FCSGX ( ch_downstream_channel )

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSADAPTOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: '--prok' // --prok || --euk
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSADAPTOR_VERSION = '0.5.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    av_screen_x \\
        -o output/ \\
        $args \\
        $assembly

    # compress and/or rename files with prefix
    num_contamination_lines=\$(cat "output/fcs_adaptor_report.txt" | wc -l)
    if [[ \$num_contamination_lines -gt 1 ]]; then
        find output/cleaned_sequences/  -type f ! -name "*.gz" -exec gzip {} \\;
        cp output/cleaned_sequences/*     "${prefix}.cleaned_sequences.fa.gz"
    fi
    cp "output/fcs_adaptor_report.txt"    "${prefix}.fcs_adaptor_report.txt"
    cp "output/fcs_adaptor.log"           "${prefix}.fcs_adaptor.log"
    cp "output/pipeline_args.yaml"        "${prefix}.pipeline_args.yaml"
    cp "output/skipped_trims.jsonl"       "${prefix}.skipped_trims.jsonl"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FCS-adaptor: $FCSADAPTOR_VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSADAPTOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSADAPTOR_VERSION = '0.5.0'

    """
    touch ${prefix}.fcs_adaptor_report.txt
    touch ${prefix}.fcs_adaptor.log
    touch ${prefix}.pipeline_args.yaml
    touch ${prefix}.skipped_trims.jsonl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FCS-adaptor: $FCSADAPTOR_VERSION
    END_VERSIONS
    """
}
