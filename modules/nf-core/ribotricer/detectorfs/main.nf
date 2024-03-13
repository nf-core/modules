process RIBOTRICER_DETECTORFS {
    tag '$bam'
    label 'process_single'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotricer:1.3.3--pyhdfd78af_0':
        'biocontainers/ribotricer:1.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(candidate_orfs)

    output:
    tuple val(meta), path('*_protocol.txt')             , emit: protocol, optional: true
    tuple val(meta), path('*_bam_summary.txt')          , emit: bam_summary
    tuple val(meta), path('*_read_length_dist.pdf')     , emit: read_length_dist
    tuple val(meta), path('*_metagene_profiles_5p.tsv') , emit: metagene_profile_5p
    tuple val(meta), path('*_metagene_profiles_3p.tsv') , emit: metagene_profile_3p
    tuple val(meta), path('*_metagene_plots.pdf')       , emit: metagene_plots
    tuple val(meta), path('*_psite_offsets.txt')        , emit: psite_offsets, optional: true
    tuple val(meta), path('*_pos.wig')                  , emit: pos_wig
    tuple val(meta), path('*_neg.wig')                  , emit: neg_wig
    tuple val(meta), path('*_translating_ORFs.tsv')     , emit: orfs
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def strandedness_cmd = ''
    if (meta.strandedness == 'forward') {
        strandedness_cmd = '--stranded yes'
    } else if (meta.strandedness == 'reverse') {
        strandedness_cmd = '--stranded revers'
    } else if (mea.strandedness == 'unstranded') {
        strandedness_cmd = '--stranded no'
    }
    """
    ribotricer detect-orfs \\
        --bam $bam \\
        --ribotricer_index $candidate_orfs \\
        --prefix $prefix \\
        $strandedness_cmd \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotricer: \$(ribotricer --version | grep ribotricer |& sed '1!d ; s/ribotricer, version //')
    END_VERSIONS
    """
}
