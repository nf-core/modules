


process TCOFFEE_ALNCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::t-coffee=13.46.0.919e8c6b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t-coffee:13.45.0.4846264--hc57179f_5':
        'biocontainers/t-coffee:13.45.0.4846264--hc57179f_5'}"

    input:
    tuple val(meta), path(msa), path(ref_msa)

    output:
    tuple val(meta), path("*.scores"), emit: scores
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // First check if the flag is given at all and if not put the argument needed, then if given and contains the argument leave it as it is
    // otherwise add at the beginning the necessary flag to the given arg by the user
    def args = task.ext.args ? ( task.ext.args.contains('-compare_mode tc') ? task.ext.args  : ('-compare_mode tc ' + task.ext.args)) : '-compare_mode tc'
    def header = meta.keySet().join(",")
    def values = meta.values().join(",")

    """

    t_coffee -other_pg aln_compare \
        -al1 ${ref_msa} \
        -al2 ${msa} \
        ${args} \
        | grep -v "seq1" | grep -v '*' | \
        awk '{ print \$4}' ORS="\t" \
        >> "scores.txt"

    # Add metadata info to output file
    echo "${header},sp,tc,column" > "${msa.baseName}.scores"

    # Add values
    scores=\$(awk '{sub(/[[:space:]]+\$/, "")} 1' scores.txt | tr -s '[:blank:]' ',')
    echo "${values},\$scores" >> "${msa.baseName}.scores"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t_coffee: \$( t_coffee -version | awk -F'ersion_' '{print \$2}' | cut -d ' ' -f 1 )
    END_VERSIONS
    """
}
