process INTERPROSCAN {
    tag "$meta.id"
    label 'process_long'

    conda "bioconda::interproscan=5.59_91.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)
    val(out_ext)

    output:
    tuple val(meta), path('*.tsv'), optional: true, emit: tsv
    tuple val(meta), path('*.xml'), optional: true, emit: xml
    tuple val(meta), path('*.gff3'), optional: true, emit: gff3
    tuple val(meta), path('*.json'), optional: true, emit: json
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def appl = "-appl TIGRFAM,FunFam,SFLD,PANTHER,Gene3D,Hamap,ProSiteProfiles,Coils,SMART,CDD,PRINTS,PIRSR,ProSitePatterns,AntiFam,Pfam,MobiDBLite"
    if ( args.contains("-appl") ) {
        appl = ""
    }
    switch ( out_ext ) {
        case "tsv": break
        case "xml": break
        case "gff3": break
        case "json": break
        default:
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting tsv as fallback");
            break
    }

    //  -dp (disable precalculation) is on so no online dependency
    """
    interproscan.sh \\
        -cpu $task.cpus \\
        -i $fasta \\
        -f ${out_ext} \\
        -dp \\
        ${appl} \\
        ${args} \\
        -o ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n 1 | sed 's/^.*InterProScan version//' | sed 's/\\s*InterProScan.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    switch ( out_ext ) {
        case "tsv": break
        case "xml": break
        case "gff3": break
        case "json": break
        default:
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting tsv as fallback");
            break
    }

    """
    touch ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n 1 | sed 's/^.*InterProScan version//' | sed 's/\\s*InterProScan.*//')
    END_VERSIONS
    """
}
