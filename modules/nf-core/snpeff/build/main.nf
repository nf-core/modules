process SNPEFF_BUILD {
    tag "$db_name"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30669e5208952f30d59d0d559928772f082830d01a140a853fff13a2283a17b0/data'
        : 'community.wave.seqera.io/library/snpeff:5.4.0a--eaf6ce30125b2b17'}"

    input:
    tuple val(meta_ref), path(fasta), path(annotation), path(cds), path(protein)    // cds and protein are optional
    // TODO: optionally use , arity: '0..*' or typed inputs? See https://github.com/nextflow-io/nextflow/issues/5111 & https://github.com/nextflow-io/nextflow/issues/1694
    // TODO: optionally use stageAs to change name immediately without relying on the bash script to do it
    val annotation_format           // 'gff', 'gtf' or empty (falls back to detection in filename)
    tuple val(meta_bed), path(bed)  // BED file for mitochondrial/apicoplast detection
    path snpeff_config_template     // TODO: instead of requiring users to supply this file, it could also be read from the snpEff install directory
    val db_name

    output:
    // TODO: `tuple val(meta_ref), path("snpeff_db"), emit: db` could serve as an alternative approach => outputs self-contained directory with snpEff config and database, but requires the use of a relative -dataDir option in the snpEff annotate command
    tuple val(meta_ref), path("snpeff_db/data"),            emit: db
    tuple val(meta_ref), path("snpeff_db/snpEff.config"),   emit: config
    tuple val("${task.process}"), val('snpeff'), eval("snpEff -version 2>&1 | cut -f 2 -d '\t'"), topic: versions, emit: versions_snpeff

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // set db_name to meta.id of reference when not supplied
    db_name = db_name ?: meta_ref.id

    // extract gff/gtf format from filename if not provided
    if (!annotation_format) {
        def anno_name = annotation.name.toLowerCase()
        if (anno_name.endsWith('.gtf') || anno_name.endsWith('.gtf.gz')) {
            annotation_format = 'gtf'
        } else if (anno_name.endsWith('.gff') || anno_name.endsWith('.gff.gz') ||
                anno_name.endsWith('.gff3') || anno_name.endsWith('.gff3.gz')) {
            annotation_format = 'gff'
        } else {
            error "Could not determine annotation format from filename: ${annotation.name}. " +
                "Please provide annotation_format parameter ('gtf' or 'gff')."
        }
    }
    if (annotation_format != 'gtf' && annotation_format != 'gff' && annotation_format != 'gff3' ) {
        error "Invalid annotation_format: '${annotation_format}'. Must be 'gtf' or 'gff(3)'."
    }
    def annotation_file = (annotation_format == 'gtf') ? 'genes.gtf' : 'genes.gff'
    def annotation_arg = (annotation_format == 'gtf') ? '-gtf22' : '-gff3'

    // add cli arguments to skip checks for cds/protein files when they are not provided
    def no_check_cds_arg = cds ? '' : '-noCheckCds'
    def no_check_protein_arg = protein ? '' : '-noCheckProtein'

    """
    # Create the directory structure snpEff expects
    mkdir -p snpeff_db/data/${db_name}

    # Copy and rename files to match snpEff naming requirements, unzipping them if necessary
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c ${fasta} > snpeff_db/data/${db_name}/sequences.fa
    else
        cp ${fasta} snpeff_db/data/${db_name}/sequences.fa
    fi

    if [[ "${annotation}" == *.gz ]]; then
        gunzip -c ${annotation} > snpeff_db/data/${db_name}/${annotation_file}
    else
        cp ${annotation} snpeff_db/data/${db_name}/${annotation_file}
    fi

    # Only copy CDS and proteins files if provided
    # Note: quotes around variables during file existence check are critical,
    # otherwise the tests will default to true when the variable is an empty string
    # e.g., `if [ -f \$undeclared_var ]` = `if [ -f ]` = true
    if [ -f "${cds}" ]; then
        if [[ "${cds}" == *.gz ]]; then
            gunzip -c ${cds} > snpeff_db/data/${db_name}/cds.fa
        else
            cp ${cds} snpeff_db/data/${db_name}/cds.fa
        fi
    else
        echo "No CDS file provided, skipping CDS check."
    fi

    if [ -f "${protein}" ]; then
        if [[ "${protein}" == *.gz ]]; then
            gunzip -c ${protein} > snpeff_db/data/${db_name}/protein.fa
        else
            cp "${protein}" "snpeff_db/data/${db_name}/protein.fa"
        fi
    else
        echo "No protein file provided, skipping protein check."
    fi

    # Create snpEff config file, starting with the template file
    # Note: dataDir in config file will be ignored since it is overridden by the CLI option
    cp ${snpeff_config_template} snpeff_db/snpEff.config

    # Append custom genome configuration to config
    cat >> snpeff_db/snpEff.config << EOF
# ${db_name} genome configuration
${db_name}.genome : ${db_name}
EOF

    # Build the database - dataDir is relative to the location of snpEff.config when provided as a relative path
    snpEff build \\
        -v \\
        -c snpeff_db/snpEff.config \\
        -dataDir ./data/ \\
        ${annotation_arg} \\
        ${no_check_cds_arg} \\
        ${no_check_protein_arg} \\
        ${args} \\
        ${db_name}
    """

    stub:
    db_name = db_name ?: meta_ref.id
    // extract gff/gtf format from filename if not provided
    if (!annotation_format) {
        def anno_name = annotation.name.toLowerCase()
        if (anno_name.endsWith('.gtf') || anno_name.endsWith('.gtf.gz')) {
            annotation_format = 'gtf'
        } else if (anno_name.endsWith('.gff') || anno_name.endsWith('.gff.gz') ||
                anno_name.endsWith('.gff3') || anno_name.endsWith('.gff3.gz')) {
            annotation_format = 'gff'
        } else {
            annotation_format = 'gtf'
        }
    }
    def annotation_file = (annotation_format == 'gtf') ? 'genes.gtf' : 'genes.gff'
    """
    # Create the expected directory structure
    mkdir -p snpeff_db/data/${db_name}

    # Create empty files matching the real process outputs
    touch snpeff_db/data/${db_name}/sequences.fa
    touch snpeff_db/data/${db_name}/${annotation_file}
    touch snpeff_db/data/${db_name}/cds.fa
    touch snpeff_db/data/${db_name}/protein.fa

    # Create config file in correct location
    touch snpeff_db/snpEff.config
    """
}
