process BLAST_CDDDOWNLOADER {
    tag "$db_prefix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    input:
    val db_prefix

    output:
    path "cdd_databases/", emit: db
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -1 | cut -d ' ' -f 3"), topic: versions, emit: versions_wget
    tuple val("${task.process}"), val('untar'), eval("tar --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1"), topic: versions, emit: versions_tar

    when:
    task.ext.when == null || task.ext.when

    script:
    def db_name = 'Cdd_NCBI_LE'
    if ( "$db_prefix" ==~ /^Cdd$/ ) {
        db_name = 'Cdd_LE'
    } else if ( "$db_prefix" ==~ /^Cog$/ ) {
        db_name = 'Cog_LE'
    } else if ( "$db_prefix" ==~ /^Kog$/ ) {
        db_name = 'Kog_LE'
    } else if ( "$db_prefix" ==~ /^Pfam$/ ) {
        db_name = 'Pfam_LE'
    } else if ( "$db_prefix" ==~ /^Prk$/ ) {
        db_name = 'Prk_LE'
    } else if ( "$db_prefix" ==~ /^Smart$/ ) {
        db_name = 'Smart_LE'
    } else if ( "$db_prefix" ==~ /^Tigr$/ ) {
        db_name = 'Tigr_LE'
    } else {
        log.warn("Unknown CDD databse name (${db_prefix}): selecting Cdd_NCBI default of downloading")
        db_prefix = 'Cdd_NCBI'
    }

    """
    mkdir -p cdd_databases/${db_prefix}
    cd cdd_databases/
    mkdir data

    echo "Downloading ${db_prefix} database into ${db_prefix} dir"

    wget https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/${db_name}.tar.gz
    tar -xzf ${db_name}.tar.gz -C ./${db_prefix}
    rm -f ${db_name}.tar.gz

    echo "Downloading metadata files"

    wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O ./data/cddid.tbl.gz && gzip -d ./data/cddid.tbl.gz
    wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt -O ./data/cdtrack.txt
    wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links -O ./data/family_superfamily_links
    wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz -O ./data/cddannot.dat.gz && gzip -d ./data/cddannot.dat.gz
    wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz -O ./data/cddannot_generic.dat.gz && gzip -d ./data/cddannot_generic.dat.gz
    wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt -O ./data/bitscore_specific.txt

    echo "Finish"

    """

    stub:
    """
    mkdir cdd_databases/
    """
}
