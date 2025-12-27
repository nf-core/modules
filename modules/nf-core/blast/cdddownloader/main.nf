process BLAST_CDDDOWNLOADER {
    tag "$db_prefix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    input:
    val db_prefix
    path download_path

    output:
    path "${download_path}/cdd_databases/", emit: db
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
    # Check if cdd_databases directory exists in download_path
    if [ -d "${download_path}/cdd_databases" ]; then
        echo "Found existing cdd_databases directory at ${download_path}"
        cd ${download_path}/cdd_databases

        # Check if db_prefix directory exists
        if [ -d "${db_prefix}" ]; then
            file_count=\$(find ${db_prefix} -type f | wc -l)
            echo "WARNING: The directory ${download_path}/cdd_databases/${db_prefix} already exists and contains \${file_count} files."
            echo "If you want to download the database ${db_prefix} again, remove the current ${db_prefix} directory first."
        else
            echo "Creating ${db_prefix} directory and downloading database..."
            mkdir ${db_prefix}
            echo "Downloading ${db_prefix} database into ${db_prefix} dir"
            wget https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/${db_name}.tar.gz
            tar -xzf ${db_name}.tar.gz -C ./${db_prefix}
            rm -f ${db_name}.tar.gz
            echo "Database ${db_prefix} downloaded successfully"
        fi

        # Check if data directory exists
        if [ -d "data" ]; then
            data_file_count=\$(find data -type f | wc -l)
            echo "The directory ${download_path}/cdd_databases/data already exists and contains \${data_file_count} files."
            echo "Skipping metadata files downloading"
            echo "If you want to download the metadata files again, remove the current data directory first."
        else
            echo "Creating data directory and downloading metadata files..."
            mkdir data
            echo "Downloading metadata files"
            wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O ./data/cddid.tbl.gz && gzip -d ./data/cddid.tbl.gz
            wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt -O ./data/cdtrack.txt
            wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links -O ./data/family_superfamily_links
            wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz -O ./data/cddannot.dat.gz && gzip -d ./data/cddannot.dat.gz
            wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz -O ./data/cddannot_generic.dat.gz && gzip -d ./data/cddannot_generic.dat.gz
            wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt -O ./data/bitscore_specific.txt
            echo "Metadata files downloaded successfully"
        fi
    else
        echo "Creating new cdd_databases directory structure at ${download_path}"
        mkdir -p ${download_path}/cdd_databases/${db_prefix}
        mkdir -p ${download_path}/cdd_databases/data
        cd ${download_path}/cdd_databases

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
    fi

    echo "Finish"

    """

    stub:
    """
    mkdir output/cdd_databases/
    """
}
