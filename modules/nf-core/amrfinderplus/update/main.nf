// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00cad0bb4c7cb85b0906f357e5197bcabb32569f9b48ff049f081b5bafd04aca/data'
        : 'community.wave.seqera.io/library/ncbi-amrfinderplus:4.2.7--01f7d2a8fad0d3aa'}"

    input:
    val db_version

    output:
    path "amrfinderdb.tar.gz", emit: db
    tuple val("${task.process}"), val('amrfinder'), eval('amrfinder --version'), emit: versions_amrfinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!db_version) {
        """
        amrfinder_update -d amrfinderdb
        tar czvf amrfinderdb.tar.gz -C amrfinderdb/\$(readlink amrfinderdb/latest) ./
        """
    } else {
        """
        # Database releases are grouped by the database format of the AMRFinderPlus version (e.g. 4.2 for 4.2.7)
        db_format=\$(amrfinder --version | cut -d '.' -f 1,2)
        db_url="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/\${db_format}/${db_version}/"

        mkdir -p amrfinderdb/${db_version}
        # wget fails on a missing version; the check below then reports it clearly
        wget -q -r -l 1 -np -nd -nH -e robots=off -R 'index.html*' -P amrfinderdb/${db_version} "\${db_url}" || true
        if [ ! -s amrfinderdb/${db_version}/version.txt ]; then
            echo "ERROR: AMRFinderPlus database version '${db_version}' not found at \${db_url}" >&2
            exit 1
        fi

        amrfinder_index amrfinderdb/${db_version}
        tar czvf amrfinderdb.tar.gz -C amrfinderdb/${db_version} ./
        """
    }

    stub:
    """
    touch amrfinderdb.tar
    gzip amrfinderdb.tar
    """
}
