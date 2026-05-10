process HIFIADAPTERFILT_DOWNLOADDB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiadapterfilt:3.0.0--hdfd78af_0':
        'quay.io/biocontainers/hifiadapterfilt:3.0.0--hdfd78af_0' }"

    output:
    path("DB")                                                                                                             , emit: db
    tuple val("${task.process}"), val('hifiadapterfilt'), eval("hifiadapterfilt.sh --version 2>&1 | head -1"), topic: versions, emit: versions_hifiadapterfilt

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p DB

    # Download all pre-built BLAST database files for PacBio adapter sequences
    # (NGB00972.1 Blunt Adapter, NGB00973.1 C2 Primer) from the HiFiAdapterFilt repo
    for file in pacbio_vectors_db pacbio_vectors_db.ndb pacbio_vectors_db.nhr \\
                pacbio_vectors_db.nin pacbio_vectors_db.nog pacbio_vectors_db.nos \\
                pacbio_vectors_db.not pacbio_vectors_db.nsq pacbio_vectors_db.ntf \\
                pacbio_vectors_db.nto; do
        wget -q -O DB/\${file} \\
            https://raw.githubusercontent.com/sheinasim-USDA/HiFiAdapterFilt/master/DB/\${file}
    done
    """

    stub:
    """
    mkdir -p DB
    touch DB/pacbio_vectors_db
    touch DB/pacbio_vectors_db.ndb
    touch DB/pacbio_vectors_db.nhr
    touch DB/pacbio_vectors_db.nin
    touch DB/pacbio_vectors_db.nog
    touch DB/pacbio_vectors_db.nos
    touch DB/pacbio_vectors_db.not
    touch DB/pacbio_vectors_db.nsq
    touch DB/pacbio_vectors_db.ntf
    touch DB/pacbio_vectors_db.nto
    """
}
