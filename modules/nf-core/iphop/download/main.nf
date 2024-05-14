process IPHOP_DOWNLOAD {
    label 'process_single'

    conda "bioconda::iphop=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iphop:1.3.2--pyhdfd78af_0':
        'biocontainers/iphop:1.3.2--pyhdfd78af_0' }"

    output:
    path "iphop_db/"        , emit: iphop_db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p download_dir
    mkdir -p iphop_db

    iphop \\
        download \\
        --db_dir download_dir \\
        --no_prompt \\
        $args

    rm download_dir/*.tar.*
    mv download_dir/*/* iphop_db

    iphop \\
        download \\
        --db_dir iphop_db \\
        --no_prompt \\
        --full_verify \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir -p iphop_db/db/
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.fna
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.ndb
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.nhr
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.nin
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.not
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.nsq
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.ntf
    touch iphop_db/db/All_CRISPR_spacers_nr_clean.nto
    touch iphop_db/db/GTDBtkr202_and_newrepr_s2_mat.pkl
    mkdir -p iphop_db/db/Host_Genomes
    touch iphop_db/db/Host_Genomes/Host_Genomes.ndb
    touch iphop_db/db/Host_Genomes/Host_Genomes.nhr
    touch iphop_db/db/Host_Genomes/Host_Genomes.nin
    touch iphop_db/db/Host_Genomes/Host_Genomes.not
    touch iphop_db/db/Host_Genomes/Host_Genomes.nsq
    touch iphop_db/db/Host_Genomes/Host_Genomes.ntf
    touch iphop_db/db/Host_Genomes/Host_Genomes.nto
    mkdir -p iphop_db/db/php_db
    mkdir -p iphop_db/db/rafah_data
    touch iphop_db/db/rafah_data/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3f
    touch iphop_db/db/rafah_data/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3i
    touch iphop_db/db/rafah_data/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3m
    touch iphop_db/db/rafah_data/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3p
    touch iphop_db/db/rafah_data/HP_Ranger_Model_3_Valid_Cols.txt
    touch iphop_db/db/rafah_data/MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData
    touch iphop_db/db/rafah_data/RaFAH_ref_cds.count.tsv
    touch iphop_db/db/rafah_data/RaFAH_ref_cds.dmnd
    mkdir -p iphop_db/db/rewish_models
    touch iphop_db/db/rewish_models/Batch_1.pkl
    mkdir -p iphop_db/db/wish_data/Decoy_db
    touch iphop_db/db/wish_data/Decoy_db/Decoy_phages.fna
    mkdir -p iphop_db/db_infos
    touch iphop_db/db_infos/All_CRISPR_array_size.tsv
    touch iphop_db/db_infos/All_CRISPR_spacers_nr_clean.metrics.csv
    touch iphop_db/db_infos/Host_Genomes.tsv
    touch iphop_db/db_infos/List_contigs_removed_blast.tsv
    touch iphop_db/db_infos/Translate_genus_to_full_taxo.tsv
    touch iphop_db/db_infos/Wish_negFits.csv
    touch iphop_db/db_infos/gtdbtk.ar122.decorated.tree
    touch iphop_db/db_infos/gtdbtk.bac120.decorated.tree
    touch iphop_db/md5checkfile.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' ))
    END_VERSIONS
    """
}
