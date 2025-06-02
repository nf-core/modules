process ULTRA_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0':
        'biocontainers/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.db")    , emit: database
    tuple val(meta), path("*.pickle"), emit: pickle
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    uLTRA \\
        index \\
        ${args} \\
        ${fasta} \\
        ${gtf} \\
        ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        gffutils: \$(python -c "import gffutils; print(gffutils.__version__)")
        sqlite: \$(python -c "import sqlite3; print(sqlite3.sqlite_version)")
    END_VERSIONS
    """

    stub:
    """
    touch database.db
    touch all_splice_pairs_annotations.pickle
    touch all_splice_sites_annotations.pickle
    touch chr_to_id.pickle
    touch exon_choordinates_to_id.pickle
    touch flank_choordinates.pickle
    touch gene_to_small_segments.pickle
    touch id_to_chr.pickle
    touch max_intron_chr.pickle
    touch parts_to_segments.pickle
    touch ref_exon_sequences.pickle
    touch ref_flank_sequences.pickle
    touch ref_part_sequences.pickle
    touch ref_segment_sequences.pickle
    touch refs_id_lengths.pickle
    touch refs_lengths.pickle
    touch segment_id_to_choordinates.pickle
    touch segment_to_gene.pickle
    touch segment_to_ref.pickle
    touch splices_to_transcripts.pickle
    touch transcripts_to_splices.pickle

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        gffutils: \$(python -c "import gffutils; print(gffutils.__version__)")
        sqlite: \$(python -c "import sqlite3; print(sqlite3.sqlite_version)")
    END_VERSIONS
    """
}
