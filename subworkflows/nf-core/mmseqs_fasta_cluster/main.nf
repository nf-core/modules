include { MMSEQS_CREATEDB  } from '../../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CLUSTER   } from '../../../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_LINCLUST  } from '../../../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV } from '../../../modules/nf-core/mmseqs/createtsv/main'

workflow MMSEQS_FASTA_CLUSTER {
    take:
    sequences       // tuple val(meta), path(fasta)
    clustering_tool // string: ["linclust", "cluster"]

    main:
    ch_versions       = Channel.empty()
    ch_clustering_tsv = Channel.empty()

    MMSEQS_CREATEDB( sequences )
    ch_versions = ch_versions.mix( MMSEQS_CREATEDB.out.versions )

    if (clustering_tool == 'cluster') {
        cluster_res = MMSEQS_CLUSTER( MMSEQS_CREATEDB.out.db )
        ch_versions = ch_versions.mix( MMSEQS_CLUSTER.out.versions )
    } else if (clustering_tool == 'linclust') {
        cluster_res = MMSEQS_LINCLUST( MMSEQS_CREATEDB.out.db )
        ch_versions = ch_versions.mix( MMSEQS_LINCLUST.out.versions )
    }

    // Join to ensure in sync in case of multiple sequence files
    ch_input_for_createtsv = MMSEQS_CREATEDB.out.db
        .join(cluster_res.db_cluster)
        .multiMap { meta, db, db_cluster ->
            db: [ meta, db ]
            db_cluster: [ meta, db_cluster ]
        }

    MMSEQS_CREATETSV( ch_input_for_createtsv.db_cluster, ch_input_for_createtsv.db, ch_input_for_createtsv.db )
    ch_versions = ch_versions.mix( MMSEQS_CREATETSV.out.versions )
    ch_clustering_tsv = MMSEQS_CREATETSV.out.tsv

    // Join to ensure in sync in case of multiple sequence files
    ch_clustering_output = sequences
        .join(MMSEQS_CREATETSV.out.tsv)
        .multiMap { meta, seqs, clusters ->
            seqs: [meta, seqs]
            clusters: [meta, clusters]
        }

    emit:
    versions = ch_versions
    seqs     = ch_clustering_output.seqs
    clusters = ch_clustering_output.clusters
}
