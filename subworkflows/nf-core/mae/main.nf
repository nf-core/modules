include { CREATESNVS        } from '../../../modules/nf-core/mae/counting/createsvns/main.nf'
include { ALLELIC_COUNTS    } from '../../../modules/nf-core/mae/counting/alleliccounts/main.nf'

include { DESEQ             } from '../../../modules/nf-core/mae/test/deseq/main.nf'
include { GENE_NAME_MAPPING } from '../../../modules/nf-core/mae/test/genemapping/main.nf'
include { MAE_RESULTS       } from '../../../modules/nf-core/mae/test/maeresults/main.nf'

include { DESEQ_QC          } from '../../../modules/nf-core/mae/dnarnamatch/deseq/main.nf'
include { MATRIX_DNA_RNA    } from '../../../modules/nf-core/mae/dnarnamatch/matrixdnarna/main.nf'

workflow MAE {
    take:
        params

    main:
        ch_versions = Channel.empty()

        preprocess_data = params.procAnnotation
            .filter {it -> it.DROP_GROUP.contains("mae")}
            .map    {it -> ["vcf": it.DNA_ID, "rnaID": it.RNA_ID, "vcf_file": it.DNA_VCF_FILE, "bam_file": it.RNA_BAM_FILE]}
            .combine (Channel.of("nQC", "QC"))
            .map     {it -> { if (it[1] == "nQC") { it[0] }
                        else {
                            a = [:]
                            a.putAll(it[0])
                            a["vcf"] = "QC"
                            a["vcf_file"] = params.mae.qc_vcf
                            a
                        }}}
        CREATESNVS(preprocess_data, params.mae).result
            .map {it -> ["vcf": it[0], "rnaID": it[1], "BAM": it[2], "SNV": it[3]]}
            .map {it -> {if (it.BAM.contains("ncbi")) {it["fasta"] = params.mae.chr1_ncbi}
                    else {
                        it['fasta'] = params.mae.chr1
                    }
                    it
                }}
            .set {create_snv}
        ALLELIC_COUNTS(create_snv, params.mae).result
            .map {it -> ["vcf": it[0], "rnaID": it[1], "counts": it[2]]}
            .set {allelic_data}

        allelic_data_qc = allelic_data.filter {it -> it.vcf == "QC"}
        allelic_data = allelic_data.filter {it -> it.vcf != "QC"}
        DESEQ(allelic_data, params)
            .map {it -> [0, it[2]]}
            .groupTuple()
            .set {deseq_data}

        DESEQ_QC(allelic_data_qc)
            .map {it -> [it[2]]}
            .collect()
            .set {deseq_qc_data}

        gtf_data = Channel.from(params.gtfData)
            .map {it -> ["annotation": it.version, "gtf": it.path]}
        GENE_NAME_MAPPING(gtf_data).result
            .map {it -> ["annotation": it[0], "tsv": it[1]]}
            .set {genemap_data}

        results_data = genemap_data.combine(deseq_data)
            .map {it -> ["annotation": it[0].annotation, "genemap": it[0].tsv, "res": it[2]]}
        MAE_RESULTS(results_data, params.mae)

        rnaIDs = params.procAnnotation
            .filter(it -> it.DROP_GROUP.contains("mae"))
            .map   (it -> it.RNA_ID)
            .collect()

        MATRIX_DNA_RNA(deseq_qc_data, rnaIDs, params)
}
