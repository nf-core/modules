include { PREPROCESSGENEANNOTATION              } from '../../../modules/nf-core/aberrantexpression/counting/preprocessgeneannotation'
include { DEFINE_DATASET_FROM_ANNO              } from '../../../modules/nf-core/aberrantsplicing/precomp/definedatasetfromanno'
include { COUNTRNA_INIT                         } from '../../../modules/nf-core/aberrantsplicing/precomp/countrnainit'

include { COUNTRNA_SPLITREADS_SAMPLEWISE        } from '../../../modules/nf-core/aberrantsplicing/counting/countrnasplitreads'
include { COUNTRNA_SPLITREADS_MERGE             } from '../../../modules/nf-core/aberrantsplicing/counting/countrnasplitreadsmerge'

include { COUNTRNA_NONSPLITREADS_SAMPLEWISE     } from '../../../modules/nf-core/aberrantsplicing/counting/countrnanonsplitreads'
include { COUNTRNA_NONSPLITREADS_MERGE          } from '../../../modules/nf-core/aberrantsplicing/counting/countrnanonsplitreadsmerge'

include { COUNTRNA_COLLECT                      } from '../../../modules/nf-core/aberrantsplicing/counting/countrnacollect'
include { PSI_VALUE_CALC                        } from '../../../modules/nf-core/aberrantsplicing/counting/psivaluecalc'

workflow ABERRANT_SPLICING {
    take:
        params

    main:
        ch_versions = Channel.empty()

        PREPROCESSGENEANNOTATION(params.gtfData).result
            .map {it -> ["annotation": it[0], "txtDb": it[1], "geneMap": it[2], "countRanges": it[3]]}
            .set {preprocess_data}

        fMapping = params.procAnnotation
            .filter (it -> it.RNA_BAM_FILE != "" && it.DNA_VCF_FILE != "")
            .multiMap {
                it ->
                    RNA_BAM_FILE: tuple(it.RNA_ID, it.RNA_BAM_FILE)
                    DNA_VCF_FILE: tuple(it.RNA_ID, it.DNA_VCF_FILE)
                    SPLICE_COUNTS_DIR: tuple(it.RNA_ID, it.SPLICE_COUNTS_DIR)
            }
        DEFINE_DATASET_FROM_ANNO (
            params.aberrantsplicing.groups,
            fMapping.RNA_BAM_FILE.collect(),
            params.annotation
        ).result
            .map {it -> ["group": it[0], "tsv": it[1]]}
            .set {groupData}

        COUNTRNA_INIT (groupData).result
            .map {it -> ["group": it[0], "fdsObj": it[1], "output": it[2]]}
            .set {initData}

        sampleIDs = params.procAnnotation.combine(initData)
            .filter (it -> it[0].RNA_BAM_FILE != "")
            .filter (it -> it[0].DROP_GROUP.contains(it[1].group))
            .map {it -> ["group": it[1].group, "sampleID": it[0].RNA_ID, "output": it[1].output]}
        COUNTRNA_SPLITREADS_SAMPLEWISE(
                sampleIDs,
                params.genomeAssembly,
                params.aberrantsplicing.recount,
                params.aberrantsplicing.keepNonStandardChrs
            )
            .groupTuple()
            .map {it -> ["group": it[0], "output": it[1]]}
            .set {mergeIDs_data}

        // mergeIDs_data.view()
        COUNTRNA_SPLITREADS_MERGE(
                mergeIDs_data,
                params.aberrantsplicing.recount,
                params.aberrantsplicing.minExpressionInOneSample
            )
            .map {it -> ["groupData": it[0], "output": it[1], \
                "grSplit": it[2], "grNonSplit": it[3], "spliceSites": it[4]]}
            .set {splitmerge_data}

        splitmerge_data = splitmerge_data.combine(sampleIDs)
            .filter (it -> it[0].groupData.group == it[1].group)
            .map {it -> ["group": it[0].groupData.group, "sampleID": it[1].sampleID, \
                "output": it[0].output, "grSplit": it[0].grSplit, \
                "grNonSplit": it[0].grNonSplit, "spliceSites": it[0].spliceSites]}
        COUNTRNA_NONSPLITREADS_SAMPLEWISE(
            splitmerge_data,
            params.aberrantsplicing.recount,
            params.aberrantsplicing.longRead
        )
            .map {it -> [it[0].group, it[0].grNonSplit, it[0].grSplit, it[0].spliceSites, it[1]]}
            .groupTuple()
            .map {it -> ["group": it[0], "grNonSplit": it[1][0], \
                "grSplit": it[2][0], "spliceSites": it[3][0], "output": it[4]]}
            .set {nonSplit_data}

        COUNTRNA_NONSPLITREADS_MERGE(
                nonSplit_data,
                params.aberrantsplicing.recount,
                params.aberrantsplicing.longRead
            )
            .map {it -> ["group": it[0].group, "output": it[1], \
                "grSplit": it[0].grSplit, "grNonSplit": it[0].grNonSplit, \
                "spliceSites": it[0].spliceSites]}
            .set {nonSplitMerge_data}


        // nonSplitMerge_data.view()
        // COUNTRNA_COLLECT(nonSplitMerge_data)


        // psiCalcGroup = groups.
        //     map {it -> tuple(it.group, it.fdsObj)}
        //  v(psiCalcGroup, COUNTRNA_COLLECT.out.counting_done.collect())
}
