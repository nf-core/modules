/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS TO EXTRACT REFERENCES FILES OR VALUES FROM THE REFERENCES YAML OR PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_references_file(references, param, attribute, basepath) {
    return references
        .map { meta, _readme ->
            if (param || meta[attribute]) {
                [meta.subMap(['id']), file(param ?: meta[attribute].replace('${params.igenomes_base}', basepath), checkIfExists: true)]
            }
            else {
                null
            }
        }
        .collect()
}

def get_references_value(references, param, attribute) {
    return references
        .map { meta, _readme ->
            if (param || meta[attribute]) {
                [meta.subMap(['id']), param ?: meta[attribute]]
            }
            else {
                null
            }
        }
        .collect()
}
