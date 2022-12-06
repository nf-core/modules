#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEAGCT } from '../../../../../modules/nf-core/custom/tabulartogseagct/main.nf'

// Example input from rnaseq workflow has an extra column for symbols we have
// to strip, and we want to test both tsv and csv. This function does the prep
// for both.

def prepareExample(sourceFilePath, targetFilePath, strip_columns = null) {

    sourceDelimiter = targetDelimiter = '\t'
    if (sourceFilePath.endsWith('.csv')){
        sourceDelimiter = ','
    }
    if (targetFilePath.endsWith('.csv')){
        targetDelimiter = ','
    }

    sourceFile = file(sourceFilePath, checkIfExists: true)
    targetFile = file(targetFilePath)

    sourceFile.withReader { source ->
        targetFile.withWriter { target ->
            String line
            while( line=source.readLine() ) {

                parts = line.tokenize(sourceDelimiter)
                if (strip_columns){
                    parts.removeAt(strip_columns)
                }
                line = parts.join(targetDelimiter)
                
                target << parts.join(targetDelimiter) + '\n'
            }
        }
    }
    println(targetFile.getClass())
    return targetFile
}

workflow test_custom_tabulartogseagct {
    
    infile = prepareExample(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], 'test.tsv', strip_columns = 1)
    input = [ [ id:'test' ], infile]    

    CUSTOM_TABULARTOGSEAGCT ( input )
}

workflow test_custom_tabulartogseagct_csv {
    
    infile = prepareExample(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], 'test.csv', strip_columns = 1)
    input = [ [ id:'test' ], infile]    

    CUSTOM_TABULARTOGSEAGCT ( input )
}
