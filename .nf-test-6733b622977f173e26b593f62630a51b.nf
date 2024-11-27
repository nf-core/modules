import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { NACHO_QC } from '/workspace/modules/modules/nf-core/nacho/qc/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                // RCC Files: Collect from sample sheet
                input[0] =
                    Channel.fromPath('https://raw.githubusercontent.com/nf-core/test-datasets/nanostring/samplesheets/samplesheet_test.csv', checkIfExists: true)
                        .splitCsv( header: true )
                        .map{ row -> return file(row.RCC_FILE, checkIfExists: true) } // Select first column: path to file // Select first column: path to file
                        .collect()
                        .map{ files ->
                            tuple( [id: 'test'], files ) // Add meta component
                        }

                // Sample sheet
                input[1] = Channel.of( [
                        [ id: 'test_samplesheet'],
                        [ file('https://raw.githubusercontent.com/nf-core/test-datasets/nanostring/samplesheets/samplesheet_test.csv', checkIfExists: true) ]
                    ] )
                
    //----

    //run process
    NACHO_QC(*input)

    if (NACHO_QC.output){

        // consumes all named output channels and stores items in a json file
        for (def name in NACHO_QC.out.getNames()) {
            serializeChannel(name, NACHO_QC.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = NACHO_QC.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
