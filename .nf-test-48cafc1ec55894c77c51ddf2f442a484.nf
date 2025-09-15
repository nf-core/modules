import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES  as VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL } from '/Users/famke/02-nf-core/modules/modules/nf-core/varlociraptor/callvariants/tests/../../estimatealignmentproperties/main.nf'

include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES  as VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_TUMOR } from '/Users/famke/02-nf-core/modules/modules/nf-core/varlociraptor/callvariants/tests/../../estimatealignmentproperties/main.nf'

include { VARLOCIRAPTOR_PREPROCESS  as VARLOCIRAPTOR_PREPROCESS_NORMAL } from '/Users/famke/02-nf-core/modules/modules/nf-core/varlociraptor/callvariants/tests/../../preprocess/main.nf'

include { VARLOCIRAPTOR_PREPROCESS  as VARLOCIRAPTOR_PREPROCESS_TUMOR } from '/Users/famke/02-nf-core/modules/modules/nf-core/varlociraptor/callvariants/tests/../../preprocess/main.nf'


// include test process
include { VARLOCIRAPTOR_CALLVARIANTS } from '/Users/famke/02-nf-core/modules/modules/nf-core/varlociraptor/callvariants/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    
    {
        def input = []
        
                input[0] = [
                    [id:'test_normal'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam', checkIfExists:true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai', checkIfExists:true),
                    ]
                input[1] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta', checkIfExists:true)
                    ]
                input[2] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta.fai', checkIfExists:true)
                    ]
                
        VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL(*input)
    }
    
    {
        def input = []
        
                input[0] = [
                    [id:'test_tumor'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam', checkIfExists:true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam.bai', checkIfExists:true),
                    ]
                input[1] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta', checkIfExists:true)
                    ]
                input[2] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta.fai', checkIfExists:true)
                    ]
                
        VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_TUMOR(*input)
    }
    
    {
        def input = []
        
                input[0] = Channel.of([
                    [id:'test_normal'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam', checkIfExists:true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai', checkIfExists:true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/gvcf/test.genome.vcf', checkIfExists:true),
                    ]).collect().join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json)
                input[1] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta', checkIfExists:true)
                    ]
                input[2] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta.fai', checkIfExists:true)
                    ]
                
        VARLOCIRAPTOR_PREPROCESS_NORMAL(*input)
    }
    
    {
        def input = []
        
                input[0] = Channel.of([
                    [id:'test_tumor'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam', checkIfExists:true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam.bai', checkIfExists:true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/gvcf/test2.genome.vcf', checkIfExists:true),
                    ]).collect().join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_TUMOR.out.alignment_properties_json)
                input[1] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta', checkIfExists:true)
                    ]
                input[2] = [
                    [id:'test'],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta.fai', checkIfExists:true)
                    ]
                
        VARLOCIRAPTOR_PREPROCESS_TUMOR(*input)
    }
    

    // process mapping
    def input = []
    
                input[0] = VARLOCIRAPTOR_PREPROCESS_NORMAL.out.bcf.map{meta1,vcf->[meta1,vcf,[]]}
                input[1] = Channel.of(file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/varlociraptor/scenario.yml', checkIfExists:true))
                input[2] = "normal"
                
    //----

    //run process
    VARLOCIRAPTOR_CALLVARIANTS(*input)

    if (VARLOCIRAPTOR_CALLVARIANTS.output){

        // consumes all named output channels and stores items in a json file
        for (def name in VARLOCIRAPTOR_CALLVARIANTS.out.getNames()) {
            serializeChannel(name, VARLOCIRAPTOR_CALLVARIANTS.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = VARLOCIRAPTOR_CALLVARIANTS.out as Object[]
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
