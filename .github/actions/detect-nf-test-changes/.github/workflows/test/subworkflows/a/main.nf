include { 
    PROCESS_A;
    PROCESS_B
} from '../../modules/a/main.nf'

workflow SUBWORKFLOW_A {
    PROCESS_A()
}