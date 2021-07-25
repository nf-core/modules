@Grab(group='org.yaml', module='snakeyaml', version='1.28')

import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions

// /**
//  * recursivly convert a groovy map to strings. Maps and Lists are descended into,
//  * all other objects are converted to string.
//  */
// def mapToString(map) {
//     return map.collectEntries {
//         key, val -> {
//             if (val instanceof Map) {
//                 [key, mapToString(val)]
//             } else if (val instanceof List) {
//                 [key, listToString(val)]
//             } else {
//                 [key, val.toString()]
//             }
//         }
//     }
// }


// /**
//  * recursivly convert a groovy list to strings. Maps and Lists are descended into,
//  * all other objects are converted to string.
//  */
// def listToString(list) {
//     return list.collect {
//         val -> {
//             if (val instanceof Map) {
//                 mapToString(val)
//             } else if (val instanceof List) {
//                 listToString(val)
//             } else {
//                 val.toString()
//             }
//         }
//     }
// }


/**
 * Create a config YAML file from a groovy map
 *
 * @params task The process' `task` variable
 * @returns a line to be inserted in the bash script.
 */
def dump_params_yml(params) {
    DumperOptions options = new DumperOptions();
    options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK);
    def yaml = new Yaml(options)
    def yaml_str = yaml.dump(params)

    // this does not work (it only works in 'exec:', but not if there is a `script:` section. )
    // task.workDir.resolve('.params.yml').text = yaml_str

    //therefore, we inject it into the bash script:
    return """cat <<"END_PARAMS_SECTION" > ./.params.yml\n${yaml_str}\nEND_PARAMS_SECTION\n\n"""
}
