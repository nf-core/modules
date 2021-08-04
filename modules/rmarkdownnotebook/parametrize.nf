import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions

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
