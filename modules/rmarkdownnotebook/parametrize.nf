import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions


/**
 * Multiline code blocks need to have the same indentation level
 * as the `script:` section. This function re-indents code to the specified level.
 */
def indent_code_block(code, n_spaces) {
    def indent_str = " ".multiply(n_spaces)
    return code.stripIndent().split("\n").join("\n" + indent_str)
}

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

    // Writing the .params.yml file directly as follows does not work.
    // It only works in 'exec:', but not if there is a `script:` section:
    // task.workDir.resolve('.params.yml').text = yaml_str

    // Therefore, we inject it into the bash script:
    return """\
        cat <<"END_PARAMS_SECTION" > ./.params.yml
        ${indent_code_block(yaml_str, 8)}
        END_PARAMS_SECTION
    """
}
