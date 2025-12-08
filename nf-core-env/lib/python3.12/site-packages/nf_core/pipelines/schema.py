"""Code to deal with pipeline JSON Schema"""

import copy
import json
import logging
import tempfile
import webbrowser
from pathlib import Path

import jinja2
import jsonschema
import markdown
import rich.console
import yaml
from rich.prompt import Confirm
from rich.syntax import Syntax

import nf_core.pipelines.list
import nf_core.utils
from nf_core.pipelines.lint_utils import dump_json_with_prettier, run_prettier_on_file

log = logging.getLogger(__name__)


class PipelineSchema:
    """Class to generate a schema object with
    functions to handle pipeline JSON Schema"""

    def __init__(self):
        """Initialise the object"""

        self.schema = {}
        self.pipeline_dir = ""
        self._schema_filename = ""
        self.schema_defaults = {}
        self.schema_types = {}
        self.schema_params = {}
        self.input_params = {}
        self.pipeline_params = {}
        self.invalid_nextflow_config_default_parameters = {}
        self.pipeline_manifest = {}
        self.schema_from_scratch = False
        self.no_prompts = False
        self.web_only = False
        self.web_schema_build_url = "https://oldsite.nf-co.re/pipeline_schema_builder"
        self.web_schema_build_web_url = None
        self.web_schema_build_api_url = None
        self.validation_plugin = None
        self.schema_draft = None
        self.defs_notation = None
        self.ignored_params = []

    # Update the validation plugin code every time the schema gets changed
    def set_schema_filename(self, schema: str) -> None:
        self._schema_filename = schema
        self._update_validation_plugin_from_config()

    def get_schema_filename(self) -> str:
        return self._schema_filename

    def del_schema_filename(self) -> None:
        del self._schema_filename

    schema_filename = property(get_schema_filename, set_schema_filename, del_schema_filename)

    def _update_validation_plugin_from_config(self) -> None:
        plugin = "nf-schema"
        if self.schema_filename:
            conf = nf_core.utils.fetch_wf_config(Path(self.schema_filename).parent)
        else:
            conf = nf_core.utils.fetch_wf_config(Path(self.pipeline_dir))

        plugins = str(conf.get("plugins", "")).strip("'\"").strip(" ").split(",")
        plugin_found = False
        for plugin_instance in plugins:
            if "nf-schema" in plugin_instance:
                plugin = "nf-schema"
                plugin_found = True
                break
            elif "nf-validation" in plugin_instance:
                plugin = "nf-validation"
                plugin_found = True
                break

        if not plugin_found:
            log.info(
                "Could not find nf-schema or nf-validation in the pipeline config. Defaulting to nf-schema notation for the JSON schema."
            )

        self.validation_plugin = plugin
        # Previous versions of nf-schema used "defs", but it's advised to use "$defs"
        if plugin == "nf-schema":
            self.defs_notation = "$defs"
            ignored_params = [
                conf.get("validation.help.shortParameter", "help"),
                conf.get("validation.help.fullParameter", "helpFull"),
                conf.get("validation.help.showHiddenParameter", "showHidden"),
                "trace_report_suffix",  # report suffix should be ignored by default as it is a Java Date object
            ]  # Help parameter should be ignored by default
            ignored_params_config_str = conf.get("validation.defaultIgnoreParams", "")
            ignored_params_config = [
                item.strip().strip("'") for item in ignored_params_config_str[1:-1].split(",")
            ]  # Extract list elements and remove whitespace

            if len(ignored_params_config) > 0:
                log.debug(f"Ignoring parameters from config: {ignored_params_config}")
                ignored_params.extend(ignored_params_config)
            self.ignored_params = ignored_params
            log.debug(f"Ignoring parameters: {self.ignored_params}")
            self.schema_draft = "https://json-schema.org/draft/2020-12/schema"

        else:
            self.defs_notation = "definitions"
            self.schema_draft = "https://json-schema.org/draft-07/schema"
            self.get_wf_params()
            self.ignored_params = self.pipeline_params.get("validationSchemaIgnoreParams", "").strip("\"'").split(",")
            self.ignored_params.append("validationSchemaIgnoreParams")

    def get_schema_path(self, path: str | Path, local_only: bool = False, revision: str | None = None) -> None:
        """Given a pipeline name, directory, or path, set self.schema_filename"""
        path = Path(path)
        # Supplied path exists - assume a local pipeline directory or schema
        if path.exists():
            log.debug(f"Path exists: {path}. Assuming local pipeline directory or schema")
            local_only = True
            if revision is not None:
                log.warning(f"Local workflow supplied, ignoring revision '{revision}'")
            if path.is_dir():
                self.pipeline_dir = path
                self.schema_filename = path / "nextflow_schema.json"
            else:
                self.pipeline_dir = path.parent
                self.schema_filename = path

        # Path does not exist - assume a name of a remote workflow
        elif not local_only:
            self.pipeline_dir = nf_core.pipelines.list.get_local_wf(path, revision=revision)
            self.schema_filename = Path(self.pipeline_dir or "", "nextflow_schema.json")
            # check if the schema file exists
            if not self.schema_filename.exists():
                self.schema_filename = None
        # Only looking for local paths, overwrite with None to be safe
        else:
            self.schema_filename = None

        # Check that the schema file exists
        if self.schema_filename is None or not Path(self.schema_filename).exists() and local_only:
            error = f"Could not find pipeline schema for '{path}': {self.schema_filename}"
            log.error(error)
            raise AssertionError(error)

    def load_lint_schema(self):
        """Load and lint a given schema to see if it looks valid"""
        try:
            self.load_schema()
            num_params = self.validate_schema()
            self.get_schema_defaults()
            self.validate_default_params()
            if len(self.invalid_nextflow_config_default_parameters) > 0:
                params = "\n --".join(
                    [f"{param}: {msg}" for param, msg in self.invalid_nextflow_config_default_parameters.items()]
                )
                log.info(
                    f"[red][✗] Invalid default parameters found:\n {params} \n\nNOTE: Use null in config for no default."
                )

            else:
                log.info(f"[green][✓] Pipeline schema looks valid[/] [dim](found {num_params} params)")
        except json.decoder.JSONDecodeError as e:
            error_msg = f"[bold red]Could not parse schema JSON:[/] {e}"
            log.error(error_msg)
            raise AssertionError(error_msg)
        except AssertionError as e:
            error_msg = f"[red][✗] Pipeline schema does not follow nf-core specs:\n {e}"
            log.error(error_msg)
            raise AssertionError(error_msg)

    def load_schema(self):
        """Load a pipeline schema from a file"""
        if self.schema_filename is None or not Path(self.schema_filename).exists():
            raise AssertionError("Pipeline schema filename could not be found.")

        with open(self.schema_filename) as fh:
            self.schema = json.load(fh)
        self.schema_defaults = {}
        self.schema_params = {}
        if "$schema" not in self.schema:
            raise AssertionError("Schema missing top-level `$schema` attribute")
        log.debug(f"JSON file loaded: {self.schema_filename}")

    def sanitise_param_default(self, param):
        """
        Given a param, ensure that the default value is the correct variable type
        """
        if "type" not in param or "default" not in param:
            return param

        # Bools
        if param["type"] == "boolean":
            if not isinstance(param["default"], bool):
                param["default"] = param["default"] == "true"
            return param

        # For everything else, an empty string is an empty string
        if isinstance(param["default"], str) and param["default"].strip() == "":
            param["default"] = ""
            return param

        # Integers
        if param["type"] == "integer":
            param["default"] = int(param["default"])
            return param

        # Numbers
        if param["type"] == "number":
            param["default"] = float(param["default"])
            return param

        if param["default"] is None:
            return param

        # Strings
        param["default"] = str(param["default"])
        return param

    def get_schema_defaults(self) -> None:
        """
        Generate set of default input parameters from schema.

        Saves defaults to self.schema_defaults
        Returns count of how many parameters were found (with or without a default value)
        """
        # Top level schema-properties (ungrouped)
        for p_key, param in self.schema.get("properties", {}).items():
            self.schema_params[p_key] = ("properties", p_key)
            if "default" in param:
                param = self.sanitise_param_default(param)
                if param["default"] is not None:
                    self.schema_defaults[p_key] = param["default"]

        # TODO add support for nested parameters
        # Grouped schema properties in subschema definitions
        for defn_name, definition in self.schema.get(self.defs_notation, {}).items():
            for p_key, param in definition.get("properties", {}).items():
                self.schema_params[p_key] = (self.defs_notation, defn_name, "properties", p_key)
                if "default" in param:
                    param = self.sanitise_param_default(param)
                    if param["default"] is not None:
                        self.schema_defaults[p_key] = param["default"]

    def get_schema_types(self) -> None:
        """Get a list of all parameter types in the schema"""
        for name, param in self.schema.get("properties", {}).items():
            if "type" in param:
                self.schema_types[name] = param["type"]
        for _, definition in self.schema.get(self.defs_notation, {}).items():
            for name, param in definition.get("properties", {}).items():
                if "type" in param:
                    self.schema_types[name] = param["type"]

    def save_schema(self, suppress_logging=False):
        """Save a pipeline schema to a file"""
        # Write results to a JSON file
        num_params = len(self.schema.get("properties", {}))
        num_params += sum(len(d.get("properties", {})) for d in self.schema.get(self.defs_notation, {}).values())
        if not suppress_logging:
            log.info(f"Writing schema with {num_params} params: '{self.schema_filename}'")
        dump_json_with_prettier(self.schema_filename, self.schema)

    def load_input_params(self, params_path):
        """Load a given a path to a parameters file (JSON/YAML)

        These should be input parameters used to run a pipeline with
        the Nextflow -params-file option.
        """
        # First, try to load as JSON
        try:
            with open(params_path) as fh:
                try:
                    params = json.load(fh)
                except json.JSONDecodeError as e:
                    raise UserWarning(f"Unable to load JSON file '{params_path}' due to error {e}")
                self.input_params.update(params)
            log.debug(f"Loaded JSON input params: {params_path}")
        except Exception as json_e:
            log.debug(f"Could not load input params as JSON: {json_e}")
            # This failed, try to load as YAML
            try:
                with open(params_path) as fh:
                    params = yaml.safe_load(fh)
                    self.input_params.update(params)
                    log.debug(f"Loaded YAML input params: {params_path}")
            except Exception as yaml_e:
                error_msg = f"Could not load params file as either JSON or YAML:\n JSON: {json_e}\n YAML: {yaml_e}"
                log.error(error_msg)
                raise AssertionError(error_msg)

    def validate_params(self):
        """Check given parameters against a schema and validate"""
        if self.schema is None:
            log.error("[red][✗] Pipeline schema not found")
            return False
        try:
            jsonschema.validate(self.input_params, self.schema)
        except jsonschema.exceptions.ValidationError as e:
            log.error(f"[red][✗] Input parameters are invalid: {e.message}")
            return False
        log.info("[green][✓] Input parameters look valid")
        return True

    def validate_default_params(self):
        """
        Check that all default parameters in the schema are valid
        Ignores 'required' flag, as required parameters might have no defaults

        Additional check that all parameters have defaults in nextflow.config and that
        these are valid and adhere to guidelines
        """
        if self.schema is None:
            log.error("[red][✗] Pipeline schema not found")
        try:
            jsonschema.validate(self.schema_defaults, strip_required(self.schema))
        except jsonschema.exceptions.ValidationError as e:
            log.debug(f"Complete error message:\n{e}")
            raise AssertionError(f"Default parameters are invalid: {e.message}")
        for param, default in self.schema_defaults.items():
            if default in ("null", "", None, "None") or default is False:
                log.warning(
                    f"[yellow][!] Default parameter '{param}' is empty, null, or False. It is advisable to remove the default from the schema"
                )
        log.info("[green][✓] Default parameters match schema validation")

        # Make sure every default parameter exists in the nextflow.config and is of correct type
        if self.pipeline_params == {}:
            self.get_wf_params()

        # Go over group keys
        for group_key, group in self.schema.get(self.defs_notation, {}).items():
            group_properties = group.get("properties")
            for param in group_properties:
                if param in self.ignored_params:
                    continue
                if param in self.pipeline_params:
                    self.validate_config_default_parameter(param, group_properties[param], self.pipeline_params[param])
                else:
                    self.invalid_nextflow_config_default_parameters[param] = (
                        "Not in pipeline parameters. Check `nextflow.config`."
                    )

        # Go over ungrouped params if any exist
        ungrouped_properties = self.schema.get("properties")
        if ungrouped_properties:
            for param in ungrouped_properties:
                if param in self.ignored_params:
                    continue
                if param in self.pipeline_params:
                    self.validate_config_default_parameter(
                        param, ungrouped_properties[param], self.pipeline_params[param]
                    )
                else:
                    self.invalid_nextflow_config_default_parameters[param] = (
                        "Not in pipeline parameters. Check `nextflow.config`."
                    )

    def validate_config_default_parameter(self, param, schema_param, config_default):
        """
        Assure that default parameters in the nextflow.config are correctly set
        by comparing them to their type in the schema
        """

        # If we have a default in the schema, check it matches the config
        if "default" in schema_param and (
            (schema_param["type"] == "boolean" and str(config_default).lower() != str(schema_param["default"]).lower())
            and (str(schema_param["default"]) != str(config_default).strip("'\""))
        ):
            # Check that we are not deferring the execution of this parameter in the schema default with squiggly brakcets
            if schema_param["type"] != "string" or "{" not in schema_param["default"]:
                self.invalid_nextflow_config_default_parameters[param] = (
                    f"Schema default (`{schema_param['default']}`) does not match the config default (`{config_default}`)"
                )
                return

        # if default is null, we're good
        if config_default == "null":
            return

        # Check variable types in nextflow.config
        if schema_param["type"] == "string":
            if str(config_default) in ["false", "true", "''"]:
                self.invalid_nextflow_config_default_parameters[param] = (
                    f"String should not be set to `{config_default}`"
                )
        if schema_param["type"] == "boolean":
            if str(config_default) not in ["false", "true"]:
                self.invalid_nextflow_config_default_parameters[param] = (
                    f"Booleans should only be true or false, not `{config_default}`"
                )
        if schema_param["type"] == "integer":
            try:
                int(config_default)
            except ValueError:
                self.invalid_nextflow_config_default_parameters[param] = (
                    f"Does not look like an integer: `{config_default}`"
                )
        if schema_param["type"] == "number":
            try:
                float(config_default)
            except ValueError:
                self.invalid_nextflow_config_default_parameters[param] = (
                    f"Does not look like a number (float): `{config_default}`"
                )

    def validate_schema(self, schema=None):
        """
        Check that the Schema is valid

        Returns: Number of parameters found
        """
        if schema is None:
            schema = self.schema

        if "$schema" not in schema:
            raise AssertionError("Schema missing top-level `$schema` attribute")
        schema_draft = schema["$schema"]
        if self.schema_draft != schema_draft:
            raise AssertionError(f"Schema is using the wrong draft: {schema_draft}, should be {self.schema_draft}")
        if self.schema_draft == "https://json-schema.org/draft-07/schema":
            try:
                jsonschema.Draft7Validator.check_schema(schema)
                log.debug("JSON Schema Draft7 validated")
            except jsonschema.exceptions.SchemaError as e:
                raise AssertionError(f"Schema does not validate as Draft 7 JSON Schema:\n {e}")
        elif self.schema_draft == "https://json-schema.org/draft/2020-12/schema":
            try:
                jsonschema.Draft202012Validator.check_schema(schema)
                log.debug("JSON Schema Draft2020-12 validated")
            except jsonschema.exceptions.SchemaError as e:
                raise AssertionError(f"Schema does not validate as Draft 2020-12 JSON Schema:\n {e}")
        else:
            raise AssertionError(
                f"Schema `$schema` should be `https://json-schema.org/draft/2020-12/schema` or `https://json-schema.org/draft-07/schema` \n Found `{schema_draft}`"
            )

        param_keys = list(schema.get("properties", {}).keys())
        num_params = len(param_keys)

        # Add a small check for older nf-schema JSON schemas
        if "defs" in schema:
            raise AssertionError(
                f'Using "defs" for schema definitions is not supported. Please use "{self.defs_notation}" instead'
            )

        for d_key, d_schema in schema.get(self.defs_notation, {}).items():
            # Check that this definition is mentioned in allOf
            if "allOf" not in schema:
                raise AssertionError("Schema has definitions, but no allOf key")
            in_allOf = False
            for allOf in schema.get("allOf", []):
                if allOf["$ref"] == f"#/{self.defs_notation}/{d_key}":
                    in_allOf = True
            if not in_allOf:
                raise AssertionError(
                    f"Definition subschema `#/{self.defs_notation}/{d_key}` not included in schema `allOf`"
                )

            # TODO add support for nested parameters
            for d_param_id in d_schema.get("properties", {}):
                # Check that we don't have any duplicate parameter IDs in different definitions
                if d_param_id in param_keys:
                    raise AssertionError(f"Duplicate parameter found in schema `{self.defs_notation}`: `{d_param_id}`")
                param_keys.append(d_param_id)
                num_params += 1

        # Check that everything in allOf exists
        for allOf in schema.get("allOf", []):
            _, allof_defs_notation, def_key = allOf["$ref"].split("/")  # "#/<defs_notation>/<def_name>"
            if allof_defs_notation not in schema:
                raise AssertionError(f"Schema has allOf, but no {allof_defs_notation}")
            if def_key not in schema.get(allof_defs_notation, {}):
                raise AssertionError(f"Subschema `{def_key}` found in `allOf` but not `{allof_defs_notation}`")

        # Check that the schema describes at least one parameter
        if num_params == 0:
            raise AssertionError("No parameters found in schema")

        return num_params

    def validate_schema_title_description(self, schema=None):
        """
        Extra validation command for linting.
        Checks that the schema "$id", "title" and "description" attributes match the pipeline config.
        """
        if schema is None:
            schema = self.schema
        if schema is None:
            log.debug("Pipeline schema not set - skipping validation of top-level attributes")
            return None

        if self.pipeline_manifest == {}:
            self.get_wf_params()

        if "name" not in self.pipeline_manifest:
            log.debug("Pipeline manifest `name` not known - skipping validation of schema id and title")
        else:
            if "$id" not in self.schema:
                raise AssertionError("Schema missing top-level `$id` attribute")
            if "title" not in self.schema:
                raise AssertionError("Schema missing top-level `title` attribute")
            # Validate that id, title and description match the pipeline manifest
            id_attr = "https://raw.githubusercontent.com/{}/main/nextflow_schema.json".format(
                self.pipeline_manifest["name"].strip("\"'")
            )
            if self.schema["$id"] not in [id_attr, id_attr.replace("/main/", "/master/")]:
                raise AssertionError(
                    f"Schema `$id` should be `{id_attr}` or {id_attr.replace('/main/', '/master/')}. \n Found `{self.schema['$id']}`"
                )

            title_attr = "{} pipeline parameters".format(self.pipeline_manifest["name"].strip("\"'"))
            if self.schema["title"] != title_attr:
                raise AssertionError(f"Schema `title` should be `{title_attr}`\n Found: `{self.schema['title']}`")

        if "description" not in self.pipeline_manifest:
            log.debug("Pipeline manifest 'description' not known - skipping validation of schema description")
        else:
            if "description" not in self.schema:
                raise AssertionError("Schema missing top-level 'description' attribute")
            desc_attr = self.pipeline_manifest["description"].strip("\"'")
            if self.schema["description"] != desc_attr:
                raise AssertionError(
                    f"Schema 'description' should be '{desc_attr}'\n Found: '{self.schema['description']}'"
                )

    def check_for_input_mimetype(self):
        """
        Check that the input parameter has a mimetype

        Common mime types: https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types

        Returns:
            mimetype (str): The mimetype of the input parameter

        Raises:
            LookupError: If the input parameter is not found or defined in the correct place
        """

        # Check that the input parameter is defined
        if "input" not in self.schema_params:
            raise LookupError("Parameter `input` not found in schema")
        # Check that the input parameter is defined in the right place
        if "input" not in self.schema.get(self.defs_notation, {}).get("input_output_options", {}).get("properties", {}):
            raise LookupError("Parameter `input` is not defined in the correct subschema (input_output_options)")
        input_entry = self.schema[self.defs_notation]["input_output_options"]["properties"]["input"]
        if "mimetype" not in input_entry:
            return None
        mimetype = input_entry["mimetype"]
        if mimetype == "" or mimetype is None:
            return None
        return mimetype

    def print_documentation(
        self,
        output_fn=None,
        format="markdown",
        force=False,
        columns=None,
    ):
        """
        Prints documentation for the schema.
        """
        if columns is None:
            columns = ["parameter", "description", "type,", "default", "required", "hidden"]

        output = self.schema_to_markdown(columns)
        if format == "html":
            output = self.markdown_to_html(output)

        with tempfile.NamedTemporaryFile(mode="w+") as fh:
            fh.write(output)
            run_prettier_on_file(fh.name)
            fh.seek(0)
            prettified_docs = fh.read()

        if not output_fn:
            console = rich.console.Console()
            console.print("\n", Syntax(prettified_docs, format, word_wrap=True), "\n")
        else:
            if Path(output_fn).exists() and not force:
                log.error(f"File '{output_fn}' exists! Please delete first, or use '--force'")
                return
            with open(output_fn, "w") as fh:
                fh.write(prettified_docs)
                log.info(f"Documentation written to '{output_fn}'")

        # Return as a string
        return output

    def schema_to_markdown(self, columns):
        """
        Creates documentation for the schema in Markdown format.
        """
        out = f"# {self.schema['title']}\n\n"
        out += f"{self.schema['description']}\n"
        # Grouped parameters
        for definition in self.schema.get(self.defs_notation, {}).values():
            out += f"\n## {definition.get('title', {})}\n\n"
            out += f"{definition.get('description', '')}\n\n"
            required = definition.get("required", [])
            properties = definition.get("properties", {})
            param_table = self.markdown_param_table(properties, required, columns)
            out += param_table

        # Top-level ungrouped parameters
        if len(self.schema.get("properties", {})) > 0:
            out += "\n## Other parameters\n\n"
            required = self.schema.get("required", [])
            properties = self.schema.get("properties", {})
            param_table = self.markdown_param_table(properties, required, columns)
            out += param_table

        return out

    def markdown_param_table(self, properties, required, columns):
        """Creates a markdown table for params from jsonschema properties section

        Args:
            properties (dict): A jsonschema properties dictionary
            required (list): A list of the required fields.
                Should come from the same level of the jsonschema as properties
            columns (list): A list of columns to write

        Returns:
            str: A string with the markdown table
        """
        out = ""
        out += "".join([f"| {column.title()} " for column in columns])
        out += "|\n"
        out += "".join(["|-----------" for _ in columns])
        out += "|\n"
        for p_key, param in properties.items():
            for column in columns:
                if column == "parameter":
                    out += f"| `{p_key}` "
                elif column == "description":
                    desc = param.get("description", "").replace("\n", "<br>")
                    if "enum" in param:
                        enum_values = "\\|".join(f"`{e}`" for e in param["enum"])
                        desc += f" (accepted: {enum_values})"
                    out += f"| {desc} "
                    if param.get("help_text", "") != "":
                        help_txt = param["help_text"].replace("\n", "<br>")
                        out += f"<details><summary>Help</summary><small>{help_txt}</small></details>"
                elif column == "type":
                    out += f"| `{param.get('type', '')}` "
                elif column == "required":
                    out += f"| {p_key in required or ''} "
                else:
                    out += f"| {param.get(column, '')} "
            out += "|\n"
        return out

    def markdown_to_html(self, markdown_str):
        """
        Convert markdown to html
        """
        return markdown.markdown(markdown_str, extensions=["tables"])

    def make_skeleton_schema(self):
        """Make a new pipeline schema from the template"""
        self.schema_from_scratch = True
        # Use Jinja to render the template schema file to a variable
        env = jinja2.Environment(
            loader=jinja2.PackageLoader("nf_core", "pipeline-template"), keep_trailing_newline=True
        )
        schema_template = env.get_template("nextflow_schema.json")

        template_vars = {
            "name": self.pipeline_manifest.get("name", Path(self.schema_filename).parent.name).strip("'"),
            "description": self.pipeline_manifest.get("description", "").strip("'"),
        }

        self.schema = json.loads(schema_template.render(template_vars))
        self.get_schema_defaults()

    def build_schema(self, pipeline_dir, no_prompts, web_only, url):
        """Interactively build a new pipeline schema for a pipeline"""

        # Check if supplied pipeline directory really is one
        nf_core.utils.is_pipeline_directory(pipeline_dir)

        if no_prompts:
            self.no_prompts = True
        if web_only:
            self.web_only = True
        if url:
            self.web_schema_build_url = url

        # Get pipeline schema filename
        try:
            self.get_schema_path(pipeline_dir, local_only=True)
        except AssertionError:
            log.info("No existing schema found - creating a new one from the nf-core template")
            self.get_wf_params()
            self.make_skeleton_schema()
            self.remove_schema_notfound_configs()
            self.remove_schema_empty_definitions()
            self.add_schema_found_configs()
            try:
                self.validate_schema()
            except AssertionError as e:
                log.error(f"[red]Something went wrong when building a new schema:[/] {e}")
                log.info("Please ask for help on the nf-core Slack")
                return False
        else:
            # Schema found - load and validate
            try:
                self.load_lint_schema()
            except AssertionError:
                log.error(f"Existing pipeline schema found, but it is invalid: {self.schema_filename}")
                log.info("Please fix or delete this file, then try again.")
                return False

        if not self.web_only:
            self.get_wf_params()
            self.remove_schema_notfound_configs()
            self.remove_schema_empty_definitions()
            self.add_schema_found_configs()
            self.save_schema()

        # If running interactively, send to the web for customisation
        if not self.no_prompts:
            if Confirm.ask(":rocket:  Launch web builder for customisation and editing?"):
                try:
                    self.launch_web_builder()
                except AssertionError as e:
                    log.error(e.args[0])
                    # Extra help for people running offline
                    if "Could not connect" in e.args[0]:
                        log.info(
                            f"If you're working offline, now copy your schema ({self.schema_filename}) and paste at https://nf-co.re/pipeline_schema_builder"
                        )
                        log.info("When you're finished, you can paste the edited schema back into the same file")
                    if self.web_schema_build_web_url:
                        log.info(
                            "To save your work, open {}\n"
                            f"Click the blue 'Finished' button, copy the schema and paste into this file: {self.web_schema_build_web_url, self.schema_filename}"
                        )
                    return False

    def get_wf_params(self):
        """
        Load the pipeline parameter defaults using `nextflow config`
        Strip out only the params. values and ignore anything that is not a flat variable
        """
        # Check that we haven't already pulled these (eg. skeleton schema)
        if len(self.pipeline_params) > 0 and len(self.pipeline_manifest) > 0:
            log.debug("Skipping get_wf_params as we already have them")
            return
        if self.schema_filename is None:
            log.error("Cannot get workflow params without a schema file")
            return
        log.debug("Collecting pipeline parameter defaults\n")
        config = nf_core.utils.fetch_wf_config(Path(self.schema_filename).parent)
        skipped_params = []
        # Pull out just the params. values
        for ckey, cval in config.items():
            if ckey.startswith("params."):
                # skip anything that's not a flat variable
                if "." in ckey[7:]:
                    skipped_params.append(ckey)
                    continue
                self.pipeline_params[ckey[7:]] = cval
            if ckey.startswith("manifest."):
                self.pipeline_manifest[ckey[9:]] = cval
        # Log skipped params
        if len(skipped_params) > 0:
            log.debug(
                "Skipped following pipeline params because they had nested parameter values:\n{}".format(
                    ", ".join(skipped_params)
                )
            )

    def remove_schema_empty_definitions(self):
        """
        Go through top-level schema remove definitions that don't have
        any property attributes
        """
        # Identify and remove empty definitions from the schema
        empty_definitions = []
        for d_key, d_schema in list(self.schema.get(self.defs_notation, {}).items()):
            if not d_schema.get("properties"):
                del self.schema[self.defs_notation][d_key]
                empty_definitions.append(d_key)
                log.warning(f"Removing empty group: '{d_key}'")

        # Remove "allOf" group with empty definitions from the schema
        for d_key in empty_definitions:
            allOf = {"$ref": f"#/{self.defs_notation}/{d_key}"}
            if allOf in self.schema.get("allOf", []):
                self.schema["allOf"].remove(allOf)

        # If we don't have anything left in "allOf", remove it
        if self.schema.get("allOf") == []:
            del self.schema["allOf"]

        # If we don't have anything left in "definitions", remove it
        if self.schema.get(self.defs_notation) == {}:
            del self.schema[self.defs_notation]

    def remove_schema_notfound_configs(self):
        """
        Go through top-level schema and all definitions sub-schemas to remove
        anything that's not in the nextflow config.
        """
        # Top-level properties
        self.schema, params_removed = self.remove_schema_notfound_configs_single_schema(self.schema)
        # Sub-schemas in definitions
        for d_key, definition in self.schema.get(self.defs_notation, {}).items():
            cleaned_schema, p_removed = self.remove_schema_notfound_configs_single_schema(definition)
            self.schema[self.defs_notation][d_key] = cleaned_schema
            params_removed.extend(p_removed)

        return params_removed

    def remove_schema_notfound_configs_single_schema(self, schema):
        """
        Go through a single schema / set of properties and strip out
        anything that's not in the nextflow config.

        Takes: Schema or sub-schema with properties key
        Returns: Cleaned schema / sub-schema
        """
        # Make a deep copy so as not to edit in place
        schema = copy.deepcopy(schema)
        params_removed = []
        # Use iterator so that we can delete the key whilst iterating
        for p_key in [k for k in schema.get("properties", {}).keys()]:
            if self.prompt_remove_schema_notfound_config(p_key):
                del schema["properties"][p_key]
                # Remove required flag if set
                if p_key in schema.get("required", []):
                    schema["required"].remove(p_key)
                # Remove required list if now empty
                if "required" in schema and len(schema["required"]) == 0:
                    del schema["required"]
                log.debug(f"Removing '{p_key}' from pipeline schema")
                params_removed.append(p_key)

        return schema, params_removed

    def prompt_remove_schema_notfound_config(self, p_key):
        """
        Check if a given key is found in the nextflow config params and prompt to remove it if note

        Returns True if it should be removed, False if not.
        """
        if p_key not in self.pipeline_params:
            if self.no_prompts or self.schema_from_scratch:
                return True
            if Confirm.ask(
                f":question: Unrecognised [bold]'params.{p_key}'[/] found in the schema but not in the pipeline config! [yellow]Remove it?"
            ):
                return True
        return False

    def add_schema_found_configs(self):
        """
        Add anything that's found in the Nextflow params that's missing in the pipeline schema
        Update defaults if they have changed
        """
        params_added = []

        for p_key, p_val in self.pipeline_params.items():
            s_key = self.schema_params.get(p_key)
            # Check if key is in schema parameters
            # Key is in pipeline but not in schema or ignored from schema
            if p_key not in self.schema_params and p_key not in self.ignored_params:
                if (
                    self.no_prompts
                    or self.schema_from_scratch
                    or Confirm.ask(
                        f":sparkles: Found [bold]'params.{p_key}'[/] in the pipeline config, but not in the schema. "
                        "[blue]Add to pipeline schema?"
                    )
                ):
                    if "properties" not in self.schema:
                        self.schema["properties"] = {}
                    self.schema["properties"][p_key] = self.build_schema_param(p_val)
                    log.debug(f"Adding '{p_key}' to pipeline schema")
                    params_added.append(p_key)
            # Param has a default that does not match the schema
            elif p_key in self.schema_defaults and (s_def := self.schema_defaults[p_key]) != (
                p_def := self.build_schema_param(p_val).get("default")
            ):
                if self.no_prompts or Confirm.ask(
                    f":sparkles: Default for [bold]'params.{p_key}'[/] in the pipeline config does not match schema. (schema: '{type(s_def)}: {s_def}'  | config: '{type(p_def)}: {p_def}'). "
                    "[blue]Update pipeline schema?"
                ):
                    s_key_def = s_key + ("default",)
                    if p_def is None:
                        nf_core.utils.nested_delitem(self.schema, s_key_def)
                        log.debug(f"Removed '{p_key}' default from pipeline schema")
                    else:
                        nf_core.utils.nested_setitem(self.schema, s_key_def, p_def)
                        log.debug(f"Updating '{p_key}' default to '{p_def}' in pipeline schema")
            # There is no default in schema but now there is a default to write
            elif (
                s_key
                and (p_key not in self.schema_defaults)
                and (p_key not in self.ignored_params)
                and (p_def := self.build_schema_param(p_val).get("default"))
            ):
                if self.no_prompts or Confirm.ask(
                    f":sparkles: Default for [bold]'params.{p_key}'[/] is not in schema (def='{p_def}'). "
                    "[blue]Update pipeline schema?"
                ):
                    s_key_def = s_key + ("default",)
                    nf_core.utils.nested_setitem(self.schema, s_key_def, p_def)
                    log.debug(f"Updating '{p_key}' default to '{p_def}' in pipeline schema")
        return params_added

    def build_schema_param(self, p_val):
        """
        Build a pipeline schema dictionary for an param interactively
        """
        p_val = p_val.strip("\"'")
        # p_val is always a string as it is parsed from nextflow config this way
        try:
            p_val = float(p_val)
            if p_val == int(p_val):
                p_val = int(p_val)
                p_type = "integer"
            else:
                p_type = "number"
        except ValueError:
            p_type = "string"

        # Anything can be "null", means that it is not set
        if p_val == "null":
            p_val = None

        # Booleans
        if p_val in ["true", "false", "True", "False"]:
            p_val = p_val in ["true", "True"]  # Convert to bool
            p_type = "boolean"

        # Don't return a default for anything false-y except 0
        if not p_val and not (p_val == 0 and p_val is not False):
            return {"type": p_type}

        return {"type": p_type, "default": p_val}

    def launch_web_builder(self):
        """
        Send pipeline schema to web builder and wait for response
        """

        content = {
            "post_content": "json_schema",
            "api": "true",
            "version": nf_core.__version__,
            "status": "waiting_for_user",
            "schema": json.dumps(self.schema),
        }
        web_response = nf_core.utils.poll_nfcore_web_api(self.web_schema_build_url, content)

        try:
            if "api_url" not in web_response:
                raise AssertionError('"api_url" not in web_response')
            if "web_url" not in web_response:
                raise AssertionError('"web_url" not in web_response')
            # DO NOT FIX THIS TYPO. Needs to stay in sync with the website. Maintaining for backwards compatibility.
            if web_response["status"] != "recieved":
                raise AssertionError(
                    f'web_response["status"] should be "recieved", but it is "{web_response["status"]}"'
                )
        except AssertionError:
            log.debug(f"Response content:\n{json.dumps(web_response, indent=4)}")
            raise AssertionError(
                f"Pipeline schema builder response not recognised: {self.web_schema_build_url}\n"
                " See verbose log for full response (nf-core -v schema)"
            )
        else:
            self.web_schema_build_web_url = web_response["web_url"]
            self.web_schema_build_api_url = web_response["api_url"]
            log.info(f"Opening URL: {web_response['web_url']}")
            webbrowser.open(web_response["web_url"])
            log.info("Waiting for form to be completed in the browser. Remember to click Finished when you're done.\n")
            nf_core.utils.wait_cli_function(self.get_web_builder_response)

    def get_web_builder_response(self):
        """
        Given a URL for a Schema build response, recursively query it until results are ready.
        Once ready, validate Schema and write to disk.
        """
        web_response = nf_core.utils.poll_nfcore_web_api(self.web_schema_build_api_url)
        if web_response["status"] == "error":
            raise AssertionError(f"Got error from schema builder: '{web_response.get('message')}'")
        if web_response["status"] == "waiting_for_user":
            return False
        if web_response["status"] == "web_builder_edited":
            log.info("Found saved status from nf-core pipelines schema builder")
            try:
                self.schema = web_response["schema"]
                self.remove_schema_empty_definitions()
                self.validate_schema()
            except AssertionError as e:
                raise AssertionError(f"Response from schema builder did not pass validation:\n {e}")
            else:
                self.save_schema()
                return True
        else:
            log.debug(f"Response content:\n{json.dumps(web_response, indent=4)}")
            raise AssertionError(
                f"Pipeline schema builder returned unexpected status ({web_response['status']}): "
                f"{self.web_schema_build_api_url}\n See verbose log for full response"
            )


def strip_required(node):
    if isinstance(node, dict):
        return {
            k: y
            for k, v in node.items()
            for y in [strip_required(v)]
            if k != "required" and (y or y is False or y == "")
        }
    elif isinstance(node, list):
        return [y for v in node for y in [strip_required(v)] if y or y is False or y == ""]
    else:
        return node
