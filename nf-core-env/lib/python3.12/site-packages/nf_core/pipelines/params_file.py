"""Create a YAML parameter file"""

import json
import logging
import textwrap
from pathlib import Path
from typing import Literal

import questionary

import nf_core.pipelines.list
import nf_core.utils
from nf_core.pipelines.schema import PipelineSchema

log = logging.getLogger(__name__)

INTRO = (
    "This is an example parameter file to pass to the `-params-file` option "
    "of nextflow run with the {pipeline_name} pipeline."
)

USAGE = "Uncomment lines with a single '#' if you want to pass the parameter to the pipeline."

H1_SEPERATOR = "## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
H2_SEPERATOR = "## ----------------------------------------------------------------------------"

ModeLiteral = Literal["both", "start", "end", "none"]


def _print_wrapped(text, fill_char="-", mode="both", width=80, indent=0, drop_whitespace=True) -> str:
    """Helper function to format text for the params-file template.

    Args:
        text (str): Text to print
        fill_char (str, optional):
            Character to use for creating dividers. Defaults to '-'.
        mode (str, optional):
            Where to place dividers. Defaults to "both".
        width (int, optional):
            Maximum line-width of the output text. Defaults to 80.
        indent (int, optional):
            Number of spaces to indent the text. Defaults to 0.
        drop_whitespace (bool, optional):
            Whether to drop whitespace from the start and end of lines.
    """
    if len(fill_char) != 1:
        raise ValueError("fill_char must be a single character")

    prefix = "## "
    out = ""

    if mode in ("both", "start"):
        out += prefix.ljust(width, fill_char) + "\n"

    wrap_indent = f"{prefix}{' ' * indent}"

    textlines = textwrap.wrap(
        text,
        width=width - len(prefix),
        initial_indent=wrap_indent,
        subsequent_indent=wrap_indent,
        drop_whitespace=drop_whitespace,
    )

    for line in textlines:
        out += line + "\n"

    if mode in ("both", "end"):
        out += prefix.ljust(width, fill_char) + "\n"

    return out


class ParamsFileBuilder:
    """Class to hold config option to launch a pipeline.

    Args:
        pipeline (str, optional):
            Path to a local pipeline path or a remote pipeline.
        revision (str, optional):
            Revision of the pipeline to use.
    """

    def __init__(
        self,
        pipeline=None,
        revision=None,
    ) -> None:
        """Initialise the ParamFileBuilder class

        Args:
            pipeline (str, optional): Path to a local pipeline path or a remote pipeline.
            revision (str, optional): Revision of the pipeline to use.
        """
        self.pipeline = pipeline
        self.pipeline_revision = revision
        self.schema_obj: PipelineSchema | None = None

        # Fetch remote workflows
        self.wfs = nf_core.pipelines.list.Workflows()
        self.wfs.get_remote_workflows()

    def get_pipeline(self) -> bool | None:
        """
        Prompt the user for a pipeline name and get the schema
        """
        # Prompt for pipeline if not supplied
        if self.pipeline is None:
            launch_type = questionary.select(
                "Generate parameter file for local pipeline or remote GitHub pipeline?",
                choices=["Remote pipeline", "Local path"],
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()

            if launch_type == "Remote pipeline":
                try:
                    self.pipeline = nf_core.utils.prompt_remote_pipeline_name(self.wfs)
                except AssertionError as e:
                    log.error(e.args[0])
                    return False
            else:
                self.pipeline = questionary.path(
                    "Path to workflow:", style=nf_core.utils.nfcore_question_style
                ).unsafe_ask()

        # Get the schema
        self.schema_obj = PipelineSchema()
        if self.schema_obj is None:
            return False
        self.schema_obj.get_schema_path(self.pipeline, local_only=False, revision=self.pipeline_revision)
        self.schema_obj.get_wf_params()
        return True

    def format_group(self, definition, show_hidden=False) -> str:
        """Format a group of parameters of the schema as commented YAML.

        Args:
            definition (dict): Definition of the group from the schema
            show_hidden (bool): Whether to include hidden parameters

        Returns:
            str: Formatted output for a group
        """
        out = ""
        title = definition.get("title", definition)
        description = definition.get("description", "")
        properties = definition.get("properties")

        hidden_props = set()
        for param_key, param_props in properties.items():
            if param_props.get("hidden", False):
                hidden_props.add(param_key)

        out += _print_wrapped(title, "=", mode="both")
        if description:
            out += _print_wrapped(description, mode="none")

        out += "\n"
        if len(hidden_props) > 0:
            out += _print_wrapped(f"({len(hidden_props)} hidden parameters are not shown)", mode="none")
            out += "\n\n"

        required_props = definition.get("required", [])

        for prop_key, props in properties.items():
            param_out = self.format_param(prop_key, props, required_props, show_hidden=show_hidden)
            if param_out is not None:
                out += param_out
                out += "\n"

        return out

    def format_param(
        self, name: str, properties: dict, required_properties: list[str] = [], show_hidden: bool = False
    ) -> str | None:
        """
        Format a single parameter of the schema as commented YAML

        Args:
            name (str): Name of the parameter
            properties (dict): Properties of the parameter
            required_properties (list): List of required properties
            show_hidden (bool): Whether to include hidden parameters

        Returns:
            str: Section of a params-file.yml for given parameter
            None: If the parameter is skipped because it is hidden and show_hidden is not set
        """
        out = ""
        hidden = properties.get("hidden", False)

        if not show_hidden and hidden:
            return None

        description = properties.get("description", "")
        if self.schema_obj is None:
            log.error("No schema object found")
            return ""
        self.schema_obj.get_schema_defaults()
        default = properties.get("default")
        type = properties.get("type")
        required = name in required_properties

        out += _print_wrapped(name, "-", mode="both")

        if description:
            out += _print_wrapped(description + "\n", mode="none", indent=4)

        if type:
            out += _print_wrapped(f"Type: {type}", mode="none", indent=4)

        if required:
            out += _print_wrapped("Required", mode="none", indent=4)

        out += _print_wrapped("\n", mode="end")
        out += f"# {name}: {json.dumps(default)}\n"

        return out

    def generate_params_file(self, show_hidden: bool = False) -> str:
        """Generate the contents of a parameter template file.

        Assumes the pipeline has been fetched (if remote) and the schema loaded.

        Args:
            show_hidden (bool): Whether to include hidden parameters

        Returns:
            str: Formatted output for the pipeline schema
        """
        if self.schema_obj is None:
            log.error("No schema object found")
            return ""

        schema = self.schema_obj.schema
        pipeline_name = self.schema_obj.pipeline_manifest.get("name", self.pipeline)
        pipeline_version = self.schema_obj.pipeline_manifest.get("version", "0.0.0")

        # Build the header section
        out = ""

        out += _print_wrapped(f"{pipeline_name} {pipeline_version}", "~", mode="both", indent=4)
        out += _print_wrapped(INTRO.format(pipeline_name=pipeline_name), " ", mode="none", indent=4)
        out += _print_wrapped("\n", " ", mode="none", indent=4, drop_whitespace=False)
        out += _print_wrapped(USAGE, "-", mode="end", indent=4)
        out += "\n"

        # Add all parameter groups
        for definition in schema.get("definitions", schema.get("$defs", {})).values():
            out += self.format_group(definition, show_hidden=show_hidden)
            out += "\n"

        return out

    def write_params_file(self, output_fn: Path = Path("nf-params.yaml"), show_hidden=False, force=False) -> bool:
        """Build a template file for the pipeline schema.

        Args:
            output_fn (str, optional): Filename to write the template to.
            show_hidden (bool, optional):
                Include parameters marked as hidden in the output
            force (bool, optional): Whether to overwrite existing output file.

        Returns:
            bool: True if the template was written successfully, False otherwise
        """

        self.get_pipeline()
        if self.schema_obj is None:
            log.error("No schema object found")
            return False
        try:
            self.schema_obj.load_schema()
            self.schema_obj.validate_schema()
        except AssertionError as e:
            log.error(f'Pipeline schema file is invalid ("{self.schema_obj.schema_filename}"): {e}')
            log.info("Please fix this file, then try again.")
            return False

        schema_out = self.generate_params_file(show_hidden=show_hidden)

        if output_fn.exists() and not force:
            log.error(f"File '{output_fn}' exists! Please delete first, or use '--force'")
            return False
        output_fn.write_text(schema_out)
        log.info(f"Parameter file written to '{output_fn}'")

        return True
