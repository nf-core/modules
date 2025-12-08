"""
Constant variables for eido package
"""

LOGGING_LEVEL = "INFO"
PKG_NAME = "eido"
INSPECT_CMD = "inspect"
VALIDATE_CMD = "validate"
CONVERT_CMD = "convert"
FILTERS_CMD = "filters"
SUBPARSER_MSGS = {
    VALIDATE_CMD: "Validate a PEP or its components",
    INSPECT_CMD: "Inspect a PEP",
    CONVERT_CMD: "Convert PEP format using filters",
}
PROP_KEY = "properties"

SAMPLES_KEY = "samples"

TANGIBLE_KEY = "tangible"
SIZING_KEY = "sizing"

# sample schema input validation key names, these values are required by looper
# to refer to the dict values
MISSING_KEY = "missing"
REQUIRED_INPUTS_KEY = "required_inputs"
ALL_INPUTS_KEY = "all_inputs"
INPUT_FILE_SIZE_KEY = "input_file_size"

# groups of constants
GENERAL = [
    "LOGGING_LEVEL",
    "PKG_NAME",
    "INSPECT_CMD",
    "VALIDATE_CMD",
    "CONVERT_CMD",
    "FILTERS_CMD",
    "SUBPARSER_MSGS",
]

SCHEMA_SECTIONS = ["PROP_KEY", "TANGIBLE_KEY", "SIZING_KEY"]

SCHEMA_VALIDAION_KEYS = [
    "MISSING_KEY",
    "REQUIRED_INPUTS_KEY",
    "ALL_INPUTS_KEY",
    "INPUT_FILE_SIZE_KEY",
]

__all__ = GENERAL + SCHEMA_SECTIONS + SCHEMA_VALIDAION_KEYS
