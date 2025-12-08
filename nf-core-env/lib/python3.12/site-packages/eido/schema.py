from logging import getLogger

from peppy.utils import load_yaml

from .const import SAMPLES_KEY, PROP_KEY

_LOGGER = getLogger(__name__)


def preprocess_schema(schema_dict):
    """
    Preprocess schema before validation for user's convenience

    Preprocessing includes:
    - renaming 'samples' to '_samples' since in the peppy.Project object
        _samples attribute holds the list of peppy.Samples objects.
    - adding array of strings entry for every string specified to accommodate
        subsamples in peppy.Project

    :param dict schema_dict: schema dictionary to preprocess
    :return dict: preprocessed schema
    """
    _LOGGER.debug(f"schema ori: {schema_dict}")
    if "project" not in schema_dict[PROP_KEY]:
        _LOGGER.debug("No project section found in schema")

    if SAMPLES_KEY in schema_dict[PROP_KEY]:
        if (
            "items" in schema_dict[PROP_KEY][SAMPLES_KEY]
            and PROP_KEY in schema_dict[PROP_KEY][SAMPLES_KEY]["items"]
        ):
            s_props = schema_dict[PROP_KEY][SAMPLES_KEY]["items"][PROP_KEY]
            for prop, val in s_props.items():
                if "type" in val and val["type"] in ["string", "number", "boolean"]:
                    s_props[prop] = {}
                    s_props[prop]["anyOf"] = [val, {"type": "array", "items": val}]
    else:
        _LOGGER.debug("No samples section found in schema")
    _LOGGER.debug(f"schema processed: {schema_dict}")
    return schema_dict


def read_schema(schema):
    """
    Safely read schema from YAML-formatted file.

    If the schema imports any other schemas, they will be read recursively.

    :param str | Mapping schema: path to the schema file
        or schema in a dict form
    :return list[dict]: read schemas
    :raise TypeError: if the schema arg is neither a Mapping nor a file path or
        if the 'imports' sections in any of the schemas is not a list
    """

    def _recursively_read_schemas(x, lst):
        if "imports" in x:
            if isinstance(x["imports"], list):
                for sch in x["imports"]:
                    lst.extend(read_schema(sch))
            else:
                raise TypeError("In schema the 'imports' section has to be a list")
        lst.append(x)
        return lst

    schema_list = []
    if isinstance(schema, str):
        _LOGGER.debug(f"Reading schema: {schema}")
        schema = load_yaml(schema)
    if not isinstance(schema, dict):
        raise TypeError(
            f"schema has to be a dict, path to an existing file or URL to a remote one. "
            f"Got: {type(schema)}"
        )
    return _recursively_read_schemas(schema, schema_list)
