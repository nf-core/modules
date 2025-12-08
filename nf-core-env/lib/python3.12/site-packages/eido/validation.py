import os
from typing import NoReturn, Mapping, Union
from copy import deepcopy as dpcpy
from logging import getLogger

from warnings import warn

from .exceptions import EidoValidationError

from pandas.core.common import flatten
from jsonschema import Draft7Validator
import peppy

from .const import PROP_KEY, SIZING_KEY, TANGIBLE_KEY, SAMPLES_KEY
from .exceptions import PathAttrNotFoundError
from .schema import preprocess_schema, read_schema

_LOGGER = getLogger(__name__)


def _validate_object(obj: Mapping, schema: Union[str, dict], sample_name_colname=False):
    """
    Generic function to validate object against a schema

    :param Mapping obj: an object to validate
    :param str | dict schema: schema dict to validate against or a path to one
        from the error. Useful when used ith large projects

    :raises EidoValidationError: if validation is unsuccessful
    """
    validator = Draft7Validator(schema)
    _LOGGER.debug(f"{obj},\n {schema}")
    if not validator.is_valid(obj):
        errors = sorted(validator.iter_errors(obj), key=lambda e: e.path)
        errors_by_type = {}

        # Accumulate and restructure error objects by error type
        for error in errors:
            if not error.message in errors_by_type:
                errors_by_type[error.message] = []

            try:
                instance_name = error.instance[sample_name_colname]
            except KeyError:
                instance_name = "project"
            except TypeError:
                instance_name = obj["samples"][error.absolute_path[1]][
                    sample_name_colname
                ]
            errors_by_type[error.message].append(
                {
                    "type": error.message,
                    "message": f"{error.message} on instance {instance_name}",
                    "sample_name": instance_name,
                }
            )

        raise EidoValidationError("Validation failed", errors_by_type)
    else:
        _LOGGER.debug("Validation was successful...")


def validate_project(project: peppy.Project, schema: Union[str, dict]) -> NoReturn:
    """
    Validate a project object against a schema

    :param peppy.Project project: a project object to validate
    :param str | dict schema: schema dict to validate against or a path to one
    from the error. Useful when used ith large projects

    :return: NoReturn
    :raises EidoValidationError: if validation is unsuccessful
    """
    sample_name_colname = project.sample_name_colname
    schema_dicts = read_schema(schema=schema)
    for schema_dict in schema_dicts:
        project_dict = project.to_dict()
        _validate_object(
            project_dict, preprocess_schema(schema_dict), sample_name_colname
        )
        _LOGGER.debug("Project validation successful")


def _validate_sample_object(sample: peppy.Sample, schemas):
    """
    Internal function that allows to validate a peppy.Sample object without
    requiring a reference to peppy.Project.

    :param peppy.Sample sample: a sample object to validate
    :param list[dict] schemas: list of schemas to validate against or a path to one
    """
    for schema_dict in schemas:
        schema_dict = preprocess_schema(schema_dict)
        sample_schema_dict = schema_dict[PROP_KEY][SAMPLES_KEY]["items"]
        _validate_object(sample.to_dict(), sample_schema_dict)
        _LOGGER.debug(
            f"{getattr(sample, 'sample_name', '')} sample validation successful"
        )


def validate_sample(
    project: peppy.Project, sample_name: Union[str, int], schema: Union[str, dict]
) -> NoReturn:
    """
    Validate the selected sample object against a schema

    :param peppy.Project project: a project object to validate
    :param str | int sample_name: name or index of the sample to validate
    :param str | dict schema: schema dict to validate against or a path to one

    :raises EidoValidationError: if validation is unsuccessful
    """
    sample = (
        project.samples[sample_name]
        if isinstance(sample_name, int)
        else project.get_sample(sample_name)
    )
    _validate_sample_object(
        sample=sample,
        schemas=read_schema(schema=schema),
    )


def validate_config(
    project: Union[peppy.Project, dict], schema: Union[str, dict]
) -> NoReturn:
    """
    Validate the config part of the Project object against a schema

    :param peppy.Project project: a project object to validate
    :param str | dict schema: schema dict to validate against or a path to one
    """
    schema_dicts = read_schema(schema=schema)
    for schema_dict in schema_dicts:
        schema_cpy = preprocess_schema(dpcpy(schema_dict))
        try:
            del schema_cpy[PROP_KEY][SAMPLES_KEY]
        except KeyError:
            pass
        if "required" in schema_cpy:
            try:
                schema_cpy["required"].remove(SAMPLES_KEY)
            except ValueError:
                pass
        if isinstance(project, dict):
            _validate_object({"project": project}, schema_cpy)

        else:
            project_dict = project.to_dict()
            _validate_object(project_dict, schema_cpy)
            _LOGGER.debug("Config validation successful")


def _get_attr_values(obj, attrlist):
    """
    Get value corresponding to each given attribute.

    :param Mapping obj: an object to get the attributes from
    :param str | Iterable[str] attrlist: names of attributes to
        retrieve values for
    :return dict: value corresponding to
        each named attribute; null if this Sample's value for the
        attribute given by the argument to the "attrlist" parameter is
        empty/null, or if this Sample lacks the indicated attribute
    """
    # If attribute is None, then value is also None.
    if not attrlist:
        return None
    if not isinstance(attrlist, list):
        attrlist = [attrlist]
    # Strings contained here are appended later so shouldn't be null.
    return list(flatten([getattr(obj, attr, "") for attr in attrlist]))


def validate_input_files(
    project: peppy.Project,
    schemas: Union[str, dict],
    sample_name: Union[str, int] = None,
):
    """
    Determine which of the required and optional files are missing.

    The names of the attributes that are required and/or deemed as inputs
    are sourced from the schema, more specifically from `required_files`
    and `files` sections in samples section:

    - If any of the required files are missing, this function raises an error.
    - If any of the optional files are missing, the function raises a warning.

    Note, this function also performs Sample object validation with jsonschema.

    :param peppy.Project project: project that defines the samples to validate
    :param str | dict schema: schema dict to validate against or a path to one
    :param str | int sample_name: name or index of the sample to validate. If None,
        validate all samples in the project
    :raise PathAttrNotFoundError: if any required sample attribute is missing
    """

    if sample_name is None:
        samples = project.samples
    else:
        samples = (
            project.samples[sample_name]
            if isinstance(sample_name, int)
            else project.get_sample(sample_name)
        )
        samples = [samples]

    if isinstance(schemas, str):
        schemas = read_schema(schemas)

    for sample in samples:
        # validate attrs existence first
        _validate_sample_object(schemas=schemas, sample=sample)

        all_inputs = set()
        required_inputs = set()
        schema = schemas[-1]  # use only first schema, in case there are imports
        sample_schema_dict = schema[PROP_KEY][SAMPLES_KEY]["items"]
        if SIZING_KEY in sample_schema_dict:
            all_inputs.update(_get_attr_values(sample, sample_schema_dict[SIZING_KEY]))
        if TANGIBLE_KEY in sample_schema_dict:
            required_inputs = set(
                _get_attr_values(sample, sample_schema_dict[TANGIBLE_KEY])
            )
            all_inputs.update(required_inputs)

        missing_required_inputs = [i for i in required_inputs if not os.path.exists(i)]
        missing_inputs = [i for i in all_inputs if not os.path.exists(i)]
        if missing_inputs:
            warn(
                f"For sample '{getattr(sample, project.sample_table_index)}'. "
                f"Optional inputs not found: {missing_inputs}"
            )
        if missing_required_inputs:
            raise PathAttrNotFoundError(
                f"For sample '{getattr(sample, project.sample_table_index)}'. "
                f"Required inputs not found: {required_inputs}"
            )
