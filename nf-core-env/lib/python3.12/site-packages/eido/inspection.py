import os
from logging import getLogger

from warnings import catch_warnings
from ubiquerg import size


from .const import (
    ALL_INPUTS_KEY,
    INPUT_FILE_SIZE_KEY,
    MISSING_KEY,
    PROP_KEY,
    REQUIRED_INPUTS_KEY,
    SIZING_KEY,
    TANGIBLE_KEY,
    SAMPLES_KEY,
)
from .schema import read_schema
from .validation import _validate_sample_object, _get_attr_values

_LOGGER = getLogger(__name__)


def inspect_project(p, sample_names=None, max_attr=10):
    """
    Print inspection info: Project or,
    if sample_names argument is provided, matched samples

    :param peppy.Project p: project to inspect
    :param Iterable[str] sample_names: list of samples to inspect
    :param int max_attr: max number of sample attributes to display
    """
    if sample_names:
        samples = p.get_samples(sample_names)
        if not samples:
            print("No samples matched by names: {}".format(sample_names))
            return
        for s in samples:
            print(s.__str__(max_attr=max_attr))
            print("\n")
        return
    print(p)
    return


def get_input_files_size(sample, schema):
    """
    Determine which of this Sample's required attributes/files are missing
    and calculate sizes of the files (inputs).

    The names of the attributes that are required and/or deemed as inputs
    are sourced from the schema, more specifically from required_input_attrs
    and input_attrs sections in samples section. Note, this function does
    perform actual Sample object validation with jsonschema.

    :param peppy.Sample sample: sample to investigate
    :param list[dict] | str schema: schema dict to validate against or a path to one
    :return dict: dictionary with validation data, i.e missing,
        required_inputs, all_inputs, input_file_size
    :raise ValidationError: if any required sample attribute is missing
    """
    if isinstance(schema, str):
        schema = read_schema(schema)

    # first, validate attrs existence using jsonschema
    _validate_sample_object(schemas=schema, sample=sample)

    all_inputs = set()
    required_inputs = set()
    schema = schema[-1]  # use only first schema, in case there are imports
    sample_schema_dict = schema[PROP_KEY][SAMPLES_KEY]["items"]
    if SIZING_KEY in sample_schema_dict:
        all_inputs.update(_get_attr_values(sample, sample_schema_dict[SIZING_KEY]))
    if TANGIBLE_KEY in sample_schema_dict:
        required_inputs = set(
            _get_attr_values(sample, sample_schema_dict[TANGIBLE_KEY])
        )
        all_inputs.update(required_inputs)
    with catch_warnings(record=True) as w:
        input_file_size = sum(
            [
                size(f, size_str=False) or 0.0
                for f in all_inputs
                if f != "" and f != None
            ]
        ) / (1024**3)
        if w:
            _LOGGER.warning(
                f"{len(w)} input files missing, job input size was "
                f"not calculated accurately"
            )

    return {
        MISSING_KEY: [i for i in required_inputs if not os.path.exists(i)],
        REQUIRED_INPUTS_KEY: required_inputs,
        ALL_INPUTS_KEY: all_inputs,
        INPUT_FILE_SIZE_KEY: input_file_size,
    }
