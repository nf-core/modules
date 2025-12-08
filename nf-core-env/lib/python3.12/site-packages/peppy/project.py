"""
Build a Project object.
"""

import os
import sys
from collections.abc import Mapping, MutableMapping
from contextlib import suppress
from logging import getLogger
from typing import Iterable, List, Tuple, Union, Literal

import numpy as np
import pandas as pd
import yaml
from pandas.core.common import flatten
from rich.console import Console
from rich.progress import track
from ubiquerg import is_url
from copy import deepcopy

from .const import (
    ACTIVE_AMENDMENTS_KEY,
    AMENDMENTS_KEY,
    ATTR_KEY_PREFIX,
    CFG_IMPORTS_KEY,
    CFG_SAMPLE_TABLE_KEY,
    CFG_SUBSAMPLE_TABLE_KEY,
    CONFIG_FILE_KEY,
    CONFIG_KEY,
    CONFIG_VERSION_KEY,
    APPEND_KEY,
    DERIVED_ATTRS_KEY,
    DERIVED_KEY,
    DERIVED_SOURCES_KEY,
    DESC_KEY,
    DUPLICATED_KEY,
    IMPLIED_COND_KEYS,
    IMPLIED_IF_KEY,
    IMPLIED_KEY,
    IMPLIED_THEN_KEY,
    MAX_PROJECT_SAMPLES_REPR,
    METADATA_KEY,
    NAME_KEY,
    PEP_LATEST_VERSION,
    PKG_NAME,
    PROJ_MODS_KEY,
    REMOVE_KEY,
    REQUIRED_VERSION,
    SAMPLE_DF_KEY,
    SAMPLE_DF_LARGE,
    SAMPLE_EDIT_FLAG_KEY,
    SAMPLE_MODIFIERS,
    SAMPLE_MODS_KEY,
    SAMPLE_NAME_ATTR,
    SAMPLE_RAW_DICT_KEY,
    SAMPLE_TABLE_FILE_KEY,
    SAMPLE_TABLE_INDEX_KEY,
    SUBSAMPLE_DF_KEY,
    SUBSAMPLE_NAME_ATTR,
    SUBSAMPLE_RAW_LIST_KEY,
    SUBSAMPLE_TABLE_INDEX_KEY,
    SUBSAMPLE_TABLES_FILE_KEY,
    ORIGINAL_CONFIG_KEY,
)
from .exceptions import (
    InvalidSampleTableFileException,
    MissingAmendmentError,
    IllegalStateException,
    InvalidConfigFileException,
)
from .parsers import select_parser
from .sample import Sample
from .utils import (
    copy,
    is_cfg_or_anno,
    load_yaml,
    make_abs_via_cfg,
    make_list,
)

_LOGGER = getLogger(PKG_NAME)


@copy
class Project(MutableMapping):
    """
    A class to model a Project (collection of samples and metadata).

    :param str cfg: Project config file (YAML) or sample table (CSV/TSV)
        with one row per sample to constitute project
    :param str | Iterable[str] sample_table_index: name of the columns to set
        the sample_table index to
    :param str | Iterable[str] subsample_table_index: name of the columns to set
        the subsample_table index to
    :param str | Iterable[str] amendments: names of the amendments to activate
    :param Iterable[str] amendments: amendments to use within configuration file
    :param bool defer_samples_creation: whether the sample creation should be skipped

    :Example:

    .. code-block:: python

        from peppy import Project
        prj = Project(cfg="ngs.yaml")
        samples = prj.samples
    """

    def __init__(
        self,
        cfg: str = None,
        amendments: Union[str, Iterable[str]] = None,
        sample_table_index: Union[str, Iterable[str]] = None,
        subsample_table_index: Union[str, Iterable[str]] = None,
        defer_samples_creation: bool = False,
    ):
        _LOGGER.debug(
            "Creating {}{}".format(
                self.__class__.__name__, " from file {}".format(cfg) if cfg else ""
            )
        )
        self._project_data = {}
        super(Project, self).__init__()
        is_cfg = is_cfg_or_anno(cfg)
        if is_cfg is None:
            # no 'cfg' provided. Empty Project will be created
            self[CONFIG_FILE_KEY] = None
            self[SAMPLE_TABLE_FILE_KEY] = None
            self[SUBSAMPLE_TABLES_FILE_KEY] = None
        elif is_cfg:
            # the provided 'cfg' is a project config file
            self[CONFIG_FILE_KEY] = cfg
            self[SAMPLE_TABLE_FILE_KEY] = None
            self[SUBSAMPLE_TABLES_FILE_KEY] = None
            self.parse_config_file(cfg, amendments)
        else:
            # the provided 'cfg' is a sample table
            self[SAMPLE_TABLE_FILE_KEY] = cfg
            self[SUBSAMPLE_TABLES_FILE_KEY] = None

        self._samples = []
        self[SAMPLE_EDIT_FLAG_KEY] = False
        self.progressbar = False

        # table indexes can be specified in config or passed to the object constructor
        # That's the priority order:
        # 1. constructor specified
        # 2. config specified (already set as Project attrs if config exists)
        # 3. defaults
        self.st_index = (
            sample_table_index or getattr(self, "st_index", None) or SAMPLE_NAME_ATTR
        )

        self.sst_index = (
            ([subsample_table_index] if subsample_table_index else None)
            or (
                [getattr(self, "sst_index", None)]
                if getattr(self, "sst_index", None)
                else None
            )
            or [SAMPLE_NAME_ATTR, SUBSAMPLE_NAME_ATTR]
        )

        if not defer_samples_creation:
            self.create_samples(modify=False if self[SAMPLE_TABLE_FILE_KEY] else True)
        self._sample_table = self._get_table_from_samples(
            index=self.st_index, initial=True
        )

    def __eq__(self, other):
        return [s.to_dict() for s in self.samples] == [
            s.to_dict() for s in other.samples
        ]

    @classmethod
    def from_pandas(
        cls,
        samples_df: pd.DataFrame,
        sub_samples_df: List[pd.DataFrame] = None,
        config: dict = None,
    ):
        """
        Init a peppy project instance from a pandas Dataframe

        :param samples_df: in-memory pandas DataFrame object of samples
        :param sub_samples_df: in-memory list of pandas DataFrame objects of sub-samples
        :param config: dict of yaml file
        """
        tmp_obj = cls()
        if not config:
            config = {CONFIG_VERSION_KEY: PEP_LATEST_VERSION}
        tmp_obj[SAMPLE_DF_KEY] = samples_df.replace(np.nan, "")
        tmp_obj[SUBSAMPLE_DF_KEY] = sub_samples_df

        tmp_obj[SAMPLE_DF_LARGE] = tmp_obj[SAMPLE_DF_KEY].shape[0] > 1000

        tmp_obj[CONFIG_KEY] = config

        tmp_obj.create_samples(modify=False if tmp_obj[SAMPLE_TABLE_FILE_KEY] else True)
        tmp_obj._sample_table = tmp_obj._get_table_from_samples(
            index=tmp_obj.st_index, initial=True
        )
        return tmp_obj

    @classmethod
    def from_pephub(cls, registry_path: str) -> "Project":
        """
        Init project from pephubclient.

        :param registry_path: PEPhub registry path
        :return: peppy Project
        """
        from pephubclient import PEPHubClient

        phc = PEPHubClient()
        return phc.load_project(project_registry_path=registry_path)

    @classmethod
    def from_dict(cls, pep_dictionary: dict):
        """
        Init a peppy project instance from a dictionary representation
        of an already processed PEP.

        :param Dict[Any] pep_dictionary: dict representation of the project {_config: dict,
                                                                             _samples: list | dict,
                                                                             _subsamples: list[list | dict]}
        """
        _LOGGER.info("Processing project from dictionary...")
        temp_obj = cls()
        return temp_obj._from_dict(pep_dictionary)

    def _from_dict(self, pep_dictionary) -> "Project":
        """
        Initiate a peppy project instance from a dictionary representation of an already processed PEP.

        # This function is needed in looper to reinit the project after it was created from a dictionary representation.

        :param Dict[Any] pep_dictionary: dict representation of the project {_config: dict,
                                                                             _samples: list | dict,
                                                                             _subsamples: list[list | dict]}
        """
        self[SAMPLE_DF_KEY] = pd.DataFrame(pep_dictionary[SAMPLE_RAW_DICT_KEY]).replace(
            np.nan, ""
        )
        self[CONFIG_KEY] = pep_dictionary[CONFIG_KEY]
        self[ORIGINAL_CONFIG_KEY] = deepcopy(self[CONFIG_KEY])

        if SUBSAMPLE_RAW_LIST_KEY in pep_dictionary:
            if pep_dictionary[SUBSAMPLE_RAW_LIST_KEY]:
                self[SUBSAMPLE_DF_KEY] = [
                    pd.DataFrame(sub_a).replace(np.nan, "")
                    for sub_a in pep_dictionary[SUBSAMPLE_RAW_LIST_KEY]
                ]
        if NAME_KEY in self[CONFIG_KEY]:
            self.name = self[CONFIG_KEY][NAME_KEY]

        if DESC_KEY in self[CONFIG_KEY]:
            self.description = self[CONFIG_KEY][DESC_KEY]

        self._set_indexes(self[CONFIG_KEY])

        self.create_samples(modify=False if self[SAMPLE_TABLE_FILE_KEY] else True)
        self._sample_table = self._get_table_from_samples(
            index=self.st_index, initial=True
        )

        return self

    @classmethod
    def from_pep_config(
        cls,
        cfg: str = None,
        amendments: Union[str, Iterable[str]] = None,
        sample_table_index: Union[str, Iterable[str]] = None,
        subsample_table_index: Union[str, Iterable[str]] = None,
        defer_samples_creation: bool = False,
    ):
        """
        Init a peppy project instance from a yaml file

        :param str cfg: Project config file (YAML) or sample table (CSV/TSV)
            with one row per sample to constitute project
        :param str | Iterable[str] sample_table_index: name of the columns to set
            the sample_table index to
        :param str | Iterable[str] subsample_table_index: name of the columns to set
            the subsample_table index to
        :param str | Iterable[str] amendments: names of the amendments to activate
        :param Iterable[str] amendments: amendments to use within configuration file
        :param bool defer_samples_creation: whether the sample creation should be skipped
        """
        # TODO: this is just a copy of the __init__ method. It should be refactored
        return cls(
            cfg=cfg,
            amendments=amendments,
            sample_table_index=sample_table_index,
            subsample_table_index=subsample_table_index,
            defer_samples_creation=defer_samples_creation,
        )

    @classmethod
    def from_sample_yaml(cls, yaml_file: str):
        """
        Init a peppy project instance from a yaml file

        :param str yaml_file: path to yaml file
        """
        _LOGGER.info("Processing project from yaml...")
        with open(yaml_file, "r") as f:
            prj_dict = yaml.safe_load(f)
        pd_df = pd.DataFrame.from_dict(prj_dict)
        return cls.from_pandas(pd_df)

    def to_dict(
        self,
        # expand: bool = False, # expand was used to expand paths. This functionality was removed, because of attmapp
        extended: bool = False,
        orient: Literal[
            "dict", "list", "series", "split", "tight", "records", "index"
        ] = "dict",
    ) -> dict:
        """
        Convert the Project object to a dictionary.

        :param bool extended: whether to produce complete project dict (used to reinit the project)
        :param Literal orient: orientation of the returned df
        :return dict: a dictionary representation of the Project object
        """
        if extended:
            if self[SUBSAMPLE_DF_KEY] is not None:
                sub_df = [
                    sub_a.to_dict(orient=orient) for sub_a in self[SUBSAMPLE_DF_KEY]
                ]
            else:
                sub_df = None

            if not self.get(ORIGINAL_CONFIG_KEY):
                self[ORIGINAL_CONFIG_KEY] = self[CONFIG_KEY]
            try:
                self[ORIGINAL_CONFIG_KEY][NAME_KEY] = self.name
            except NotImplementedError:
                self[ORIGINAL_CONFIG_KEY][NAME_KEY] = "unnamed"
            self[ORIGINAL_CONFIG_KEY][DESC_KEY] = self.description
            p_dict = {
                SAMPLE_RAW_DICT_KEY: self[SAMPLE_DF_KEY].to_dict(orient=orient),
                CONFIG_KEY: dict(self[ORIGINAL_CONFIG_KEY]),
                SUBSAMPLE_RAW_LIST_KEY: sub_df,
            }
        else:
            p_dict = {
                "project": self.config.copy(),
                "samples": [s.to_dict() for s in self.samples],
            }

        return p_dict

    def create_samples(self, modify: bool = False):
        """
        Populate Project with Sample objects
        """
        self._samples: List[Sample] = self.load_samples()
        if self.samples is None:
            _LOGGER.debug("No samples found in the project.")

        if modify:
            self.modify_samples()
        else:
            self._assert_samples_have_names()
            self._auto_merge_duplicated_names()

    def _reinit(self):
        """
        Clear all object attributes and initialize again
        """
        if hasattr(self, "config_file"):
            cfg_path = getattr(self, "config_file")
        else:
            cfg_path = None
        obj_attributes = self.__dict__.copy().keys()
        for attr in obj_attributes:
            delattr(self, attr)
        self.__init__(cfg=cfg_path)

    def _get_table_from_samples(self, index, initial=False):
        """
        Generate a data frame from samples. Excludes private
        attrs (prepended with an underscore)

        :param str | Iterable[str] index: name of the columns to set the index to
        :return pandas.DataFrame: a data frame with current samples attributes
        """
        if initial and not self._modifier_exists():
            # if the sample table is generated for the first time
            # (there is no chance of manual sample edits)
            # and no sample_modifiers section is defined in the config,
            # then we can simply reuse the previously read anno sheet.
            if SAMPLE_DF_KEY in self:
                df = self[SAMPLE_DF_KEY]
            else:
                df = pd.DataFrame()
        else:
            df = pd.DataFrame.from_dict([s.to_dict() for s in self.samples])
        index = [index] if isinstance(index, str) else index
        if not all([i in df.columns for i in index]):
            _LOGGER.debug(
                f"Could not set {CFG_SAMPLE_TABLE_KEY} index. At least one of the "
                f"requested columns does not exist: {index}"
            )
            return df
        _LOGGER.debug(f"Setting sample_table index to: {index}")
        df.set_index(keys=index, drop=False, inplace=True)
        return df

    def parse_config_file(
        self,
        cfg_path: str = None,
        amendments: Iterable[str] = None,
    ):
        """
        Parse provided yaml config file and check required fields exist.

        :param str cfg_path: path to the config file to read and parse
        :param Iterable[str] amendments: Name of amendments to activate
        :raises KeyError: if config file lacks required section(s)
        """
        if CONFIG_KEY not in self:
            self[CONFIG_KEY] = {}
            self[ORIGINAL_CONFIG_KEY] = {}
        if not os.path.exists(cfg_path) and not is_url(cfg_path):
            raise OSError(f"Project config file path does not exist: {cfg_path}")
        if os.path.islink(cfg_path):
            cfg_path = os.path.realpath(
                cfg_path
            )  # due to some problems with symlinks in Nextflow
        config = load_yaml(cfg_path)

        assert isinstance(
            config, Mapping
        ), "Config file parse did not yield a Mapping; got {} ({})".format(
            config, type(config)
        )

        _LOGGER.debug(f"Raw ({cfg_path}) config data: {config}")

        self._set_indexes(config)
        # recursively import configs
        if (
            PROJ_MODS_KEY in config
            and CFG_IMPORTS_KEY in config[PROJ_MODS_KEY]
            and config[PROJ_MODS_KEY][CFG_IMPORTS_KEY]
        ):
            _make_sections_absolute(config[PROJ_MODS_KEY], [CFG_IMPORTS_KEY], cfg_path)
            _LOGGER.info(
                "Importing external Project configurations: {}".format(
                    ", ".join(config[PROJ_MODS_KEY][CFG_IMPORTS_KEY])
                )
            )
            for i in config[PROJ_MODS_KEY][CFG_IMPORTS_KEY]:
                _LOGGER.debug("Processing external config: {}".format(i))
                if os.path.exists(i):
                    self.parse_config_file(cfg_path=i)
                else:
                    _LOGGER.warning(
                        "External Project configuration does not" " exist: {}".format(i)
                    )

        self[CONFIG_KEY].update(**config)
        self[ORIGINAL_CONFIG_KEY] = deepcopy(self[CONFIG_KEY])
        # Parse yaml into the project.config attributes
        _LOGGER.debug("Adding attributes: {}".format(", ".join(config)))
        # Overwrite any config entries with entries in the amendments
        amendments = [amendments] if isinstance(amendments, str) else amendments
        if amendments:
            for amendment in amendments:
                c = self[CONFIG_KEY]
                if (
                    PROJ_MODS_KEY in c
                    and AMENDMENTS_KEY in c[PROJ_MODS_KEY]
                    and c[PROJ_MODS_KEY][AMENDMENTS_KEY] is not None
                ):
                    _LOGGER.debug("Adding entries for amendment '{}'".format(amendment))
                    try:
                        amends = c[PROJ_MODS_KEY][AMENDMENTS_KEY][amendment]
                    except KeyError:
                        raise MissingAmendmentError(
                            amendment, c[PROJ_MODS_KEY][AMENDMENTS_KEY]
                        )
                    _LOGGER.debug("Updating with: {}".format(amends))
                    self[CONFIG_KEY].update(**amends)
                    _LOGGER.info("Using amendments: {}".format(amendment))
                else:
                    raise MissingAmendmentError(amendment)
            self[ACTIVE_AMENDMENTS_KEY] = amendments

        # determine config version and reformat it, if needed
        self[CONFIG_KEY][CONFIG_VERSION_KEY] = self.pep_version
        # here specify cfg sections that may need expansion
        relative_vars = [CFG_SAMPLE_TABLE_KEY, CFG_SUBSAMPLE_TABLE_KEY]
        _make_sections_absolute(self[CONFIG_KEY], relative_vars, cfg_path)

    def _set_indexes(self, config: Mapping) -> None:
        """
        Set sample and subsample indexes if they are different then Default

        :param config: project config
        """
        self.st_index = (
            config[SAMPLE_TABLE_INDEX_KEY]
            if SAMPLE_TABLE_INDEX_KEY in config
            else SAMPLE_NAME_ATTR
        )
        self.sst_index = (
            config[SUBSAMPLE_TABLE_INDEX_KEY]
            if SUBSAMPLE_TABLE_INDEX_KEY in config
            else SUBSAMPLE_NAME_ATTR
        )
        return None

    def load_samples(self):
        """
        Read the sample_table and subsample_tables into dataframes
        and store in the object root. The values sourced from the
        project config can be overwritten by the optional arguments.
        """
        # To initiate project from pandas or dictionary we shouldn't run
        # this function otherwise it will cause errors
        if SAMPLE_DF_KEY not in self or self.amendments is not None:
            self._read_sample_data()

        samples_list = []
        if SAMPLE_DF_KEY not in self:
            return []

        if CONFIG_KEY not in self:
            self[CONFIG_KEY] = {CONFIG_VERSION_KEY: PEP_LATEST_VERSION}
            self[CONFIG_FILE_KEY] = None

        elif len(self[CONFIG_KEY]) < 1:
            self[CONFIG_KEY][CONFIG_VERSION_KEY] = PEP_LATEST_VERSION
            self[CONFIG_FILE_KEY] = None

        if SUBSAMPLE_DF_KEY not in self:
            self[SUBSAMPLE_DF_KEY] = None

        for _, r in self[SAMPLE_DF_KEY].iterrows():
            samples_list.append(Sample(r, prj=self))
        return samples_list

    def modify_samples(self):
        """
        Perform any sample modifications defined in the config.
        """
        if self._modifier_exists():
            # check for unrecognizable modification keys
            mod_diff = set(self[CONFIG_KEY][SAMPLE_MODS_KEY].keys()) - set(
                SAMPLE_MODIFIERS
            )
            if len(mod_diff) > 0:
                _LOGGER.warning(
                    f"Config '{SAMPLE_MODS_KEY}' section contains unrecognized "
                    f"subsections: {mod_diff}"
                )
        self.attr_remove()
        self.attr_constants()
        self.attr_synonyms()
        self.attr_imply()
        self._assert_samples_have_names()
        self._auto_merge_duplicated_names()
        self.attr_merge()
        self.attr_derive()

    def _modifier_exists(self, modifier_key=None):
        """
        Check whether a specified sample modifier is defined and can be applied

        If no modifier is specified, only the sample_modifiers section's
        existence is checked

        :param str modifier_key: modifier key to be checked
        :return bool: whether the requirements are met
        """
        _LOGGER.debug("Checking existence: {}".format(modifier_key))
        if CONFIG_KEY not in self or SAMPLE_MODS_KEY not in self[CONFIG_KEY]:
            return False
        if (
            modifier_key is not None
            and modifier_key not in self[CONFIG_KEY][SAMPLE_MODS_KEY]
        ):
            return False
        return True

    def attr_remove(self):
        """
        Remove declared attributes from all samples that have them defined
        """

        def _del_if_in(obj, attr):
            if attr in obj:
                del obj[attr]

        if self._modifier_exists(REMOVE_KEY):
            to_remove = self[CONFIG_KEY][SAMPLE_MODS_KEY][REMOVE_KEY]
            _LOGGER.debug(f"Removing attributes: {to_remove}")
            for s in track(
                self.samples,
                description="Removing sample attributes",
                disable=not (self.is_sample_table_large and self.progressbar),
                console=Console(file=sys.stderr),
            ):
                for attr in to_remove:
                    _del_if_in(s, attr)

    def attr_constants(self):
        """
        Update each Sample with constants declared by a Project.
        If Project does not declare constants, no update occurs.
        """
        if self._modifier_exists(APPEND_KEY):
            to_append = self[CONFIG_KEY][SAMPLE_MODS_KEY][APPEND_KEY]
            _LOGGER.debug("Applying constant attributes: {}".format(to_append))

            for s in track(
                self.samples,
                description="Applying constant sample attributes",
                disable=not (self.is_sample_table_large and self.progressbar),
                console=Console(file=sys.stderr),
            ):
                for attr, val in to_append.items():
                    if attr not in s:
                        s.update({attr: val})

    def attr_synonyms(self):
        """
        Copy attribute values for all samples to a new one
        """
        if self._modifier_exists(DUPLICATED_KEY):
            synonyms = self[CONFIG_KEY][SAMPLE_MODS_KEY][DUPLICATED_KEY]
            _LOGGER.debug(f"Applying synonyms: {synonyms}")
            for sample in track(
                self.samples,
                description="Applying synonymous sample attributes",
                disable=not (self.is_sample_table_large and self.progressbar),
                console=Console(file=sys.stderr),
            ):
                for attr, new in synonyms.items():
                    if attr in sample:
                        sample[new] = sample[attr]
                    else:
                        _LOGGER.warning(
                            f"The sample attribute to duplicate not found: {attr}"
                        )

    def _assert_samples_have_names(self):
        """
        Make sure samples have sample_name attribute specified.
        Try to derive this attribute first.

        :raise InvalidSampleTableFileException: if names are not specified
        """
        with suppress(KeyError):
            # before merging, which requires sample_name attribute to map
            # sample_table rows to subsample_table rows,
            # perform only sample_name attr derivation
            if (
                SAMPLE_NAME_ATTR
                in self[CONFIG_KEY][SAMPLE_MODS_KEY][DERIVED_KEY][DERIVED_ATTRS_KEY]
            ):
                self.attr_derive(attrs=[SAMPLE_NAME_ATTR])

        for sample in self.samples:
            if self.st_index not in sample:
                message = (
                    f"{CFG_SAMPLE_TABLE_KEY} is missing '{self.st_index}' column; "
                    f"you must specify {CFG_SAMPLE_TABLE_KEY}s in {self.st_index} or derive them"
                )
                raise InvalidSampleTableFileException(message)

    def _auto_merge_duplicated_names(self):
        """
        If sample_table specifies samples with non-unique names, try to merge these samples

        :raises IllegalStateException: if both duplicated samples are detected and subsample_table is
            specified in the config
        """
        sample_names_list = [s[self.st_index] for s in self.samples]
        duplicated_sample_ids = self._get_duplicated_sample_ids(sample_names_list)

        if not duplicated_sample_ids:
            return

        _LOGGER.info(
            f"Found {len(duplicated_sample_ids)} samples with non-unique names: {duplicated_sample_ids}. "
            f"Attempting to auto-merge."
        )
        if SUBSAMPLE_DF_KEY in self and self[SUBSAMPLE_DF_KEY] is not None:
            raise IllegalStateException(
                f"Duplicated sample names found and subsample_table is specified in the config; "
                f"you may use either auto-merging or subsample_table-based merging. "
                f"Duplicates: {duplicated_sample_ids}"
            )

        for duplicated_id in duplicated_sample_ids:
            (
                duplicated_samples,
                non_duplicated_samples,
            ) = self._get_duplicated_and_not_duplicated_samples(
                duplicated_id, self.st_index, self.samples
            )
            sample_attributes = [
                attr
                for attr in duplicated_samples[0].keys()
                if not attr.startswith("_")
            ]
            merged_attrs = self._get_merged_attributes(
                sample_attributes, duplicated_samples
            )
            self._samples = non_duplicated_samples

            # make single element lists scalars
            for attribute_name, values in merged_attrs.items():
                if isinstance(values, list) and len(list(set(values))) == 1:
                    merged_attrs[attribute_name] = values[0]

            self.add_samples(Sample(series=merged_attrs))

    @staticmethod
    def _get_duplicated_and_not_duplicated_samples(
        duplication: str, st_index: str, samples: List[Sample]
    ) -> Tuple[List, List]:
        """
        Iterates over list of samples and splits them into list of duplicated samples and list of not duplicated.

        Args:
            duplication: Sample ID.
            st_index: Sample table index.
            samples: List of Sample instances.

        Returns:
            List of Sample objects of duplicated samples and not duplicated samples.
        """
        duplicated_samples, not_duplicated_samples = [], []
        for sample in samples:
            if sample[st_index] == duplication:
                duplicated_samples.append(sample)
            else:
                not_duplicated_samples.append(sample)
        return duplicated_samples, not_duplicated_samples

    @staticmethod
    def _all_values_in_the_list_are_the_same(list_of_values: List) -> bool:
        return all(value == list_of_values[0] for value in list_of_values)

    def _get_duplicated_sample_ids(self, sample_names_list: List) -> set:
        duplicated_ids = set()
        seen = set()

        for sample_id in sample_names_list:
            if sample_id in seen:
                duplicated_ids.add(sample_id)
            else:
                seen.add(sample_id)

        return duplicated_ids

    @staticmethod
    def _get_merged_attributes(
        sample_attributes: List[str], duplicated_samples: List[Sample]
    ) -> dict:
        merged_attributes = {}
        for attr in sample_attributes:
            attribute_values = []
            for sample in duplicated_samples:
                attribute_value_for_sample = sample.get(attr, "")
                attribute_values.append(attribute_value_for_sample)

            merged_attributes[attr] = list(flatten(attribute_values))

        return merged_attributes

    def attr_merge(self):
        """
        Merge sample subannotations (from subsample table) with
        sample annotations (from sample_table)
        """
        if SUBSAMPLE_DF_KEY not in self or self[SUBSAMPLE_DF_KEY] is None:
            _LOGGER.debug("No {} found, skipping merge".format(CFG_SUBSAMPLE_TABLE_KEY))
            return
        for subsample_table in self[SUBSAMPLE_DF_KEY]:
            for sample_name in list(subsample_table[self.st_index]):
                if sample_name not in [s[self.st_index] for s in self.samples]:
                    _LOGGER.warning(
                        ("Couldn't find matching sample for subsample: {}").format(
                            sample_name
                        )
                    )
            for sample in track(
                self.samples,
                description=f"Merging subsamples, adding sample attrs: {', '.join(subsample_table.keys())}",
                disable=not (self.is_sample_table_large and self.progressbar),
                console=Console(file=sys.stderr),
            ):
                sample_colname = self.st_index
                if sample_colname not in subsample_table.columns:
                    raise KeyError(
                        "Subannotation requires column '{}'.".format(sample_colname)
                    )
                _LOGGER.debug(
                    "Using '{}' as sample name column from "
                    "subannotation table".format(sample_colname)
                )
                sample_indexer = (
                    subsample_table[sample_colname] == sample[self.st_index]
                )
                this_sample_rows = subsample_table[sample_indexer]
                if len(this_sample_rows) == 0:
                    _LOGGER.debug(
                        "No merge rows for sample '%s', skipping",
                        sample[self.st_index],
                    )
                    continue
                _LOGGER.debug("%d rows to merge", len(this_sample_rows))
                _LOGGER.debug("Merge rows dict: {}".format(this_sample_rows.to_dict()))

                merged_attrs = {key: list() for key in this_sample_rows.columns}
                _LOGGER.debug(this_sample_rows)
                for subsample_row_id, row in this_sample_rows.iterrows():
                    try:
                        row[SUBSAMPLE_NAME_ATTR]
                    except KeyError:
                        row[SUBSAMPLE_NAME_ATTR] = str(subsample_row_id)
                    rowdata = row.to_dict()

                    def _select_new_attval(merged_attrs, attname, attval):
                        """
                        Select new attribute value for the merged columns
                        dictionary
                        """
                        if attname in merged_attrs:
                            return merged_attrs[attname] + [attval]
                        return [str(attval).rstrip()]

                    for attname, attval in rowdata.items():
                        if attname == sample_colname or not attval:
                            _LOGGER.debug(f"Skipping KV: {attname}={attval}")
                            continue
                        _LOGGER.debug(
                            f"merge: sample '{sample[self.st_index]}'; '{attname}'='{attval}'"
                        )
                        merged_attrs[attname] = _select_new_attval(
                            merged_attrs, attname, attval
                        )

                # remove sample name from the data with which to update sample
                merged_attrs.pop(sample_colname, None)

                _LOGGER.debug(
                    f"Updating Sample {sample[self.st_index]}: {merged_attrs}"
                )
                sample.update(merged_attrs)

    def attr_imply(self):
        """
        Infer value for additional field(s) from other field(s).

        Add columns/fields to the sample based on values in those already-set
        that the sample's project defines as indicative of implications for
        additional data elements for the sample.
        """
        if not self._modifier_exists(IMPLIED_KEY):
            return
        implications = self[CONFIG_KEY][SAMPLE_MODS_KEY][IMPLIED_KEY]
        if not isinstance(implications, list):
            raise InvalidConfigFileException(
                f"{SAMPLE_MODS_KEY}.{IMPLIED_KEY} has to be a list of key-value pairs"
            )
        _LOGGER.debug(f"Sample attribute implications: {implications}")
        for implication in implications:
            if not all([key in implication for key in IMPLIED_COND_KEYS]):
                raise InvalidConfigFileException(
                    f"{SAMPLE_MODS_KEY}.{IMPLIED_KEY} section is invalid: {implication}"
                )
        for sample in track(
            self.samples,
            description="Implying sample attributes",
            disable=not (self.is_sample_table_large and self.progressbar),
            console=Console(file=sys.stderr),
        ):
            for implication in implications:
                implier_attrs = list(implication[IMPLIED_IF_KEY].keys())
                implied_attrs = list(implication[IMPLIED_THEN_KEY].keys())
                _LOGGER.debug(f"Setting Sample attributes implied by '{implier_attrs}'")
                for implier_attr in implier_attrs:
                    implier_val = implication[IMPLIED_IF_KEY][implier_attr]
                    if implier_attr not in sample:
                        _LOGGER.debug(
                            f"Sample lacks implier attr ({implier_attr}), skipping:"
                        )
                        break
                    sample_val = sample[implier_attr]
                    if sample_val not in implier_val:
                        _LOGGER.debug(
                            "Sample attr value does not match any of implier "
                            f"requirements ({sample_val} not in {implier_val}), skipping"
                        )
                        break
                else:
                    # only executed if the inner loop did NOT break
                    for implied_attr in implied_attrs:
                        imp_val = implication[IMPLIED_THEN_KEY][implied_attr]
                        _LOGGER.debug(
                            f"Setting implied attr: '{implied_attr}={imp_val}'"
                        )
                        sample.__setitem__(implied_attr, imp_val)

    def attr_derive(self, attrs=None):
        """
        Set derived attributes for all Samples tied to this Project instance
        """
        if not self._modifier_exists(DERIVED_KEY):
            return
        da = self[CONFIG_KEY][SAMPLE_MODS_KEY][DERIVED_KEY][DERIVED_ATTRS_KEY]
        ds = self[CONFIG_KEY][SAMPLE_MODS_KEY][DERIVED_KEY][DERIVED_SOURCES_KEY]
        derivations = attrs or (da if isinstance(da, list) else [da])
        _LOGGER.debug("Derivations to be done: {}".format(derivations))
        for sample in track(
            self.samples,
            description="Deriving sample attributes",
            disable=not (self.is_sample_table_large and self.progressbar),
            console=Console(file=sys.stderr),
        ):
            for attr in derivations:
                if attr not in sample:
                    _LOGGER.debug(f"sample lacks '{attr}' attribute")
                    continue
                elif attr in sample._derived_cols_done:
                    _LOGGER.debug(f"'{attr}' has been derived")
                    continue
                _LOGGER.debug(
                    f"Deriving '{attr}' attribute for '{sample[self.st_index]}'"
                )

                # Set {atr}_key, so the original source can also be retrieved
                sample[ATTR_KEY_PREFIX + attr] = sample[attr]

                derived_attr = sample.derive_attribute(ds, attr)
                if derived_attr:
                    _LOGGER.debug("Setting '{}' to '{}'".format(attr, derived_attr))
                    sample[attr] = derived_attr
                else:
                    _LOGGER.debug(
                        f"Not setting null/empty value for data source '{attr}': {type(derived_attr)}"
                    )
                sample._derived_cols_done.append(attr)

    def activate_amendments(self, amendments):
        """
        Update settings based on amendment-specific values.

        This method will update Project attributes, adding new values
        associated with the amendments indicated, and in case of collision with
        an existing key/attribute the amendments' values will be favored.

        :param Iterable[str] amendments: A string with amendment
            names to be activated
        :return peppy.Project: Updated Project instance
        :raise TypeError: if argument to amendment parameter is null
        :raise NotImplementedError: if this call is made on a project not
            created from a config file
        """
        amendments = [amendments] if isinstance(amendments, str) else amendments
        if amendments is None:
            raise TypeError(
                "The amendment argument can not be null. To deactivate an "
                "amendment use the deactivate_amendments method."
            )
        if not self[CONFIG_FILE_KEY]:
            raise NotImplementedError(
                "amendment activation isn't supported on a project not "
                "created from a config file"
            )
        prev = [(k, v) for k, v in self.items() if not k.startswith("_")]
        conf_file = self[CONFIG_FILE_KEY]
        self.__init__(cfg=conf_file, amendments=amendments)
        for k, v in prev:
            if k.startswith("_"):
                continue
            if k not in self or (self.is_null(k) and v is not None):
                _LOGGER.debug("Restoring {}: {}".format(k, v))
                self[k] = v
        self[ACTIVE_AMENDMENTS_KEY] = amendments
        return self

    def deactivate_amendments(self):
        """
        Bring the original project settings back.

        :return peppy.Project: Updated Project instance
        :raise NotImplementedError: if this call is made on a project not
            created from a config file
        """
        if ACTIVE_AMENDMENTS_KEY not in self or self[ACTIVE_AMENDMENTS_KEY] is None:
            _LOGGER.warning("No amendments have been activated.")
            return self
        if not self[CONFIG_FILE_KEY]:
            raise NotImplementedError(
                "amendments deactivation isn't supported on a project that "
                "lacks a config file."
            )
        self._reinit()
        return self

    def add_samples(self, samples):
        """
        Add list of Sample objects

        :param peppy.Sample | Iterable[peppy.Sample] samples: samples to add
        """
        samples = [samples] if isinstance(samples, Sample) else samples
        for sample in samples:
            if not isinstance(sample, Sample):
                _LOGGER.warning("Not a peppy.Sample object, not adding")
                continue
            self._samples.append(sample)
            self[SAMPLE_EDIT_FLAG_KEY] = True

    def remove_samples(self, sample_names):
        """
        Remove Samples from Project

        :param Iterable[str] sample_names: sample names to remove
        """
        sample_names = [sample_names] if isinstance(sample_names, str) else sample_names
        samples_keep = [
            s for s in self.samples if s[self.sample_name_colname] not in sample_names
        ]
        if len(self._samples) != len(samples_keep):
            self._samples = samples_keep
            self[SAMPLE_EDIT_FLAG_KEY] = True

    def infer_name(self):
        """
        Infer project name from config file path.

        First assume the name is the folder in which the config file resides,
        unless that folder is named "metadata", in which case the project name
        is the parent of that folder.

        :return str: inferred name for project.
        :raise InvalidConfigFileException: if the project lacks both a name and
            a configuration file (no basis, then, for inference)
        :raise InvalidConfigFileException: if specified Project name is invalid
        """
        if CONFIG_KEY not in self:
            return
        if NAME_KEY in self[CONFIG_KEY]:
            if " " in self[CONFIG_KEY][NAME_KEY]:
                raise InvalidConfigFileException(
                    "Specified Project name ({}) contains whitespace".format(
                        self[CONFIG_KEY][NAME_KEY]
                    )
                )
            return self[CONFIG_KEY][NAME_KEY].replace(" ", "_")
        if not self[CONFIG_FILE_KEY]:
            raise NotImplementedError(
                "Project name inference isn't supported "
                "on a project that lacks a config file."
            )
        config_folder = os.path.dirname(self[CONFIG_FILE_KEY])
        project_name = os.path.basename(config_folder)
        if project_name == METADATA_KEY:
            project_name = os.path.basename(os.path.dirname(config_folder))
        return project_name.replace(" ", "_")

    def get_description(self):
        """
        Infer project description from config file.

        The provided description has to be of class coercible to string

        :return str: inferred name for project.
        :raise InvalidConfigFileException: if description is not of class
            coercible to string
        """
        if CONFIG_KEY not in self:
            return
        if DESC_KEY in self[CONFIG_KEY]:
            desc_str = str(self[CONFIG_KEY][DESC_KEY])
            if not isinstance(desc_str, str):
                try:
                    desc_str = str(desc_str)
                except Exception as e:
                    raise InvalidConfigFileException(
                        "Could not convert the specified Project description "
                        "({}) to string. Caught exception: {}".format(
                            desc_str, getattr(e, "message", repr(e))
                        )
                    )
            return desc_str

    def __str__(self):
        """Representation in interpreter."""
        if len(self) == 0:
            return "{}"
        msg = "Project"
        if NAME_KEY in self and self[NAME_KEY] is not None:
            msg += f" '{self[NAME_KEY]}'"
        if CONFIG_FILE_KEY in self and self[CONFIG_FILE_KEY] is not None:
            msg += f" ({self[CONFIG_FILE_KEY]})"
        if DESC_KEY in self and self[DESC_KEY] is not None:
            msg += f"\n{DESC_KEY}: {self[DESC_KEY]}"
        try:
            num_samples = len(self._samples)
        except (AttributeError, TypeError):
            _LOGGER.debug("No samples established on project")
            num_samples = 0
        if num_samples > 0:
            msg = f"{msg}\n{num_samples} samples"
            sample_names = [s[self.st_index] for s in self.samples]
            repr_names = sample_names[:MAX_PROJECT_SAMPLES_REPR]
            context = (
                f" (showing first {MAX_PROJECT_SAMPLES_REPR})"
                if num_samples > MAX_PROJECT_SAMPLES_REPR
                else ""
            )
            msg = f"{msg}{context}: {', '.join(repr_names)}"
        else:
            msg = f"{msg} 0 samples"
        if CONFIG_KEY not in self:
            return msg
        msg = f"{msg}\nSections: {', '.join([s for s in self[CONFIG_KEY].keys()])}"
        if (
            PROJ_MODS_KEY in self[CONFIG_KEY]
            and AMENDMENTS_KEY in self[CONFIG_KEY][PROJ_MODS_KEY]
        ):
            msg = f"{msg}\nAmendments: {', '.join(self[CONFIG_KEY][PROJ_MODS_KEY][AMENDMENTS_KEY].keys())}"
        if self.amendments:
            msg = (
                f"{msg}\nActivated amendments: {', '.join(self[ACTIVE_AMENDMENTS_KEY])}"
            )
        return msg

    @property
    def amendments(self):
        """
        Return currently active list of amendments or None if none was activated

        :return Iterable[str]: a list of currently active amendment names
        """
        return self[ACTIVE_AMENDMENTS_KEY] if ACTIVE_AMENDMENTS_KEY in self else None

    @property
    def list_amendments(self):
        """
        Return a list of available amendments or None if not declared

        :return Iterable[str]: a list of available amendment names
        """
        try:
            return list(self[CONFIG_KEY][PROJ_MODS_KEY][AMENDMENTS_KEY].keys())
        except Exception as e:
            _LOGGER.debug(
                "Could not retrieve available amendments: {}".format(
                    getattr(e, "message", repr(e))
                )
            )
            return None

    @property
    def config(self):
        """
        Get the config mapping

        :return Mapping: config. May be formatted to comply with the most
            recent version specifications
        """
        return self[CONFIG_KEY] if CONFIG_KEY in self else {}

    @property
    def config_file(self):
        """
        Get the config file path

        :return str: path to the config file
        """
        return self[CONFIG_FILE_KEY]

    @property
    def samples(self):
        """
        Generic/base Sample instance for each of this Project's samples.

        :return Iterable[Sample]: Sample instance for each
            of this Project's samples
        """
        if self._samples is not None:
            return self._samples
        if SAMPLE_DF_KEY not in self or self[SAMPLE_DF_KEY] is None:
            _LOGGER.debug("No samples are defined")
            return []

    @property
    def sample_table_index(self):
        """
        The effective sample table index.

        It is `sample_name` by default, but could be overwritten by the selected sample table index,
        defined on the object instantiation stage or in the project configuration file
        via `sample_table_index` field.

        That's the sample table index selection priority order:

        1. Constructor specified
        2. Config specified
        3. Deafult: `sample_table`

        :return str: name of the column that consist of sample identifiers
        """
        # this property is used solely for documentation purposes
        return self.st_index

    @property
    def subsample_table_index(self):
        """
        The effective subsample table indexes.

        It is `[subasample_name, sample_name]` by default,
        but could be overwritten by the selected subsample table indexes,
        defined on the object instantiation stage or in the project configuration file
        via `subsample_table_index` field.

        That's the subsample table indexes selection priority order:

        1. Constructor specified
        2. Config specified
        3. Deafult: `[subasample_name, sample_name]`

        :return List[str]: names of the columns that consist of sample and subsample identifiers
        """
        # this property is used solely for documentation purposes
        return self.sst_index

    @property
    def sample_name_colname(self):
        """
        **Deprecated, please use `Project.sample_table_index` instead**

        Name of the effective sample name containing column in the sample table.

        It is "sample_name" by default, but when it's missing it could be
        replaced by the selected sample table index, defined on the
        object instantiation stage.

        :return str: name of the column that consist of sample identifiers
        """
        return SAMPLE_NAME_ATTR if SAMPLE_NAME_ATTR == self.st_index else self.st_index

    @property
    def sample_table(self):
        """
        Get sample table. If any sample edits were performed,
        it will be re-generated

        :return pandas.DataFrame: a data frame with current samples attributes
        """
        if self[SAMPLE_EDIT_FLAG_KEY]:
            _LOGGER.debug("Generating new sample_table DataFrame")
            self[SAMPLE_EDIT_FLAG_KEY] = False
            new_df = self._get_table_from_samples(index=self.st_index)
            self._sample_table = new_df
            return new_df

        _LOGGER.debug("Returning stashed sample_table DataFrame")
        return self._sample_table

    @property
    def subsample_table(self):
        """
        Get subsample table

        :return pandas.DataFrame: a data frame with subsample attributes
        """

        if not self[SUBSAMPLE_DF_KEY]:
            return None

        subsample_dataframes_array = make_list(self[SUBSAMPLE_DF_KEY], pd.DataFrame)
        for subsample_table in subsample_dataframes_array:
            if not all([i in subsample_table.columns for i in self.sst_index]):
                _LOGGER.info(
                    "Could not set {} index. At least one of the"
                    " requested columns does not exist: {}".format(
                        CFG_SUBSAMPLE_TABLE_KEY, self.sst_index
                    )
                )
                return subsample_table
            subsample_table.set_index(keys=self.sst_index, drop=False, inplace=True)
            _LOGGER.info("Setting subsample_table index to: {}".format(self.sst_index))
            subsample_table.index = self._fix_subsample_table_index(subsample_table)
        return (
            subsample_dataframes_array
            if len(subsample_dataframes_array) > 1
            else subsample_dataframes_array[0]
        )

    @staticmethod
    def _fix_subsample_table_index(subsample_table):
        new_levels = []
        current_levels = getattr(subsample_table.index, "levels", [])

        for index in current_levels:
            new_levels.append(index.astype(str))

        if current_levels:
            subsample_table.index = subsample_table.index.set_levels(new_levels)

        return subsample_table.index

    @property
    def is_sample_table_large(self):
        return getattr(self, SAMPLE_DF_LARGE, False)

    def _read_sample_data(self):
        """
        Read the sample_table and subsample_table into dataframes
        and store in the object root. The values sourced from the
        project config can be overwritten by the optional arguments

        :param str sample_table: a path to a sample table
        :param List[str] sample_table: a list of paths to sample tables
        """

        no_metadata_msg = "No {} specified"
        if self[SAMPLE_TABLE_FILE_KEY] is not None:
            st = self[SAMPLE_TABLE_FILE_KEY]
        else:
            if CONFIG_KEY not in self:
                _LOGGER.info("No config key in Project, or reading project from dict")
                return
            if CFG_SAMPLE_TABLE_KEY not in self[CONFIG_KEY]:
                _LOGGER.debug(f"No {CFG_SAMPLE_TABLE_KEY} found in config file")
                return
            st = self[CONFIG_KEY][CFG_SAMPLE_TABLE_KEY]

        if self[SUBSAMPLE_TABLES_FILE_KEY] is not None:
            sst = self[SUBSAMPLE_TABLES_FILE_KEY]
        else:
            if (
                CONFIG_KEY in self
                and self[CONFIG_KEY].get(CFG_SUBSAMPLE_TABLE_KEY) is not None
            ):
                sst = make_list(self[CONFIG_KEY][CFG_SUBSAMPLE_TABLE_KEY], str)
            else:
                sst = None

        if st is not None:
            parser_class = select_parser(path=st)
            self[SAMPLE_DF_KEY] = parser_class(path=st).table.replace(np.nan, "")
            self[SAMPLE_DF_LARGE] = self[SAMPLE_DF_KEY].shape[0] > 1000
        else:
            _LOGGER.warning(no_metadata_msg.format(CFG_SAMPLE_TABLE_KEY))
            self[SAMPLE_DF_KEY] = None
        if sst is not None:
            ssts = []
            for subsample_table in sst:
                parser_class = select_parser(path=subsample_table)
                ssts.append(parser_class(path=subsample_table).table)
            self[SUBSAMPLE_DF_KEY] = ssts
        else:
            _LOGGER.debug(no_metadata_msg.format(CFG_SUBSAMPLE_TABLE_KEY))
            self[SUBSAMPLE_DF_KEY] = None

    @property
    def pep_version(self):
        """
        The declared PEP version string

        It is validated to make sure it is a valid PEP version string

        :raise InvalidConfigFileException: in case of invalid PEP version
        :return str: PEP version string
        """
        req_version_str = ".".join(REQUIRED_VERSION)
        # if no config file was passed (init with sample table csv)
        # then just auto-assign the required version
        if CONFIG_KEY not in self:
            return req_version_str
        # if the pep_version was omitted from the config file
        # then just warn the user and auto assign
        if CONFIG_VERSION_KEY not in self[CONFIG_KEY]:
            _LOGGER.warning(
                f"Config file does not have version key. Defaulting to {req_version_str}"
            )
            return req_version_str

        v_str = self[CONFIG_KEY][CONFIG_VERSION_KEY]
        if not isinstance(v_str, str):
            raise InvalidConfigFileException(f"{CONFIG_VERSION_KEY} must be a string")
        v_bundle = v_str.split(".")
        if len(v_bundle) != 3:
            raise InvalidConfigFileException(
                f"Version string is not tripartite: {v_str}"
            )
        try:
            v_bundle = list(map(int, v_bundle))
        except ValueError:
            raise InvalidConfigFileException(
                f"Version string elements are not coercible to integers: {v_str}"
            )
        if v_bundle[0] < 2:
            raise InvalidConfigFileException(
                f"PEP version is invalid: {v_str}. "
                f"Please use version {req_version_str} or previous versions of peppy package."
            )
        return ".".join(list(map(str, v_bundle)))

    def get_sample(self, sample_name):
        """
        Get an individual sample object from the project.

        Will raise a ValueError if the sample is not found.
        In the case of multiple samples with the same name (which is not
        typically allowed), a warning is raised and the first sample is returned

        :param str sample_name: The name of a sample to retrieve
        :raise ValueError: if there's no sample with the specified name defined
        :return peppy.Sample: The requested Sample object
        """
        samples = self.get_samples([sample_name])
        if len(samples) > 1:
            _LOGGER.warning(
                f"{len(samples)} samples matched the name: {sample_name}. Returning the first one."
            )
        try:
            return samples[0]
        except IndexError:
            raise ValueError(f"Project has no sample named {sample_name}.")

    def get_samples(self, sample_names):
        """
        Returns a list of sample objects given a list of sample names

        :param list sample_names: A list of sample names to retrieve
        :return list[peppy.Sample]: A list of Sample objects
        """
        return [s for s in self.samples if s[self.st_index] in sample_names]

    @property
    def description(self):
        return self.get_description()

    @description.setter
    def description(self, value):
        self[CONFIG_KEY][DESC_KEY] = str(value)

    @property
    def name(self):
        return self.infer_name()

    @name.setter
    def name(self, value):
        self[CONFIG_KEY][NAME_KEY] = str(value)

    def __setitem__(self, key, value):
        self._project_data[key] = value

    def __getitem__(self, item):
        """
        Fetch the value of given key.

        :param hashable item: key for which to fetch value
        :return object: value mapped to given key, if available
        :raise KeyError: if the requested key is unmapped.
        """
        return self._project_data[item]

    def __iter__(self):
        return iter(self._project_data)

    def __len__(self):
        return len(self._project_data)

    def __delitem__(self, key):
        value = self[key]
        del self._project_data[key]
        self.pop(value, None)

    def __repr__(self):
        return str(self)

    def __reduce__(self):
        return (
            self.__class__.from_dict,
            (self.to_dict(extended=True, orient="records"),),
        )


def infer_delimiter(filepath):
    """
    From extension infer delimiter used in a separated values file.

    :param str filepath: path to file about which to make inference
    :return str | NoneType: extension if inference succeeded; else null
    """
    ext = os.path.splitext(filepath)[1][1:].lower()
    return {"txt": "\t", "tsv": "\t", "csv": ","}.get(ext)


def _make_sections_absolute(object, sections, cfg_path):
    for key in sections:
        try:
            relpath = object[key]
        except KeyError:
            _LOGGER.debug(
                "No '{}' section in configuration file: {}".format(key, cfg_path)
            )
            continue
        if relpath is None:
            continue
        _LOGGER.debug("Ensuring absolute path for '{}'".format(relpath))
        # Parsed from YAML, so small space of possible datatypes
        if isinstance(relpath, list):
            absolute = [
                make_abs_via_cfg(maybe_relpath, cfg_path) for maybe_relpath in relpath
            ]
        else:
            absolute = make_abs_via_cfg(relpath, cfg_path)
        _LOGGER.debug("Setting '{}' to '{}'".format(key, absolute))
        object[key] = absolute
