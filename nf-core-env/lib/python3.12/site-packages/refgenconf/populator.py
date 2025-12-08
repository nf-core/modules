# refgenie looper plugin

import logging
from collections.abc import Mapping

from attmap import AttMap
from ubiquerg import parse_registry_path as prp

import refgenconf

_LOGGER = logging.getLogger(__name__)


def looper_refgenie_populate(namespaces):
    """
    A looper plugin that populates refgenie references in a PEP from
    refgenie://genome/asset:tag registry paths. This can be used to convert
    all refgenie references into their local paths at the looper stage, so the
    final paths are passed to the workflow. This way the workflow does not
    need to depend on refgenie to resolve the paths.
    This is useful for example for CWL pipelines, which are built to have
    paths resolved outside the workflow.

    The namespaces structure required to run the plugin is:
    `namespaces["pipeline"]["var_templates"]["refgenie_config"]`

    :param Mapping namespaces: a nested variable namespaces dict
    :return dict: sample namespace dict
    :raises TypeError: if the input namespaces is not a mapping
    :raises KeyError: if the namespaces mapping does not include 'pipeline'
    :raises NotImplementedError: if 'var_templates' key is missing in the 'pipeline' namespace or
        'refgenie_config' is missing in 'var_templates' section.
    """
    if not isinstance(namespaces, Mapping):
        raise TypeError("Namespaces must be a Mapping")
    if "pipeline" not in namespaces:
        raise KeyError(
            "Namespaces do not include 'pipeline'. The job is misconfigured."
        )
    if (
        "var_templates" in namespaces["pipeline"]
        and "refgenie_config" in namespaces["pipeline"]["var_templates"]
    ):
        rgc_path = namespaces["pipeline"]["var_templates"]["refgenie_config"]
        rgc = refgenconf.RefGenConf(rgc_path)

        complete_sk_dict = rgc.list_seek_keys_values()
        paths_dict = {}

        # This function allows you to specify tags for specific assets to use
        # in the project config like:
        # refgenie_asset_tags:
        #   genome:
        #     asset_name: tag_name
        def get_asset_tag(genome, asset):
            try:
                return namespaces["project"]["refgenie"]["tag_overrides"][genome][asset]
            except KeyError:
                default_tag = rgc.get_default_tag(genome=genome, asset=asset)
                _LOGGER.info(
                    f"Refgenie asset ({genome}/{asset}) tag not specified in `refgenie.tag_overrides` section. "
                    f"Using the default tag: {default_tag}"
                )
                return default_tag
            except TypeError:
                default_tag = rgc.get_default_tag(genome=genome, asset=asset)
                _LOGGER.warn(f"tag_overrides section is malformed. Using default.")
                return default_tag

        # Restructure the seek key paths to make them accessible with
        # {refgenie.asset_name.seek_key} in command templates
        for g, gdict in complete_sk_dict.items():
            _LOGGER.debug(f"Processing genome {g}")
            paths_dict[g] = {}
            for k, v in gdict.items():
                tag = get_asset_tag(genome=g, asset=k)
                # print(k,v)
                try:
                    paths_dict[g][k] = v[tag]
                except KeyError:
                    _LOGGER.warn(
                        f"Can't find tag '{tag}' for asset '{g}/{k}', as specified in your project config. Using default."
                    )
                    paths_dict[g][k] = v[rgc.get_default_tag(genome=g, asset=k)]

        if "project" in namespaces and "refgenie" in namespaces["project"]:
            try:
                for po in namespaces["project"]["refgenie"]["path_overrides"]:
                    rp = prp(po["registry_path"])
                    _LOGGER.debug(
                        f"Overriding {po['registry_path']} with {po['value']}."
                    )
                    if not rp["subitem"]:
                        rp["subitem"] = rp["item"]
                    _LOGGER.debug(rp)
                    paths_dict[rp["namespace"]][rp["item"]][rp["subitem"]] = po["value"]
            except KeyError:
                _LOGGER.debug("Did not find path_overrides section")
            except TypeError:
                _LOGGER.warn("Warning: path_overrides is not iterable")

        # print(paths_dict)
        # Provide these values under the 'refgenie' namespace
        namespaces["refgenie"] = AttMap(paths_dict)
        return rgc.populate(namespaces)
    else:
        msg = """
        var_templates:
          refgenie_config: "$REFGENIE"
        """
        _LOGGER.error(
            f"refgenie_config not specified in pipeline interface. Do like so: {msg}"
        )
        raise NotImplementedError
