import platform
import sys
from collections import OrderedDict

import yaml as pyyaml


_items = "viewitems" if sys.version_info < (3,) else "items"
_std_dict_is_order_preserving = sys.version_info >= (3, 7) or (
    sys.version_info >= (3, 6) and platform.python_implementation() == "CPython"
)


def map_representer(dumper, data):
    return dumper.represent_dict(getattr(data, _items)())


def map_constructor(loader, node):
    loader.flatten_mapping(node)
    pairs = loader.construct_pairs(node)
    try:
        return OrderedDict(pairs)
    except TypeError:
        loader.construct_mapping(node)  # trigger any contextual error
        raise


_loaders = [getattr(pyyaml.loader, x) for x in pyyaml.loader.__all__]
_dumpers = [getattr(pyyaml.dumper, x) for x in pyyaml.dumper.__all__]
try:
    _cyaml = pyyaml.cyaml.__all__
except AttributeError:
    pass
else:
    _loaders += [getattr(pyyaml.cyaml, x) for x in _cyaml if x.endswith("Loader")]
    _dumpers += [getattr(pyyaml.cyaml, x) for x in _cyaml if x.endswith("Dumper")]

Dumper = None
for Dumper in _dumpers:
    pyyaml.add_representer(dict, map_representer, Dumper=Dumper)
    pyyaml.add_representer(OrderedDict, map_representer, Dumper=Dumper)

Loader = None
if not _std_dict_is_order_preserving:
    for Loader in _loaders:
        pyyaml.add_constructor("tag:yaml.org,2002:map", map_constructor, Loader=Loader)


# Merge PyYAML namespace into ours.
# This allows users a drop-in replacement:
#   import oyaml as yaml
del map_constructor, map_representer, Loader, Dumper
from yaml import *
