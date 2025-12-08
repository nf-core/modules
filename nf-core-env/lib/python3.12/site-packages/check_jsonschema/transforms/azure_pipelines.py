"""
A data transform which unpacks "compile-time expressions" from Azure Pipelines files.

For the original source which inspired this transform, see the YAML parser used in the
Azure Pipelines Language Server:
    https://github.com/microsoft/azure-pipelines-language-server/blob/71b20f92874c02dfe82ad2cc2dcc7fa64996be91/language-service/src/parser/yamlParser.ts#L182

That source is licensed under the MIT License.
The original license can be found in
    src/check_jsonschema/builtin_schemas/vendor/licenses/LICENSE.azure-pipelines


The transform does not deeply interpret the expressions. It just "unnests" them.

It will turn this input

    jobs:
      - ${{ each val in parameter.vals }}:
        - job: foo
          steps:
            - bash: echo ${{ val }}

into

    jobs:
      - job: foo
        steps:
          - bash: echo ${{ val }}
"""

from __future__ import annotations

import typing as t

from .base import Transform


class AzurePipelinesDataError(ValueError):
    def __init__(self, message: str) -> None:
        super().__init__(f"azure-pipelines transform: {message}")


def is_expression(s: str) -> bool:
    return s.startswith("${{") and s.endswith("}}")


def traverse_data(data: t.Any) -> t.Any:
    if isinstance(data, dict):
        return traverse_dict(data)
    if isinstance(data, list):
        return traverse_list(data)
    return data


def traverse_list(data: list) -> list:
    ret = []
    for item in data:
        # is the current item a single-value dict with an expression as its key?
        item_is_expr = (
            isinstance(item, dict)
            and len(item) == 1
            and is_expression(tuple(item)[0])  # tuple() gets keys
        )

        if item_is_expr:
            # unpack the expression item and recurse over the value
            item_key, item_value = list(item.items())[0]
            item_value = traverse_data(item_value)

            if isinstance(item_value, list):
                ret.extend(item_value)
            else:
                ret.append(item_value)
        # not expression? process the item and append
        else:
            ret.append(traverse_data(item))
    return ret


def traverse_dict(data: dict) -> dict:
    newdata = {}
    for key, value in data.items():
        newvalue = traverse_data(value)
        if is_expression(key):
            # WARNING -- correctness unclear
            #
            # "lift" any dict by moving its attributes up into the object being evaluated
            #
            # e.g.
            #    parent:
            #      ${{ each x in xs }}:
            #        - k: v-${{ x }}
            #
            # becomes
            #
            #    parent:
            #      - k: v-${{ x }}
            if isinstance(newvalue, dict):
                for add_k, add_v in newvalue.items():
                    newdata[add_k] = add_v
            # In all other cases, drop the content from the data. This is based on the
            # azure-pipelines-language server behavior:
            # https://github.com/microsoft/azure-pipelines-language-server/blob/71b20f92874c02dfe82ad2cc2dcc7fa64996be91/language-service/src/parser/yamlParser.ts#L185
            #
            # earlier versions would raise an error here, but this caused issues with
            # data in which expressions were mapped to simple strings
            #
            # e.g.
            #
            #    parent:
            #      ${{ x }}: ${{ y }}
            #
            # which occurs naturally *after* a lifting operation, as in
            #
            #    parent:
            #      ${{ each x, y in attrs }}:
            #        ${{ x }}: ${{ y }}
            else:
                continue
        else:
            newdata[key] = newvalue
    return newdata


def azure_main(data: dict | list) -> dict | list:
    if isinstance(data, list):
        raise AzurePipelinesDataError(
            "this transform requires that the data be an object, got list"
        )
    return traverse_dict(data)


AZURE_TRANSFORM = Transform(on_data=azure_main)
