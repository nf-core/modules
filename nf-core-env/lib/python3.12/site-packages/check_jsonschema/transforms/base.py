from __future__ import annotations

import typing as t

import ruamel.yaml


class Transform:
    def __init__(
        self,
        *,
        on_data: t.Callable[[list | dict], list | dict] | None = None,
    ) -> None:
        self.on_data = on_data

    def modify_yaml_implementation(self, implementation: ruamel.yaml.YAML) -> None:
        pass

    def __call__(self, data: list | dict) -> list | dict:
        if self.on_data is not None:
            return self.on_data(data)
        return data
