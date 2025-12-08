from abc import ABC, abstractmethod
from typing import Iterable, List, Union

from peppy.sample import Sample


class BaseOutputFormatter(ABC):
    @staticmethod
    @abstractmethod
    def format(samples: List[Sample]):
        """
        Convert the samples to correct format.
        """
        pass


class MultilineOutputFormatter(BaseOutputFormatter):
    @staticmethod
    def format(samples: List[Sample]) -> str:
        output_rows = []
        sample_attributes = [
            attribute
            for attribute in samples[0].keys()
            if not attribute.startswith("_") and not attribute == "subsample_name"
        ]
        header = MultilineOutputFormatter._get_header(sample_attributes)

        for sample in samples:
            attribute_with_multiple_properties = MultilineOutputFormatter._get_the_name_of_the_first_attribute_with_multiple_properties(
                sample, sample_attributes
            )
            if attribute_with_multiple_properties:
                sample_rows = MultilineOutputFormatter._split_sample_to_multiple_rows(
                    sample, sample_attributes, attribute_with_multiple_properties
                )
                output_rows.extend(sample_rows)
            else:
                one_sample_row = MultilineOutputFormatter._convert_sample_to_row(
                    sample, sample_attributes
                )
                output_rows.append(one_sample_row)

        return "\n".join(header + output_rows) + "\n"

    @staticmethod
    def _get_header(header_column_names: List[str]):
        return [",".join(header_column_names)]

    @staticmethod
    def _get_the_name_of_the_first_attribute_with_multiple_properties(
        sample: Sample, sample_attributes: List[str]
    ) -> Union[str, None]:
        for attribute in sample_attributes:
            if MultilineOutputFormatter._sample_attribute_is_list(sample, attribute):
                return attribute

    @staticmethod
    def _split_sample_to_multiple_rows(
        sample: Sample, sample_attributes: List, attribute_with_multiple_properties: str
    ) -> Iterable[str]:
        """
        If one sample object contains array properties instead of single value, then it will be converted
        to multiple rows.

        Args:
            sample: Sample from project.
            sample_attributes: List of all sample properties names (name of columns from sample_table).

        Returns:
            List of rows created from given sample object.
        """
        number_of_samples_after_split = len(
            getattr(sample, attribute_with_multiple_properties)
        )
        sample_rows_after_split = []

        for sample_index in range(number_of_samples_after_split):
            sample_row = MultilineOutputFormatter._convert_sample_to_row(
                sample, sample_attributes, sample_index
            )
            sample_rows_after_split.append(sample_row)

        return sample_rows_after_split

    @staticmethod
    def _convert_sample_to_row(
        sample: Sample, sample_attributes: List, sample_index: int = 0
    ) -> str:
        """
        Converts single sample object to CSV row.

        Some samples have a list of values instead of single value for given attribute (column), and
        sample_index indicates index of the value that will be used to create a row. For samples that don't
        have any attributes with given names this will always be zero.

        Args:
            sample: Single sample object.
            sample_attributes: Array of all attributes names (column names) for given sample.
            sample_index: Number indicating which value will be used to create row. Some samples

        Returns:
            Representation of sample as a CSV row.
        """
        sample_row = []

        for attribute in sample_attributes:
            if (
                MultilineOutputFormatter._sample_attribute_is_list(sample, attribute)
                and sample[attribute]
            ):
                value = sample[attribute][sample_index]
            else:
                value = sample.get(attribute)

            sample_row.append(value or "")

        return ",".join(sample_row)

    @staticmethod
    def _sample_attribute_is_list(sample: Sample, attribute: str) -> bool:
        return isinstance(getattr(sample, attribute, ""), list)


class SampleSubsampleOutputFormatter(BaseOutputFormatter):
    def format(self, samples: List[Sample]):
        pass
