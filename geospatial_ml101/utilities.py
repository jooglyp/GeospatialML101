"""Utilities for constructing geopadas dataset."""

import functools
import logging

import pandas

LOGGER = logging.getLogger(__name__)


class GeoDataConstructor:
    def __init__(self, geodata: list):
        self.geodata = geodata
        self.model_covariates = None

    def _concatenate_dataframes(self) -> pandas.DataFrame:
        """
        All dataframes must be of the same size.
        Args:
            dataframes: a list of dataframes to perform dataframe concatenation on.

        Returns:

        """
        LOGGER.info("Concatenating Dataframes")
        datas = functools.reduce(
            lambda left, right: pandas.merge(
                left, right, left_index=True, right_index=True, how="outer"
            ),
            self.geodata,
        )
        LOGGER.info("Final Dataset:")
        self.model_covariates = datas.columns
        return datas
