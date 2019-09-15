"""The main entrypoints live here."""
import logging
import os

import numpy
import geopandas

from . import log, model

LOGGER = logging.getLogger(__name__)

__folder__ = os.path.dirname(os.path.realpath(__file__))
__data_folder__ = os.path.abspath(os.path.join(os.path.dirname( __folder__ ), '..', 'chlor_a'))


def main():
    log.init()

    # Assignment Simulation:
    LOGGER.info("Loading data from disk...")

    for layer in __data_folder__:
        print(layer)

    return

    with open("/tmp/data.csv", "r") as fileobj:
        raw_data = geopandas.read_csv(fileobj)
        X = raw_data[raw_data.columns.difference(["is_bad"])]
        y = numpy.array(raw_data[["is_bad"]])

    xgboost_model = model.XGBoostModel()

    LOGGER.info(xgboost_model.evaluate(X, y))
    # Assuming X is a pandas dataframe...
    # TODO: Make note that dask is not talking to local threads very well. Not an actual error.
    LOGGER.info(
        xgboost_model.predict_proba(geopandas.DataFrame(xgboost_model.model.X.compute()))
    )
