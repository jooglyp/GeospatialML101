"""The main entrypoints live here."""
import logging
import os

import csv
import pandas
import geopandas
from shapely.geometry import Point

from . import log, model

LOGGER = logging.getLogger(__name__)

__folder__ = os.path.dirname(os.path.realpath(__file__))
__data_folder__ = os.path.abspath(os.path.join(os.path.dirname( __folder__ ), 'chlor_a/'))


def main():
    log.init()

    # Assignment Simulation:
    LOGGER.info("Loading data from disk...")
    geo_files = {}
    for file in os.listdir(__data_folder__):
        if file.endswith(".csv"):
            geo_files[file] = os.path.join(__data_folder__, file)
    print(geo_files)

    for file in geo_files.items():
        print(file)
        with open(file[1], "r") as fileobj:
            flattened_dict = flatten_data(fileobj)
            df = pandas.DataFrame(flattened_dict)
            print(df)
            gdf = geopandas.GeoDataFrame(
                df.drop(['latitude', 'longitude'], axis=1),
                crs={'init': 'epsg:4326'},
                geometry=[Point(xy) for xy in zip(df.latitude, df.longitude)])
            print(gdf)
            out_dir = os.path.join(__folder__, "Test.shp")
            gdf.to_file(out_dir)
            break


def flatten_data(csvfile) -> dict:
    """

    :param fileobj: csv file handled in context.
    :return: transformed csv file to dictionary that can be read-into geopandas.
    """
    print("================")
    df = pandas.read_csv(csvfile)
    # latitudes = list(df.columns[1:])
    # print(latitudes)
    # longitudes = list(df[df.columns[0]])
    # print(longitudes)
    # print("================")

    transformed_data = {'latitude': [], 'longitude': [], 'value': []}
    for x in df.iterrows():
        for i, pixel in enumerate(x[1].iteritems()):
            if i == 0:
                longitude = pixel[1]
                print(longitude)
                continue
            else:
                transformed_data['latitude'].append(float(pixel[0]))
                transformed_data['longitude'].append(float(longitude))
                transformed_data['value'].append(float(pixel[1]))
    return transformed_data
