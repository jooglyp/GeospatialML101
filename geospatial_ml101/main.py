"""The main entrypoints live here."""
import logging
import os

import pandas
import geopandas
from shapely.geometry import Point

from . import log, model, utilities

LOGGER = logging.getLogger(__name__)

__folder__ = os.path.dirname(os.path.realpath(__file__))
__data_folder__ = os.path.abspath(os.path.join(os.path.dirname( __folder__ ), 'chlor_a/'))


def main():
    log.init()
    try:
        os.mkdir(os.path.join(os.path.dirname( __folder__ ), "outputs"))
    except FileExistsError:
        LOGGER.info("Output directory already exists, continuing...")

    # Assignment Data Processing:
    LOGGER.info("Loading data from disk...")
    gdf_list = []
    geo_files = scan_data()
    for file in geo_files.items():
        print(file)
        with open(file[1], "r") as fileobj:
            flattened_dict = flatten_data(fileobj)
            df = pandas.DataFrame(flattened_dict)
            name = "chlor_" + str(file[0].split(".")[0].split("_")[2])
            gdf = geopandas.GeoDataFrame(
                df.drop(['latitude', 'longitude'], axis=1),
                crs={'init': 'epsg:4326'},
                geometry=[Point(xy) for xy in zip(df.latitude, df.longitude)])
            gdf.rename(columns={"value": name}, inplace=True)
            """
            out_dir = os.path.join(__folder__, "..", "outputs",
                                   "chlor_" + str(file[0].split(".")[0].split("_")[2]) + ".shp")
            gdf.to_file(out_dir)
            break
            """
            gdf_list.append(gdf)

    geodata = utilities.GeoDataConstructor(gdf_list)
    concatenated_gdf = geodata._concatenate_dataframes()
    print(concatenated_gdf)
    return


def scan_data() -> dict:
    geo_files = {}
    for file in os.listdir(__data_folder__):
        if file.endswith(".csv"):
            geo_files[file] = os.path.join(__data_folder__, file)
    print(geo_files)
    return geo_files


def flatten_data(csvfile) -> dict:
    """

    :param fileobj: csv file handled in context.
    :return: transformed csv file to dictionary that can be read-into geopandas.
    """
    df = pandas.read_csv(csvfile)

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
