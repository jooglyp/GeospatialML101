"""The main entrypoints live here."""
import logging
import os

import pandas
import geopandas
from shapely.geometry import Point

from . import log, model, utilities, hexes

LOGGER = logging.getLogger(__name__)

__folder__ = os.path.dirname(os.path.realpath(__file__))
__data_folder__ = os.path.abspath(os.path.join(os.path.dirname( __folder__ ), 'chlor_a/'))


def data_generation():
    out_dir = os.path.join(__folder__, "..", "outputs")
    data_generator = model.GenerateCovariates(os.path.join(out_dir, "Hex_Values_Points.shp"))
    high_winter, seasonal_yoys, interannual_variance, seasonal_variance = data_generator.generate_covariates()
    data_generator.generate_hotspots(high_winter)
    data_generator.generate_spatial_association(seasonal_variance, ['summer_var', 'winter_var'])
    data_generator.generate_spatial_association(interannual_variance, ['2003_var', '2004_var',
                      '2005_var', '2006_var', '2007_var', '2008_var', '2009_var', '2010_var', '2011_var',
                      '2012_var', '2013_var', '2014_var', '2015_var', '2016_var', '2017_var', '2018_var',
                      '2019_var'])
    data_generator.generate_spatial_association(seasonal_yoys, ['2019_winter_yoy', '2019_summer_yoy',
       '2018_winter_yoy', '2018_summer_yoy', '2017_winter_yoy',
       '2017_summer_yoy', '2016_winter_yoy', '2016_summer_yoy',
       '2015_winter_yoy', '2015_summer_yoy', '2014_winter_yoy',
       '2014_summer_yoy', '2013_winter_yoy', '2013_summer_yoy',
       '2012_winter_yoy', '2012_summer_yoy', '2011_winter_yoy',
       '2011_summer_yoy', '2010_winter_yoy', '2010_summer_yoy',
       '2009_winter_yoy', '2009_summer_yoy', '2008_winter_yoy',
       '2008_summer_yoy', '2007_winter_yoy', '2007_summer_yoy',
       '2006_winter_yoy', '2006_summer_yoy', '2005_winter_yoy',
       '2005_summer_yoy', '2004_winter_yoy', '2004_summer_yoy',
       '2003_winter_yoy', '2003_summer_yoy'])


def data_processing():
    log.init()
    try:
        os.mkdir(os.path.join(os.path.dirname( __folder__ ), "outputs"))
    except FileExistsError:
        LOGGER.info("Output directory already exists, continuing...")

    # Assignment Data Processing:
    LOGGER.info("Loading data from disk...")
    df_list = []
    geo_files = scan_data()
    i = 0
    for file in geo_files.items():
        if i == 3:
            break
        with open(file[1], "r") as fileobj:
            flattened_dict = flatten_data(fileobj)
            df = pandas.DataFrame(flattened_dict)
            name = "chlor_" + str(file[0].split(".")[0].split("_")[2])
            print(name)
            df.rename(columns={"value": name}, inplace=True)
            df_list.append(df)
            i += 1
    geoconstructor = utilities.GeoDataConstructor(df_list)
    concatenated_df = geoconstructor._concatenate_dataframes()
    gdf = geopandas.GeoDataFrame(
        concatenated_df.drop(['latitude', 'longitude'], axis=1),
        crs={'init': 'epsg:4326'},
        geometry=[Point(xy) for xy in zip(concatenated_df.latitude, concatenated_df.longitude)])

    del concatenated_df
    del df_list
    hex_attribute_points = simplify_data(gdf)


def simplify_data(gdf: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
    out_dir = os.path.join(__folder__, "..", "outputs")
    # default_hex_size = 40000
    default_hex_size = 100000  # testing only
    grid_obj = hexes.Hexify(gdf, out_dir)
    hex_grid = grid_obj.create_hex_grid(default_hex_size, "atlantic_grid", saveout=True)
    LOGGER.info(hex_grid)
    hex_attributes = grid_obj.point_to_avg_hex_scores(gdf, hex_grid)
    hex_attribute_points = grid_obj.hex_centroid_values(hex_attributes)
    return hex_attribute_points


def write_data(gdf: geopandas.GeoDataFrame):
    """Writes data to directory as shapefile for manual inspection."""
    out_dir = os.path.join(__folder__, "..", "outputs", "chlor_timeseries.gpkg")
    gdf.to_file(out_dir, driver='GPKG')


def scan_data() -> dict:
    """Generates file directory structure."""
    geo_files = {}
    for file in os.listdir(__data_folder__):
        if file.endswith(".csv"):
            geo_files[file] = os.path.join(__data_folder__, file)
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
