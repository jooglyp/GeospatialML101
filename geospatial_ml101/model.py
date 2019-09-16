"""Generates seasonal statistics from data."""

import logging
import os

import geopandas
import pandas
import numpy
from .hotspots import HotSpots

LOGGER = logging.getLogger(__name__)

__folder__ = os.path.dirname(os.path.realpath(__file__))


class GenerateCovariates:
    def __init__(self, hex_path: str):
        # Designate commonly used projection systems as class attributes:
        self.us_albers_equal_area_prj = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic",\
        GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",\
        SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],\
        PROJECTION["Albers"],PARAMETER["False_Easting",0.0],\
        PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-96.0],\
        PARAMETER["Standard_Parallel_1",29.5],PARAMETER["Standard_Parallel_2",45.5],\
        PARAMETER["Latitude_Of_Origin",37.5],UNIT["Meter",1.0]]'
        self.wgs84_prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],\
        PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'
        self.us_albers_equal_area = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 \
        +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
        self.wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

        self.path = hex_path
        self.hex_point_layer = geopandas.GeoDataFrame.from_file(hex_path)

    def get_season(self, day) -> str:
        seasons = {(0, 104): 'winter', (104, 104 + 184): 'summer', (104 + 184, 104 + 184 + 77): 'winter'}
        for season_range, season in seasons.items():
            if day in range(*season_range):
                return season

    def unique(self, list1):
        x = numpy.array(list1)
        return numpy.unique(x)

    def seasonal_average(self, hex_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        hex_data_copy = hex_data.copy()
        new_hex_data = hex_data_copy[['polyid', 'geometry']]
        unique_seasons = hex_data_copy.columns
        unique_seasons = set(unique_seasons) - set(list(["geometry", "polyid"]))
        seasons = self.unique([i.split("_")[1] for i in unique_seasons])

        for season in seasons:
            relevant_fields = [x for x in unique_seasons if x.endswith('{}'.format(season))]
            season_df = hex_data_copy[list(relevant_fields)]
            variance = season_df.mean(axis=1)
            new_hex_data['{}_mean'.format(season)] = variance
        print(new_hex_data)
        return new_hex_data

    def seasonal_variance(self, hex_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        hex_data_copy = hex_data.copy()
        new_hex_data = hex_data_copy[['polyid', 'geometry']]
        unique_seasons = hex_data_copy.columns
        unique_seasons = set(unique_seasons) - set(list(["geometry", "polyid"]))
        seasons = self.unique([i.split("_")[1] for i in unique_seasons])

        for season in seasons:
            relevant_fields = [x for x in unique_seasons if x.endswith('{}_yoy'.format(season))]
            season_df = hex_data_copy[list(relevant_fields)]
            variance = season_df.var(axis=1)
            new_hex_data['{}_var'.format(season)] = variance
        print(new_hex_data)
        return new_hex_data

    def interannual_variance(self, hex_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        hex_data_copy = hex_data.copy()
        new_hex_data = hex_data_copy[['polyid', 'geometry']]
        unique_seasons = hex_data_copy.columns
        unique_seasons = set(unique_seasons) - set(list(["geometry", "polyid"]))
        years = self.unique([i.split("_")[0] for i in unique_seasons])

        for year in years:
            relevant_fields = [x for x in unique_seasons if x.startswith(year)]
            year_df = hex_data_copy[list(relevant_fields)]
            variance = year_df.var(axis=1)
            new_hex_data['{}_var'.format(year)] = variance
        print(new_hex_data)
        return new_hex_data

    def deseasonalization(self, hex_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        hex_data = hex_data.copy()
        deseasonalized = hex_data.copy()
        deseasonalized = deseasonalized[['polyid', 'geometry']]
        year_seasons = set(list(hex_data.columns.sort_values())) - set(list(["geometry", "polyid"]))
        year_seasons = list(year_seasons)
        year_seasons.sort(reverse=True)
        for season in year_seasons:
            try:
                season_ = season.split("_")[1]
                current_year = int(season.split("_")[0])
                last_year = int(season.split("_")[0]) - 1
                deseasonalized['{}_{}_yoy'.format(str(current_year), season_)] = \
                    (hex_data['{}_{}'.format(str(current_year), season_)] -
                     hex_data['{}_{}'.format(str(last_year), season_)]) / \
                    hex_data['{}_{}'.format(str(last_year), season_)]
            except:
                continue
        return deseasonalized

    def average_year_seasons(self, hex_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        hex_data_copy = hex_data.copy()
        new_hex_data = hex_data_copy[['polyid', 'geometry']]
        unique_seasons = self.unique(hex_data_copy.columns)
        unique_seasons = set(unique_seasons) - set(list(["geometry", "polyid"]))
        for year_season in unique_seasons:
            df = hex_data_copy.copy(deep=True)[year_season]
            if type(df) == pandas.Series:
                new_hex_data[year_season] = df
                continue
            df['polyid'] = hex_data_copy['polyid']
            df = df.set_index(['polyid'])
            df = df.groupby(by=df.columns, axis=1).mean()
            df = df.reset_index()
            new_hex_data[year_season] = df[year_season]
        print(new_hex_data)
        return new_hex_data

    def coerce_to_numeric(self, hex_data: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        for column in hex_data.columns:
            if column[-3:] == 'avg':
                series = self.hex_point_layer[column]
                hex_data[column] = pandas.to_numeric(series)
        return hex_data

    def analyze_bias(self, df: geopandas.GeoDataFrame):
        return df.isna().sum()

    def generate_covariates(self):
        print(self.hex_point_layer)
        self.hex_point_layer = self.coerce_to_numeric(self.hex_point_layer)

        new_name_dict = {}
        for column in self.hex_point_layer.columns:
            if column[-3:] == 'avg':
                year = column[0:4]
                day = column[4:7]
                season = self.get_season(int(day))
                rename = str(year) + "_" + str(season)
                new_name_dict[column] = rename
        self.hex_point_layer.rename(columns=new_name_dict, inplace=True)

        seasonals = self.average_year_seasons(self.hex_point_layer.copy())  # Answers Question 1
        seasonal_yoys = self.deseasonalization(seasonals.copy())
        seasonal_bias = self.analyze_bias(seasonals)  # Answers Question 2
        seasonal_yoys_bias = self.analyze_bias(seasonal_yoys)  # Answers Question 2
        interannual_variance = self.interannual_variance(seasonal_yoys)  # Answers Question 3
        seasonal_variance = self.seasonal_variance(seasonal_yoys)  # Answers Question 3
        seasonal_mean = self.seasonal_average(seasonals)

        out_dir = os.path.join(__folder__, "..", "outputs")
        seasonals.to_file(os.path.join(out_dir, "metrics_seasonals.shp"))
        seasonal_yoys.to_file(os.path.join(out_dir, "metrics_seasonal_yoy.shp"))
        interannual_variance.to_file(os.path.join(out_dir, "metrics_interannual_variance.shp"))
        seasonal_variance.to_file(os.path.join(out_dir, "metrics_seasonal_variance.shp"))

        high_winter = seasonal_mean[['polyid', 'winter_mean', 'geometry']]
        return high_winter

    def generate_hotspots(self, gdf: geopandas.GeoDataFrame):
        out_dir = os.path.join(__folder__, "..", "outputs")
        hex_dir = os.path.join(__folder__, "..", "outputs", 'HexGrid.shp')
        hotspots = HotSpots(gdf, hex_dir, self.path)
        scored_points = hotspots.create_scores()
        scored_points = scored_points[['id', 'geometry', 'winter_mean', 'Z_Scores', 'Heat_Scores']]
        scored_points.to_file(os.path.join(out_dir, "winter_mean_hotspots.shp"))  # Answers Decision Making Problem
