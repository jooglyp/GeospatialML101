"""Generates seasonal statistics from data."""

import geopandas
import pandas


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
        self.hex_attributes = geopandas.GeoDataFrame.from_file(hex_path)

    def __stuff__(self):
        return None