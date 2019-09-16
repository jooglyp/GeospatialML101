"""Module for creating a hexgrid and aggregating point data within that hexgrid."""

import logging
import os

import geopandas as gpd
import shapely
from shapely.geometry import Polygon, Point, box
import math
import shapefile
import pandas as pd
import numpy as np

LOGGER = logging.getLogger(__name__)


class Hexify:
    def __init__(self, gdf: gpd.GeoDataFrame, output_dir):
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

        self.output_dir = output_dir
        self.points_shp = gdf
        self.points_shp = self.points_shp[self.points_shp.geometry.notnull()]
        if self.points_shp.crs != self.us_albers_equal_area:
            print("converting crs of points layer")
            self.points_shp.to_crs(self.us_albers_equal_area, inplace=True)

    def calculate_polygons(self, startx, starty, endx, endy, radius):
        """
        Calculate a grid of hexagon coordinates of the given radius
        given lower-left and upper-right coordinates
        Returns a list of lists containing 6 tuples of x, y point coordinates
        These can be used to construct valid regular hexagonal polygons

        You will probably want to use projected coordinates for this
        """
        # calculate side length given radius
        sl = (2 * radius) * math.tan(math.pi / 6)
        # calculate radius for a given side-length
        # (a * (math.cos(math.pi / 6) / math.sin(math.pi / 6)) / 2)
        # see http://www.calculatorsoup.com/calculators/geometry-plane/polygon.php

        # calculate coordinates of the hexagon points
        # sin(30)
        p = sl * 0.5
        b = sl * math.cos(math.radians(30))
        w = b * 2
        h = 2 * sl

        # offset start and end coordinates by hex widths and heights to guarantee coverage
        startx = startx - w
        starty = starty - h
        endx = endx + w
        endy = endy + h

        origx = startx
        origy = starty

        # offsets for moving along and up rows
        xoffset = b
        yoffset = 3 * p

        polygons = []
        row = 1
        counter = 0

        while starty < endy:
            if row % 2 == 0:
                startx = origx + xoffset
            else:
                startx = origx
            while startx < endx:
                p1x = startx
                p1y = starty + p
                p2x = startx
                p2y = starty + (3 * p)
                p3x = startx + b
                p3y = starty + h
                p4x = startx + w
                p4y = starty + (3 * p)
                p5x = startx + w
                p5y = starty + p
                p6x = startx + b
                p6y = starty
                poly = [
                    [p1x, p1y],
                    [p2x, p2y],
                    [p3x, p3y],
                    [p4x, p4y],
                    [p5x, p5y],
                    [p6x, p6y],
                    [p1x, p1y]]
                polygons.append(poly)
                counter += 1
                startx += w
            starty += yoffset
            row += 1
        return polygons

    def geoms_to_hex_shp(self, in_geoms, name: str, hex_area: float, projection):
        prj_name = '{}.prj'.format('{}/{}_{}_Hex'.format(self.output_dir, name, hex_area))
        with open(prj_name, 'w') as prj:
            prj.write(projection)
        shp_writer = shapefile.Writer(os.path.join(self.output_dir, "HexGrid"))

        out_fields = [
            ['id', 'N']
        ]
        for name in out_fields:
            shp_writer.field(*name)

        for in_id, geom in enumerate(in_geoms, start=1):
            shp_writer.record(*[str(in_id)])
            shp_writer.poly([geom])

    def create_hex_grid(self, hex_area, name: str, saveout=True):
        """

        :param hex_area: should be around 5000-8000 meters squared
        :param saveout: if True, will save-out the hex-grid as a shp to output_dir
        :return: a hex-grid stored in a dictionary
        """

        # hex_area = 86602540.378  # 86602540.378 is the area of a hexagon with radius = 5000 meters
        bbox = self.points_shp['geometry'].total_bounds

        geom_box = box(*bbox)
        extent_area = geom_box.area

        if extent_area < hex_area:  # 86602540.378 is the area of a hexagon with radius = 5000 meters
            # begin rescaling the extent_area
            # (1) The following condition, box(a,b,c,d).area/box(w,s,e,n).area = xfact*yfact.
            # (2) want to know --> (86602540.378*1.1)/box(w,s,e,n).area --> sqrt(*)
            req_growth = (hex_area * 1.1) / extent_area
            scale_factor = req_growth ** 0.5
            nshp_bounds = shapely.affinity.scale(geom_box, xfact=scale_factor, yfact=scale_factor).bounds
            LOGGER.info(nshp_bounds)
            # obtain new bounds
            bbox = nshp_bounds

        # # # Create the hexgrid # # #
        w, s, e, n = bbox
        vertices = self.calculate_polygons(w, s, e, n, hex_area)  # 5000 meter radius (UTM)

        # Create the hex geodataframe
        v_dict = {ndx: Polygon(v) for ndx, v in enumerate(vertices)}
        hex_gdf = gpd.GeoDataFrame.from_dict({'geometry': v_dict})
        hex_gdf.crs = self.points_shp.crs
        hex_gdf = hex_gdf.reset_index().rename(columns={'index': 'polyid'})

        # store in a zipped list (file_name will keep track of which hex-region is being stored)
        LOGGER.info("Generating hex grid for: %s" % name)
        LOGGER.info("================================")
        # write hex to shapefile for debugging
        if saveout:
            self.geoms_to_hex_shp([list(s.exterior.coords) for s in hex_gdf['geometry'].tolist()],
                                  name, hex_area, self.us_albers_equal_area_prj)
        return hex_gdf

    def point_to_avg_hex_scores(self, points, hex_grid, saveout=True):
        """

        :param points: points to perform aggregations on
        :return: a new hex-list with a z-score
        """

        if hex_grid.crs != self.us_albers_equal_area:
            hex_grid.to_crs(self.us_albers_equal_area, inplace=True)

        if points.crs != self.us_albers_equal_area:
            points.to_crs(self.us_albers_equal_area, inplace=True)

        points_in_hexes = gpd.sjoin(points, hex_grid, how='left')
        # Perform Math:
        attributes = points_in_hexes.columns
        attributes = [i for i in attributes if i.startswith('chlor')]
        hex_avgs = {}
        for attribute in attributes:
            point_hex_avg = points_in_hexes.groupby(['polyid'])[attribute].mean()
            point_hex_avg_rc = pd.DataFrame({attribute.split("_")[1] + "avg": point_hex_avg}).reset_index()
            point_hex_avg_rc.set_index('polyid', inplace=True)
            LOGGER.info(point_hex_avg_rc)
            hex_avgs[attribute.split("_")[1] + "avg"] = point_hex_avg_rc
        for item in hex_avgs.items():
            hex_grid = pd.merge(hex_grid, item[1], left_on='polyid', right_index=True,
                                how='left', sort=False)
        if saveout:
            out_dir = os.path.join(self.output_dir, "Hex_Values.shp")
            # note importing in qgis seems to screw up the column names.
            hex_grid.to_file(out_dir)
        return hex_grid

    def hex_centroid_values(self, hex_grid: gpd.GeoDataFrame, saveout=True) -> gpd.geodataframe:
        cg_long = hex_grid.geometry.centroid.x
        cg_lat = hex_grid.geometry.centroid.y
        cg_df = pd.DataFrame(
            {
                'Latitude': list(cg_lat),
                'Longitude': list(cg_long)
            }
        )
        cg_df['Coordinates'] = list(zip(cg_df.Longitude, cg_df.Latitude))
        cg_df['Coordinates'] = cg_df['Coordinates'].apply(Point)
        cg_gdf = gpd.GeoDataFrame(cg_df, geometry='Coordinates')
        # cg_gdf = gpd.GeoDataFrame(cg_df, geometry=gpd.points_from_xy(cg_df.Longitude, cg_df.Latitude))
        if cg_gdf.crs != self.us_albers_equal_area:
            print("converting crs of points layer")
            hex_grid.to_crs(self.us_albers_equal_area, inplace=True)
        if saveout:
            out_dir = os.path.join(self.output_dir, "Hex_Values_Points.shp")
            hex_grid.to_file(out_dir)
        return cg_gdf