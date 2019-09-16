"""Module to create hotspot scores for a particular covariate.
"""

import logging

import pandas
import geopandas
import pysal
from pysal.explore.esda.getisord import G_Local
from shapely.geometry import Point
import numpy as np

# ignore numpy floating point errors
np.seterr(all='ignore')

LOGGER = logging.getLogger(__name__)


class HotSpots:
    def __init__(self, sample_points, hex_path, output_dir):
        # Designate commonly used projection systems as global attributes:
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
        self.output_dict = None
        self.sample_points = sample_points
        self.hex_template = geopandas.GeoDataFrame.from_file(hex_path)

    def min_threshold_dist_from_points(self, in_geodataframe, radius=None, p=2):
        """
        Source Code:
        http://pysal.readthedocs.io/en/latest/_modules/pysal/weights/Distance.html

        MODIFIED VERSION OF 'pysal.min_threshold_dist_from_shapefile'
        NOW TAKES GEODATAFRAME INSTEAD OF SHAPEFILE
        ----------
        Kernel weights with adaptive bandwidths.
        Parameters
        ----------
        shapefile  : string
                    shapefile name with shp suffix.
        radius     : float
                    If supplied arc_distances will be calculated
                    based on the given radius. p will be ignored.
        p          : float
                    Minkowski p-norm distance metric parameter:
                    1<=p<=infinity
                    2: Euclidean distance
                    1: Manhattan distance
        Returns
        -------
        d          : float
                    Maximum nearest neighbor distance between the n
                    observations.
        """
        points = np.vstack([np.array(shape.centroid) for shape in in_geodataframe['geometry']])
        if radius is not None:
            kdt = pysal.lib.cg.kdtree.Arc_KDTree(points, radius=radius)
            nn = kdt.query(kdt.data, k=2)
            nnd = nn[0].max(axis=0)[1]
            return nnd
        return pysal.lib.weights.min_threshold_distance(points, p)

    def max_distance_points(self, points):
        sub_gdf = points['geometry'].copy()
        out_distances = []
        for ndx_1, geom_1 in sub_gdf.iteritems():
            for ndx_2, geom_2 in sub_gdf.iteritems():
                dist = geom_1.distance(geom_2)
                out_distances.append(dist)
        max_dist = max(out_distances)
        return max_dist

    def transform_dataframe(self, hex_grid: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
        hex_grid = geopandas.GeoDataFrame(hex_grid, geometry='geometry', crs=self.us_albers_equal_area)
        cg_long = hex_grid.geometry.centroid.x
        cg_lat = hex_grid.geometry.centroid.y
        cg_df = pandas.DataFrame(
            {
                'Latitude': list(cg_lat),
                'Longitude': list(cg_long)
            }
        )
        cg_df['Coordinates'] = list(zip(cg_df.Longitude, cg_df.Latitude))
        cg_df['geometry'] = cg_df['Coordinates'].apply(Point)
        cg_gdf = geopandas.GeoDataFrame(cg_df, geometry='geometry', crs=self.us_albers_equal_area)
        cg_gdf['winter_mean'] = hex_grid['winter_mean']
        return cg_gdf

    def create_scores(self):
        self.output_dict = {}
        point = self.sample_points.copy()
        point = self.transform_dataframe(point)

        region = 'atlantic'
        if len(point) < 200:  # low number of points may have sparse r-trees
            grid_thresh = self.max_distance_points(point)
        else:
            grid_thresh = self.min_threshold_dist_from_points(point)

        # the spatial weights matrix will yield greater agglomerations as the parameter p increases
        # (see min_threshold_dist_from_points(point)
        print("pass 2: resetting indices")
        point = point.reset_index()
        print("pass 3: calculating inverse distance weights matrix")
        w = pysal.lib.weights.DistanceBand.from_dataframe(point, threshold=grid_thresh, alpha=-1.5, binary=False)
        print("pass 4: applying row standardization to inverse distance weights matrix")
        w.transform = "R"  # row standardized weights (see pysal handbook pg. 171)
        # GETIS-ORD Z STATISTIC:
        print("pass 5: extracting winter_mean as y-attribute")
        heat_scores_array = point.filter(['winter_mean'], axis=1)  # note that this score is >= 1
        # G() will only accept a flattened 1-D array. .Ravel() transforms a dictionary data-structure to np-array [x]
        # G(i): g = G_Local(heat_scores_array, w, transform='R')
        print("pass 6: calculating getis-ord z-scores")
        g = G_Local(heat_scores_array, w, transform='R', star=True)
        # Important Note: Large dimensions of w as args to G() will yield a MemoryError during processing.
        # Solution: will need to process the national footprint's w and heat_score matrices in segments!
        print("pass 7: outputting results to dictionary")
        # Append the z-scores and heat_scores back to point
        point["Z_Scores"] = g.z_sim
        point["Heat_Scores"] = heat_scores_array
        # Assign this iteration of g and heat_scores to the Outputs dictionary
        self.output_dict[region] = g, heat_scores_array, point
        if self.hex_template.crs != self.us_albers_equal_area:
            self.hex_template.crs = self.us_albers_equal_area
        if point.crs != self.us_albers_equal_area:
            point.crs = self.us_albers_equal_area
        hex_scores = geopandas.sjoin(self.hex_template, point, how='left')
        return hex_scores
