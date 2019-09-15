
import geopandas as gpd
import shapely
from shapely.geometry import Polygon, box, asPolygon, asPoint, MultiPoint, shape, base
import math
import shapefile

class Hexify:
    def __init__(self, region_dir, point_dir, output_dir):
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

        self.region_dir = region_dir
        self.region_shp = gpd.GeoDataFrame.from_file(region_dir)
        if self.region_shp.crs != self.us_albers_equal_area:
            print("converting crs of region layer")
            self.region_shp.to_crs(self.us_albers_equal_area, inplace=True)

        self.point_dir = point_dir
        self.points_shp = gpd.GeoDataFrame.from_file(point_dir)
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

    def geoms_to_hex_shp(self, in_geoms, out_shp, projection):
        # algorithm pulled from here:
        # https://github.com/DigitalGlobe/gbdxtools/blob/master/gbdxtools/catalog_search_aoi.py
        # pip install pyshp

        prj_name = '{}.prj'.format(out_shp.split('.')[0])
        with open(prj_name, 'w') as prj:
            prj.write(projection)
        shp_writer = shapefile.Writer(shapefile.POLYGON)

        out_fields = [
            ['id', 'N']
        ]
        # out_fields_names = [x[0] for x in out_fields]
        for name in out_fields:
            shp_writer.field(*name)

        for in_id, geom in enumerate(in_geoms, start=1):
            print(in_id)
            print(geom)
            shp_writer.record(*[str(in_id)])
            shp_writer.poly(parts=[geom])
        shp_writer.save(out_shp)

    # ------------------------------------------------------------------------------------------------
    def create_subpoints(self, points, region, us_albers_equal_area):
        """
        a spatial join function for points and a region shapefile
        :param points: point shapefile
        :param region: polygon shapefile
        :param us_albers_equal_area: projection
        :return: a spatially joined gdf
        """
        sub_gdf = points.copy()
        if sub_gdf.crs != us_albers_equal_area:
            print("converting crs of points layer")
            sub_gdf.to_crs(us_albers_equal_area, inplace=True)

        regions_gdf = region.copy()

        points_polgon_join = gpd.sjoin(sub_gdf, regions_gdf, how='left')

        return points_polgon_join

    # ------------------------------------------------------------------------------------------------
    def create_hex_grid_new(self, hex_area, overlap, saveout=False):
        """

        :param hex_area: should be around 5000-8000 meters squared
        :param saveout: if True, will save-out the hex-grid as a shp to output_dir
        :return: a hex-grid stored in a dictionary
        """

        hex_region_field = 'NAME'
        # hex_area = 86602540.378  # 86602540.378 is the area of a hexagon with radius = 5000 meters

        points_region = self.create_subpoints(self.points_shp, self.region_shp, self.us_albers_equal_area)
        s_join_n = points_region.dropna(subset=['NAME'])
        # For each gdf in self.sample_points_region, determine the bounding box extent, and generate a hex grid
        all_bboxes = s_join_n.groupby(hex_region_field)['geometry'].agg(lambda g: list(g.total_bounds))
        # gpd update: https://stackoverflow.com/questions/27439023/pandas-groupby-agg-function-does-not-reduce
        # https://stackoverflow.com/questions/39840546/must-produce-aggregated-value-i-swear-that-i-am

        hex_list = {}

        for i, b in enumerate(all_bboxes.iteritems()):
            name = b[0]
            print(name)
            geom = b[1]
            print(geom)
            geom_box = box(*geom)
            print(geom_box)
            extent_area = geom_box.area
            print(extent_area)

            if extent_area < hex_area:  # 86602540.378 is the area of a hexagon with radius = 5000 meters
                # begin rescaling the extent_area
                # (1) The following condition, box(a,b,c,d).area/box(w,s,e,n).area = xfact*yfact.
                # (2) want to know --> (86602540.378*1.1)/box(w,s,e,n).area --> sqrt(*)
                req_growth = (hex_area * 1.1) / extent_area
                scale_factor = req_growth ** 0.5
                nshp_bounds = shapely.affinity.scale(geom_box, xfact=scale_factor, yfact=scale_factor).bounds
                # obtain new bounds
                geom = nshp_bounds

            # # # Create the fishnet # # #
            w, s, e, n = geom
            vertices = self.calculate_polygons(w, s, e, n, hex_area)  # 5000 meter radius (UTM)

            # Create the hex geodataframe
            v_dict = {ndx: Polygon(v) for ndx, v in enumerate(vertices)}
            hex_gdf = gpd.GeoDataFrame.from_dict({'geometry': v_dict})
            hex_gdf.crs = s_join_n.crs
            hex_gdf = hex_gdf.reset_index().rename(columns={'index': 'polyid'})
            hex_gdf = hex_gdf[
                hex_gdf.intersects(self.region_shp.geometry[self.region_shp[self.region_shp['NAME'] == name].index[0]])]
            # print(hex_gdf)

            # store in a zipped list (file_name will keep track of which hex-region is being stored)
            hex_list[name] = hex_gdf
            print("Generating hex grid for: %s" % name)
            print("================================")
            # write hex to shapefile for debugging
            if saveout == True:
                self.geoms_to_hex_shp([list(s.exterior.coords) for s in hex_list[name]['geometry'].tolist()], \
                                      '{}/{}_{}_{}p_Hex.shp'.format(self.output_dir, name, hex_area, overlap),
                                      self.us_albers_equal_area_prj)
            del hex_gdf
        return hex_list
