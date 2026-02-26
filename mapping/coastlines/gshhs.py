"""
Modified from: https://github.com/aragong/coastline-loader
Made easier to use as a script, rather than requiring a whole package

Takes GSHHG Shapefile coastline data and exports coastline data in a given
lat lon box to a GeoPandas dataframe so that it can be plotted using external
tools and not requiring something like PyGMT. Makes it easier to make maps
outside of dedicated software.

Data Hosted At: https://www.soest.hawaii.edu/pwessel/gshhg/
Download Link: http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip

.. note:: dependencies
    - pandas
    - geopandas
    - shapely
    - matplotlib
"""
import os
import sys
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

# SHP_DATABASE_DIR = os.environ.get("GSHHG_SHP_DATABASE_PATH")
SHP_DATABASE_DIR = "/Users/chow/Data/cartography/coastlines/gshhg-shp-2.3.7"


class GetCoastline:
    def __init__(
        self,
        resolution: str = "f",
        lonlatbox: tuple = None,
        output_epsg: int = None,
        layer: str = "L1",
    ) -> None:
        """
        Load GSHHS coastline database to GeoDataFrame.
        Clip and convert coordinates to required lonlatbox and epsg code, 
        if required.

        Args:
            resolution (str, optional): GSHHS resolution 
                "f", "h", "i", "l" or "c" (full, high, intermediate, low 
                or crude). Defaults to "c".
            lonlatbox (tuple, optional): box coordinates, (lonmin, lonmax, 
                latmin, latmax). Defaults to None.
            output_epsg (int, optional): epsg-code for output transformation. 
                Defaults to None.
            layer (str, optional): layer from GSHHS database. Defaults to "L1".
        """
        assert(resolution  in ["c", "f", "h", "i", "l"]), "incorrect resolution"
        self.gdf = self._load_from_shp(resolution, lonlatbox, 
                                       output_epsg, layer)

    def _load_from_shp(
        self, resolution: str, lonlatbox: tuple, output_epsg: int, layer: str
    ) -> gpd.GeoDataFrame:
        """
        Load and clip GSHHS shapefile (layer "L1") and convert to required 
        epsg, if needed

        Args:
            resolution (str): GSHHS database resolution
            lonlatbox (tuple): box coordinates, (lonmin, lonmax, latmin, latmax)
            output_epsg (int): epsg-code for output transformation
            layer (str): GSHHS database layer

        Returns:
            gpd.GeoDataFrame: costlines polygons
        """
        if lonlatbox:
            polygon = Polygon(
                [(lonlatbox[0], lonlatbox[2]),
                 (lonlatbox[0], lonlatbox[3]),
                 (lonlatbox[1], lonlatbox[3]),
                 (lonlatbox[1], lonlatbox[2]),
                 (lonlatbox[0], lonlatbox[2]),
                 ]
            )
        path = os.path.join(SHP_DATABASE_DIR, "GSHHS_shp", resolution, 
                            f"GSHHS_{resolution}_{layer}.shp")
        if polygon:
            gdf = gpd.GeoDataFrame.from_file(
                path, bbox=(lonlatbox[0], lonlatbox[2], 
                            lonlatbox[1], lonlatbox[3])
            ).set_crs(4326)
            gdf = gdf.clip(polygon).explode(index_parts=False)
        else:
            gdf = gpd.GeoDataFrame.from_file(path).set_crs(4326)

        if output_epsg:
            gdf = gdf.to_crs(output_epsg)

        return gdf.reset_index(drop=True)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert GeoDataFrame to standard DataFrame

        Returns:
            pd.DataFrame: costlines polygons
        """
        df = pd.DataFrame([])
        dataframes = []
        for index, (feature, area) in enumerate(zip(self.gdf.geometry, 
                                                    self.gdf.area)):
            lon, lat = feature.exterior.coords.xy
            poly_id = [index] * len(feature.exterior.coords)
            area = [area] * len(feature.exterior.coords)
            tmp = pd.DataFrame(
                {"polygon_id": poly_id, "longitude": lon, 
                 "latitude": lat, "area": area}
            )
            dataframes.append(tmp)
        df = pd.concat(dataframes, ignore_index=True)
        return df.reset_index(drop=True)

def main():
    """
    Main processing function, following coastline-loader example

    Args
        lonlatbox (list): [lon_min, lon_max, lat_min, lat_max]
        resolution (str): "f", "h", "i", "l" or "c" (full, high, intermediate, 
            low or crude). Defaults to "c".
    """
    lonlatbox = [128, 130, 40, 43]
    resolution = "f"
    output_epsg = None  # UTM52N = 32562
    output_file = (f"coastline_{lonlatbox[0]}_{lonlatbox[1]}_"
                   f"{lonlatbox[2]}_{lonlatbox[3]}.csv")
    plot = True 

    lonlatbox = [float(_) for _ in lonlatbox]
    coast = GetCoastline(resolution=resolution, lonlatbox=lonlatbox,
                         output_epsg=output_epsg)
    coast.gdf.head()

    # if plot:
    #     coast.gdf.plot()
    #     plt.show()

    df = coast.to_dataframe()
    df.head()

    if plot:
        f, ax = plt.subplots()
        for i in df["polygon_id"].unique():
            df[df["polygon_id"] == i].plot(x="longitude", y="latitude", 
                                           ax=ax, legend=False)
        plt.show()

    df.to_csv(output_file)

if __name__ == "__main__":
    main()
