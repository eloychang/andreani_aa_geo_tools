import geopandas as gpd
import h3
from shapely.geometry import Polygon
import pandas as pd

def generate_hexagons_list(gdf, resolution=9):
    """Dada la geometria de las zonas con sus datos, devuelve un listado con los hexágonos pertenecientes a cada una.

    Args:
        gdf (GeoDataFrame): Geometria con sus campos de descripcion.
        resolution (int, optional): Resolución (tamaño) de hexágonos. Default: 9.

    Returns:
        DataFrame: Listado de hexágonos con su asignación de zona.
    """
    hexagons_list = [
        (hexagon, gdf['codigo'][index], gdf['idgla'][index], gdf['idsucursal'][index], gdf['zona'][index])
        for index, polylist in gdf.geometry.iteritems()
        for polygon in polylist.geoms
        for hexagon in h3.polyfill(polygon.__geo_interface__, resolution, geo_json_conformant=True)
    ]
    return pd.DataFrame(hexagons_list, columns=['h3_id', 'codigo', 'idgla', 'idsucursal', 'zona'])


def generate_branch_limits(df, filename = None):
    """Genera los bordes de las zonas delimitados por los hexágonos.

    Args:
        df (DataFrame): Listado de hexágonos con su asignación de zona.
        filename (String, optional): Nombre de archivo de salida en formato GeoJSON
    
    Returns:
        GeoDataFrame: Bordes de las zonas agrupados por código.
    """
    geometry = [Polygon(h3.h3_to_geo_boundary(h, geo_json=True)) for h in df['h3_id']]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
    gdf = gdf.dissolve(by='codigo').reset_index()
    if filename:
        gdf.to_file(f'{filename}.geojson', driver='GeoJSON')
    return gdf