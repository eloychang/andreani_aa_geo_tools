import geopandas as gpd
import h3
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import cascaded_union
import pandas as pd

def generate_hexagons_list(gdf, resolution=9):
    """Dada la geometría de las zonas con sus datos, devuelve un listado con los hexágonos pertenecientes a cada una.

    Args:
        gdf (GeoDataFrame): Geometría con sus campos de descripcion.
        resolution (int, optional): Resolución (tamaño) de hexágonos. Default: 9.

    Returns:
        DataFrame: Listado de hexágonos con su asignación de zona.
    """
    hexagons_list = []
    
    for _, row in gdf.iterrows():
        geometry = row['geometry']
        
        if isinstance(geometry, Polygon):
            polygons = [geometry]
        elif isinstance(geometry, MultiPolygon):
            polygons = geometry.geoms
        else:
            raise ValueError("Invalid geometry type. Expected Polygon or MultiPolygon.")
        
        for polygon in polygons:
            hexagons = h3.polyfill(polygon.__geo_interface__, resolution, geo_json_conformant=True)
            hexagons_list.extend([(hexagon, row['codigo'], row['idgla'], row['idsucursal'], row['zona']) for hexagon in hexagons])
    
    return pd.DataFrame(hexagons_list, columns=['h3_id', 'codigo', 'idgla', 'idsucursal', 'zona'])


def fill_voids(geom):
    """Rellena los huecos del MultiPolygon y lo convierte en Polygon si es necesario."""
    if not isinstance(geom, MultiPolygon):
        return geom
    
    polygons = list(geom)
    return cascaded_union(polygons)

def generate_branch_limits(df, filename=None):
    """Genera los bordes de las zonas delimitados por los hexágonos.

    Args:
        df (DataFrame): Listado de hexágonos con su asignación de zona.
        filename (String, optional): Nombre de archivo de salida en formato GeoJSON
    
    Returns:
        GeoDataFrame: Bordes de las zonas agrupados por código.
    """
    geometry = [Polygon(h3.h3_to_geo_boundary(h, geo_json=True)) for h in df['h3_id']]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
    gdf = gdf.dissolve(by='codigo', aggfunc='first').reset_index()
    gdf['geometry'] = gdf['geometry'].apply(fill_voids)
    
    if filename:
        gdf.to_file(f'{filename}.geojson', driver='GeoJSON')
        
    return gdf