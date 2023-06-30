import geopandas as gpd
import h3
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, Point, mapping, LineString
from shapely.ops import unary_union, transform
import pyproj
import warnings
warnings.filterwarnings("ignore", message="CRS not set for some of the concatenation inputs.")

# Proyecciones para calcular las distancias en metros
source_crs = pyproj.CRS("EPSG:4326")
target_crs = pyproj.CRS("EPSG:3857")
projection = pyproj.Transformer.from_crs(
    source_crs, target_crs, always_xy=True).transform

def update_zone(row_group):
    """Actualiza datos de hexágonos duplicados, quedándose con las del borde del polígono más lejano al centroide del hexágono."""
    proyected_poly1 = transform(projection, row_group.iloc[0]['geometry'])
    proyected_poly2 = transform(projection, row_group.iloc[1]['geometry'])
    
    # Obtengo centroide
    lat, lon = h3.h3_to_geo(row_group.iloc[0]['h3_id'])
    # Proyecto para poder obtener distancia en metros
    projected_point = transform(projection, Point(lon, lat))
    
    distances = [
        proyected_poly1.distance(projected_point),
        proyected_poly2.distance(projected_point)
    ]
    
    updated_row = row_group.iloc[0].copy()
    
    if distances[0] > distances[1]:
        updated_row[['codigo', 'idgla', 'idsucursal', 'zona']] = row_group.iloc[0][['codigo', 'idgla', 'idsucursal', 'zona']]
    else:
        updated_row[['codigo', 'idgla', 'idsucursal', 'zona']] = row_group.iloc[1][['codigo', 'idgla', 'idsucursal', 'zona']]
    
    return updated_row


def order_by_nearest_polygon(set_hex_border, list_polygons):
    """Divido los hexágonos del borde por conjuntos, asignandoles el polígono más cercano.

    Args:
        set_hex_border (set(string)): conjunto de hexágonos del borde de los polígonos
        list_polygons (list(geometry)): lista de las geometrías de los polígonos

    Returns:
        list(set(string)): lista de conjuntos de hexágonos,
        donde el índice en la lista indica la zona a la que pertenecen
    """
    # Proyecto todos los polígonos
    list_polygons_proj = [transform(projection, poly) for poly in list_polygons]

    # Asigno borde a cada zona dependiendo de la cercanía con el centroide
    list_border = [set() for _ in range(len(list_polygons))]
    for h in set_hex_border:
        # Obtengo centroide
        lat, lon = h3.h3_to_geo(h)
        # Proyecto para poder obtener distancia en metros
        projected_point = transform(projection, Point(lon, lat))

        min_distance = float('inf')
        min_index = 0
        # Recorro los polígonos para ver el más cercano
        for i, projected_polygon in enumerate(list_polygons_proj):
            distance_to_polygon = projected_polygon.distance(projected_point)
            if distance_to_polygon < min_distance:
                min_distance = distance_to_polygon
                min_index = i
        list_border[min_index].add(h)

    return list_border


def get_set_border_not_in_polygon(gdf):
    """Tomando todos los hexágonos del borde de los polígonos, 
    obtengo los que no estén incluidos dentro de alguna de las áreas.

    Args:
        gdf (GeoDataFrame): debe contener la columna 'hex_border' con los bordes de
        cada área, y la columna 'geometry' con la geometría de los polígonos

    Returns:
        set(string): conjunto de los hexágonos del borde
    """
    # Acumulo todos los hexágonos de los bordes
    set_hex_border = {h for row in gdf['hex_border'] for h in row}

    set_hex_border_copy = set_hex_border.copy()
    polygons = gdf['geometry']
    for h in set_hex_border:
        lat, lon = h3.h3_to_geo(h)
        for poly in polygons:
            # Elimino el hexágono del borde si está dentro de una zona
            if poly.contains(Point(lon, lat)):
                set_hex_border_copy.discard(h)
    return set_hex_border_copy


def interpolate_coords_between(start, end, num_points):
    # Create a LineString object between the start and stop points
    segment = LineString([start, end])
    interpolated_coords = [start]
    for i in range(num_points):
        fraction = (i + 1) / (num_points + 1)
        point = segment.interpolate(fraction, normalized=True)
        interpolated_coords.append((point.x, point.y))
    interpolated_coords.append(end)
    return interpolated_coords


def get_hexes_traversed_by_borders(polygon, res):
    """Obtiene los hexágonos del borde del polígono, en la resolución dada.

    Args:
        polygon (geometry): polígono del área
        res (int): resolución de los hexágonos

    Returns:
        set(string): conjunto de hexágonos obtenidos
    """
    set_traversed_hexes = set()

    # Obtengo coordenadas del borde del polígono
    coords = polygon.boundary.coords
    for j in range(len(coords)-1):
        # Para cada segmento del linestring, obtengo el punto inicial y final
        start_leg = coords[j]
        stop_leg = coords[j+1]
        interpolated = interpolate_coords_between(start_leg, stop_leg, 4)
        for start_i, end_i in zip(interpolated[:-1], interpolated[1:]):
            # Obtengo el hexágono de cada coordenada (están como lon/lat)
            start_hexid = h3.geo_to_h3(lat=start_i[1],
                                    lng=start_i[0],
                                    resolution=res)
            stop_hexid = h3.geo_to_h3(lat=end_i[1],
                                    lng=end_i[0],
                                    resolution=res)
            # Genero una línea como sucesión de hexágonos entre el inicial y final
            traversed_hexes = h3.h3_line(start=start_hexid,
                                        end=stop_hexid)
            set_traversed_hexes |= set(traversed_hexes)
            
    return set_traversed_hexes


def include_border_into_polygon(gdf):
    """Tomando los hexágonos del borde de los polígonos, incluyo los que estén dentro de las áreas

    Args:
        gdf (GeoDataFrame): debe contener la columna 'hex_border' con los bordes de
        cada área, y la columna 'geometry' con la geometría de los polígonos

    Returns:
        GeoDataFrame: GeoDataFrame agregando los hexágonos a las zonas
    """
    # Acumulo todos los hexágonos de los bordes
    set_hex_border = {h for row in gdf['hex_border'] for h in row}

    list_add_hex = []
    for _, row in gdf.iterrows():
        set_hex = set()
        for h in set_hex_border:
            lat, lon = h3.h3_to_geo(h)
            # Agrego el hexágono del borde si está dentro de una zona
            if row['geometry'].contains(Point(lon, lat)):
                set_hex.add(h)
        list_add_hex.append(set_hex)

    gdf['add_in_zone'] = list_add_hex
    gdf['hex_filled'] = gdf.apply(
        lambda row: set(row['hex_filled']) | row['add_in_zone'], axis=1)
    return gdf[['codigo', 'idgla', 'idsucursal', 'zona', 'geometry', 'hex_filled', 'hex_border']].copy()


def fill_hex_res(gdf, res):
    """
    Rellena los polígonos con hexágonos de resolución dada y delimita los bordes con hexágonos
    
    Args:
        gdf (GeoDataFrame): polígonos de las áreas con sus datos
        res (int): resolución (tamaño) de los hexágonos

    Returns:
        GeoDataFrame: mismo que el input, con los hexágonos resultado de rellenar las zonas y sus bordes.
    """
    gdf["geom_geojson"] = gdf["geometry"].apply(
        lambda x: mapping(x))
    gdf["hex_filled"] = gdf["geom_geojson"].apply(
        lambda x: list(h3.polyfill(geojson=x, res=res, geo_json_conformant=True)))
    gdf['hex_border'] = gdf.apply(lambda row : get_hexes_traversed_by_borders(row['geometry'], res), axis = 1)
    return gdf[['codigo', 'idgla', 'idsucursal', 'zona', 'hex_filled', 'hex_border', 'geometry']].copy()


def generate_hexagons_list(gdf, res):
    """Dada la geometría de las zonas con sus datos, devuelve un listado con los hexágonos pertenecientes a cada una.

    Args:
        gdf (GeoDataFrame): geometría con sus campos de descripcion.
        resolution (int, optional): resolución (tamaño) de hexágonos

    Returns:
        DataFrame: listado de hexágonos con su asignación de zona.
    """ 
    gdf = fill_hex_res(gdf, res)
    # Genero la asignación de los bordes    
    gdf = include_border_into_polygon(gdf)
    set_hex_border = get_set_border_not_in_polygon(gdf)
    gdf['hex_border'] = order_by_nearest_polygon(set_hex_border, gdf['geometry'])
    gdf['h3_id'] = gdf.apply(lambda row: row['hex_border'] | row['hex_filled'], axis=1)
    # Listado de hexágonos
    df = gdf.explode('h3_id').reset_index(drop=True)
    # Asigno zona más cercana para eliminar duplicados
    duplicated_rows = df.duplicated(subset=['h3_id'], keep=False)
    selected_rows = df[duplicated_rows].groupby('h3_id').apply(update_zone)
    result = pd.concat([df[~duplicated_rows], selected_rows])
    return result[['h3_id', 'codigo', 'idgla', 'idsucursal', 'zona']].reset_index(drop=True).copy()

def generate_branch_limits(df, filename=None):
    """Genera los bordes de las zonas delimitados con hexágonos.

    Args:
        df (DataFrame): listado de hexágonos con su asignación.
        filename (String, optional): nombre de archivo para guardarse como GeoJSON 
    
    Returns:
        GeoDataFrame: bordes de zonas agrupados por código.
    """
    def fill_voids(geom):
        if not isinstance(geom, MultiPolygon):
            return geom
        polygons = geom.geoms
        return unary_union(polygons)

    geometry = [Polygon(h3.h3_to_geo_boundary(h, geo_json=True)) for h in df['h3_id']]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
    gdf = gdf.dissolve(by='codigo', aggfunc='first').reset_index()
    gdf = gdf.to_crs("EPSG:3857") 
    gdf['geometry'] = gdf['geometry'].buffer(0.1) 
    gdf['geometry'] = gdf['geometry'].apply(fill_voids)
    gdf = gdf.to_crs("EPSG:4326")
    result = gdf[['codigo', 'idgla', 'idsucursal', 'zona', 'geometry']].copy()
    
    if filename:
        result.to_file(filename, driver='GeoJSON')
        
    return result
