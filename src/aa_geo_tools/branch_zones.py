import geopandas as gpd
import h3
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.ops import cascaded_union, transform
import pyproj

# Proyecciones para calcular las distancias en metros
source_crs = pyproj.CRS("EPSG:4326")
target_crs = pyproj.CRS("EPSG:3857")
project = pyproj.Transformer.from_crs(
    source_crs, target_crs, always_xy=True).transform

def update_cp_and_zone(row, df):
    """Actualiza cp y zona de hexágonos duplicados, tomando las del borde del polígono más lejano al centroide del hexágono."""
    grouped = df.groupby('h3_id').get_group(row['h3_id'])
    
    if len(grouped) == 2:
        proyected_poly1 = transform(project, grouped.iloc[0]['geometry'])
        proyected_poly2 = transform(project, grouped.iloc[1]['geometry'])
        
        # Obtengo centroide
        lat, lon = h3.h3_to_geo(row['h3_id'])
        # Proyecto para poder obtener distancia en metros
        projected_point = transform(project, Point(lon, lat))
        
        distance_to_poly1 = proyected_poly1.distance(projected_point)
        distance_to_poly2 = proyected_poly2.distance(projected_point)
        
        if distance_to_poly1 > distance_to_poly2:
            codigo = grouped.iloc[0]['codigo']
            idgla = grouped.iloc[0]['idgla']
            idsucursal = grouped.iloc[0]['idsucursal']
            zona = grouped.iloc[0]['zona']
        else:
            codigo = grouped.iloc[1]['codigo']
            idgla = grouped.iloc[1]['idgla']
            idsucursal = grouped.iloc[1]['idsucursal']
            zona = grouped.iloc[1]['zona']
    else:
        codigo = row['codigo']
        idgla = row['idgla']
        idsucursal = row['idsucursal']
        zona = row['zona']
    
    updated_row = row.copy()
    updated_row['codigo'] = codigo
    updated_row['idgla'] = idgla
    updated_row['idsucursal'] = idsucursal
    updated_row['zona'] = zona
    
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
    list_polygons_proj = [transform(project, poly) for poly in list_polygons]

    # Asigno borde a cada zona dependiendo de la cercanía con el centroide
    list_index_cp_zone = []
    for h in set_hex_border:
        # Obtengo centroide
        lat, lon = h3.h3_to_geo(h)
        # Proyecto para poder obtener distancia en metros
        projected_point = transform(project, Point(lon, lat))

        min_distance = float('inf')
        min_index = 0
        # Recorro los polígonos para ver el más cercano
        for i, projected_polygon in enumerate(list_polygons_proj):
            distance_to_polygon = projected_polygon.distance(projected_point)
            if distance_to_polygon < min_distance:
                min_distance = distance_to_polygon
                min_index = i
        list_index_cp_zone.append(min_index)

    # Genero una lista y asigno al conjunto mediante el índice del más cercano
    list_border = [set() for _ in range(len(list_polygons))]
    for i, h in enumerate(set_hex_border):
        idx_in_gdf = list_index_cp_zone[i]
        list_border[idx_in_gdf].add(h)
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

def get_hexes_traversed_by_borders(polygon, res):
    """Obtiene los hexágonos del borde del polígono, en la resolución dada.

    Args:
        polygon (geometry): polígono del área
        res (int): resolución de los hexágonos

    Returns:
        set(string): conjunto de hexágonos obtenidos
    """
    set_traversed_hexes = set()
    coords = polygon.boundary.coords
    for j in range(len(coords)-1):
        # Para cada segmento del linestring, obtengo el punto inicial y final
        start_leg = coords[j]
        stop_leg = coords[j+1]
        # Obtengo el hexágono de cada coordenada (están como lon/lat)
        start_hexid = h3.geo_to_h3(lat=start_leg[1],
                                   lng=start_leg[0],
                                   resolution=res)
        stop_hexid = h3.geo_to_h3(lat=stop_leg[1],
                                  lng=stop_leg[0],
                                  resolution=res)
        # Genero una línea como sucesión de hexágonos entre el inicial y final
        traversed_hexes = h3.h3_line(start=start_hexid,
                                     end=stop_hexid)
        set_traversed_hexes |= set(traversed_hexes)
    return set_traversed_hexes

def include_border_into_polygon(gdf):
    """Tomando todos los hexágonos del borde de los polígonos, 
    incluyo los que estén dentro de las áreas

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
    gdf['hex_filled_voids'] = gdf.apply(
        lambda row: row['hex_filled_voids'] | row['add_in_zone'], axis=1)
    return gdf[['codigo', 'idgla', 'idsucursal', 'zona', 'geometry', 'hex_filled_voids', 'hex_border']].copy()

def generate_hexagons_list(gdf, resolution=9):
    """Dada la geometría de las zonas con sus datos, devuelve un listado con los hexágonos pertenecientes a cada una.

    Args:
        gdf (GeoDataFrame): Geometría con sus campos de descripcion.
        resolution (int, optional): Resolución (tamaño) de hexágonos. Default: 9.

    Returns:
        DataFrame: Listado de hexágonos con su asignación de zona.
    """ 
    gdf['hex_border'] = None
    gdf['hex_filled_voids'] = None
    
    for i, row in gdf.iterrows():
        geometry = row['geometry']
        
        if isinstance(geometry, Polygon):
            polygons = [geometry]
        elif isinstance(geometry, MultiPolygon):
            polygons = geometry.geoms
        else:
            raise ValueError("Invalid geometry type. Expected Polygon or MultiPolygon.")
        
        hexagons_list = []
        borders_list = []
        for polygon in polygons:
            hexagons_list.extend(h3.polyfill(polygon.__geo_interface__, resolution, geo_json_conformant=True))
            borders_list.extend(get_hexes_traversed_by_borders(polygon, 9))

        gdf.at[i,'hex_border'] = borders_list
        gdf.at[i, 'hex_filled_voids'] = hexagons_list
        
    gdf = include_border_into_polygon(gdf)
    set_hex_border = get_set_border_not_in_polygon(gdf)
    gdf['hex_border'] = order_by_nearest_polygon(set_hex_border, gdf['geometry'])
    gdf['h3_id'] = gdf.apply(lambda row: row['hex_border'] | row['hex_filled_voids'], axis=1)
    df = gdf.explode('h3_id').reset_index(drop=True)
    df = df.apply(update_cp_and_zone, args=(df,), axis=1).drop_duplicates(subset='h3_id')
    return df[['h3_id', 'codigo', 'idgla', 'idsucursal', 'zona']].reset_index(drop=True).copy()


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