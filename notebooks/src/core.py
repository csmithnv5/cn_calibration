import numpy as np
from matplotlib import pyplot as plt
import pathlib as pl
import os
import csv
import geopandas as gpd
import rasterio
import pandas as pd
# import gdal
import requests
import urllib.parse
import urllib.request
import ssl
import json
from shapely.geometry import shape
from shapely.geometry import Polygon
import ipyleaflet as ipy
from IPython.display import clear_output
import numpy
from shapely.wkt import loads
import datetime
import fiona
import matplotlib.colors as colors
from src.modules.curve_number.utils_CN import *
import rioxarray
import xarray as xr

def df_to_gdf(df,geo_field,crs):
    df[geo_field] = df[geo_field].apply(lambda x :shape(x))
    gdf = gpd.GeoDataFrame(df).set_geometry(geo_field)
    gdf = gdf.set_crs(epsg=crs)
    return gdf

def df_to_gdf_polygon(df,geo_field,crs):
    df['geometry'] = df['geometry'].apply(lambda x :Polygon(x['rings'][0]))
    gdf = gpd.GeoDataFrame(df).set_geometry(geo_field)
    gdf = gdf.set_crs(epsg=crs)
    return gdf

def usgs_api_gage_to_df(site:str,crs) -> dict:
    '''
    queries the nld api server to gather feature data
    '''
    url = 'https://labs.waterdata.usgs.gov/api/observations/collections/monitoring-locations/items/USGS-'+site+'?'
    parameters = {'f':'json'}
    r = requests.get(url,params=parameters)
    response = json.loads(r.content)
    gdf = gpd.GeoDataFrame(response['properties'],index =[0])
    gdf['geometry'] = shape(response['geometry'])
    gdf = gdf.set_crs(epsg=crs)
    gdf['id'] = response['id']
    return gdf

def usgs_api_associated_gage_geometry(feature:str,site:str) -> dict:
    '''
    queries the nld api server to gather feature data
    '''
    url = 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-'+site+'/'+feature
    parameters = {'f':'json'}
    r = requests.get(url,params=parameters)
    response = json.loads(r.content)
    return response

def json_to_df_huc(features):
    df = pd.DataFrame(features,columns = features[0].keys())
    df2 = pd.json_normalize(df['attributes'])
    df3 = pd.concat([df.drop(['attributes'], axis=1), df2], axis=1)
    
    return df3

def json_to_df(features):
    df = pd.DataFrame(features,columns = features[0].keys())
    return df

def esri_rest_query(base_url,parameters):
    import json
    r = requests.get(base_url,parameters)
    try:
        response = json.loads(r.content.decode("utf-8"))
        return response
    except:
        print('service error')
        
def get_huc_12_bounds(huc12):
    '''
    takes a HUC12 and then queries USGS HUC-12 services
    to identify all the geometry extent of the
    HUC12 feature.
    Use WWF function if not in the US
    '''
    #intialize lists
    huc12s = []
    #current USGS wbd service
    nhd_service = 'https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6/query?' #base url
    #set parameters for query
    where_c_init = "huc12 = '"+ huc12+"'"
    nhd_param = {'outFields':'*','where':where_c_init,
                 'f':'json','returnExtentOnly':'true','outSR':4326}
    #get query results
    print('obtaining huc12 basin information')
    local_huc12 = esri_rest_query(nhd_service,nhd_param)
    
    if not local_huc12:
        return ['service error']
    if local_huc12['extent']:
        extent_list = list(local_huc12['extent'].values())[:4]
        bounds = ','.join([str(x)[:str(x).find('.')+7] for x in extent_list])
        return bounds
    else:
        print('Not in USA')
        return ['not in USA']

def gdf_of_local_usgs_gages(bounds):
    url = 'https://waterservices.usgs.gov/nwis/iv/?'
    parameters = {'format':'json', 'parameterCd':'00045,00060','siteStatus':'active','bBox':bounds}
    r = requests.get(url,params=parameters)
    response = json.loads(r.content)
    df = pd.json_normalize(response['value']['timeSeries'], sep='_')
    if df.empty is False:
        df['gage_type'] = df['variable_variableCode'].apply(lambda x: x[0]['value'])
        df['gage_number'] = df['sourceInfo_siteCode'].apply(lambda x: x[0]['value'])
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.sourceInfo_geoLocation_geogLocation_longitude, df.sourceInfo_geoLocation_geogLocation_latitude))
        return gdf
    else:
        print('No gages nearby')
        return df

def gdf_of_local_precip_gages(bounds,sr,search_distance):
    desired_stations = ['GHCN Daily']
    rain_gages = []
    base_url = 'https://gis.ncdc.noaa.gov/arcgis/rest/services/cdo/stations/MapServer/'
    parameters = {'geometry':bounds,'geometryType':'esriGeometryEnvelope','f':'pjson','inSR':sr,'outFields':'*','distance':search_distance,'units':'esriSRUnit_Meter','outSR':4269}
    for station in desired_stations:
        level = layer_indexer(base_url,{'f':'pjson'},station)
        results = esri_query(base_url, parameters, level)
        if 'features' in results.keys():
            if results['features']:
                for point in results['features']:
                    rain_gages.append(point)
    df = pd.DataFrame.from_records(rain_gages)
    if df.empty is False:
        print('Found gages nearby')
        df['STATION_ID'] = df['attributes'].apply(lambda x: x['STATION_ID'])
        df['STATION_NAME'] = df['attributes'].apply(lambda x: x['STATION_NAME'])
        df['DATA_BEGIN_DATE'] = df['attributes'].apply(lambda x: pd.to_datetime(x['DATA_BEGIN_DATE'],unit='ms'))
        df['DATA_END_DATE'] = df['attributes'].apply(lambda x: pd.to_datetime(x['DATA_END_DATE'],unit='ms'))
        df['lat'] = df['attributes'].apply(lambda x: x['LATITUDE'])
        df['long'] = df['attributes'].apply(lambda x: x['LONGITUDE'])
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.long, df.lat))
        gdf.set_crs(epsg=4269,inplace=True)
        return gdf
    else:
        print(f'No gages nearby at {search_distance} meter search distance')
        return df

def add_to_gage(df):
    df['Period of Record'] = df[['DATA_END_DATE','DATA_BEGIN_DATE']].apply(lambda x: (pd.to_datetime(x[0])-pd.to_datetime(x[1])).days/365.25,axis=1)
    df['link_works'] = df['STATION_ID'].apply(lambda x: "yes" if requests.get(f"https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/{x[x.find(':')+1:]}.csv").status_code == 200 else "no")
    return df

def download_file(url,out_dir):
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(out_dir/local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                #if chunk: 
                f.write(chunk)
    return out_dir/local_filename

def ncei_api_gage_to_df(bounding_box:str,crs) -> gpd.GeoDataFrame:
    '''
    queries the nld api server to gather feature data
    '''
    url = 'https://www.ncei.noaa.gov/access/services/data/v1?'
    parameters = {'format':'json', 'boundingBox':'49.795,-2.073,49.183,-0.992'}
    r = requests.get(url,params=parameters)
    response = json.loads(r.content)
    gdf = gpd.GeoDataFrame(response['properties'],index =[0])
    gdf['geometry'] = shape(response['geometry'])
    gdf = gdf.set_crs(epsg=crs)
    gdf['id'] = response['id']
    return gdf

def get_huc_12_gdf(huc12):
    '''
    takes a HUC12 and then queries USGS HUC-12 services
    to identify all the geometry extent of the
    HUC12 feature.
    Use WWF function if not in the US
    '''
    #intialize lists
    huc12s = []
    #current USGS wbd service
    nhd_service = 'https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6/query?' #base url
    #set parameters for query
    where_c_init = "huc12 = '"+ huc12+"'"
    nhd_param = {'outFields':'*','where':where_c_init,
                 'f':'json','returnGeometry':'true','outSR':4326}
    #get query results
    print('obtaining huc12 basin information')
    local_huc12 = esri_rest_query(nhd_service,nhd_param)
    
    if not local_huc12:
        return ['service error']
    if local_huc12['features']:
        df = json_to_df_huc(local_huc12['features'])
        gdf = df_to_gdf_polygon(df,'geometry',4326)
        return gdf
    
def get_huc_12_gdf_from_bigger_huc(huc):
    '''
    takes a HUC12 and then queries USGS HUC-12 services
    to identify all the geometry extent of the
    HUC12 feature.
    Use WWF function if not in the US
    '''
    #intialize lists
    huc12s = []
    #current USGS wbd service
    nhd_service = 'https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6/query?' #base url
    #set parameters for query
    where_c_init = "huc12 LIKE '"+ huc+"%'"
    nhd_param = {'outFields':'*','where':where_c_init,
                 'f':'json','returnGeometry':'true','outSR':4326}
    #get query results
    print('obtaining huc12 basin information')
    local_huc12 = esri_rest_query(nhd_service,nhd_param)
    
    if not local_huc12:
        return ['service error']
    if local_huc12['features']:
        df = json_to_df_huc(local_huc12['features'])
        gdf = df_to_gdf_polygon(df,'geometry',4326)
        return gdf