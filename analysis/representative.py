""" Compare representativeness of WorldPop data to home locations for January 2022 based on Veraset data"""

import os 

import geopandas as gpd
import pandas as pd
import numpy as np

import folium as folium

from matplotlib import pyplot as plt

import rasterio 

import rasterio
from rasterio.features import shapes
import rasterio.plot
import geopandas as gpd
import pandas as pd
from shapely import wkt
from shapely.ops import unary_union

import scipy.stats

from pyspark import SparkContext
from pyspark.sql import SparkSession, SQLContext
from pyspark.sql import functions as F
from pyspark.sql.functions import pandas_udf, PandasUDFType, col, asc
from pyspark.sql.types import (
    StructType,
    StructField,
    DoubleType,
    IntegerType,
    StringType,
)
spark = SparkSession.builder.appName(f"spark").getOrCreate()



dir_main = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..'))
dir_data = os.path.join(dir_main, 'data', 'raw', 'geo', 'worldpop')

# set seaborn plot style
import seaborn as sns
sns.set_style("white")


def make_center_line_coords(df,x,y):
  # This creates two coords to define the hypothetical line in which x and y are perfectly correlatd
  expected_slope = df[y].sum()/df[x].sum()
  xcoord = np.array([df[x].min(), df[x].max()])
  ycoord = expected_slope*xcoord
  return((xcoord,ycoord))


def getHomesData():
    # load homes data
    homes = spark.read.parquet(os.path.join(dir_main, 'data', 'raw', 'homes_2022-01-01_2022-12-31'))
    homes = homes.filter(F.col("start_date") == "2022-01-02") # only month of january
    homes = homes.toPandas() # 3 million unique homes

    # remove sinkholes 
    sinkholes = homes.groupby(['longitude', 'latitude']).agg({'caid' : "count"}).reset_index()
    sinkholes = sinkholes[sinkholes.caid > 3]
    homes = homes[~homes.index.isin(sinkholes.index)]

    # join homes and raion shapefile
    homes_gdf = gpd.GeoDataFrame(homes, geometry=gpd.points_from_xy(homes.longitude, homes.latitude), crs = "epsg:4326")
    homes_gdf = homes_gdf.to_crs("epsg:2587")
    homes_gdf.drop(['longitude', 'latitude', 'start_date'], axis=1, inplace=True)
    
    return homes_gdf

def prepWorldPopData():
    # read wp data

    reread_wp = False

    if reread_wp:
        
        fn = os.path.join(dir_data, "ukr_ppp_2020_UNadj_constrained.tif")
        #fn = bsdir + "data/in/uae/are_ppp_2020_UNadj.tif"

        # load worldpop data from tif file
        wp = rasterio.open(fn)

        mask = None
        with rasterio.Env():
            with rasterio.open(fn) as src:
                image = src.read(1) # first band
                results = (
                {'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v) 
                in enumerate(
                    shapes(image, mask=mask, transform=src.transform)))

        geoms = list(results)
        # first feature
        geoms[0]
        wp = gpd.GeoDataFrame.from_features(geoms)
        # drop last row as it is a bounding box for the whole country
        wp = wp.drop(wp.tail(1).index)
        wp = wp.set_crs('epsg:4326')
        wp = wp.to_crs('epsg:2587')

        wp.to_file(os.path.join(dir_data, 'worldpop.shp'), driver='ESRI Shapefile')
    else:
        wp = gpd.read_file(os.path.join(dir_data, 'worldpop.shp'))

    # aggregate to raion
    raion = gpd.read_file(os.path.join(dir_main, 'data', 'raw', 'geo', 'ukr_admin_level2'))
    raion = raion[['geometry', 'name_2']]
    raion = raion.to_crs("epsg:2587")

    # do overlay and then only retain grid cells with largest overlap with raion
    wp['wp_index'] = wp.index
    wpr = gpd.overlay(raion, wp, how='intersection')
    # Sort by area so largest area is last
    wpr.sort_values(by='geometry', inplace=True, key =lambda col: np.array([x.area for x in col]))

    wpr.drop_duplicates(subset= "wp_index", keep='last', inplace=True)
    wpr.loc[wpr.raster_val < 0, 'raster_val'] = 0

    # aggregate to raion
    raipop = wpr.groupby(['name_2']).agg({'raster_val' : 'sum'}).reset_index()

    raipop = raipop.join(raion.set_index('name_2'), on='name_2')
    raipop = gpd.GeoDataFrame(raipop, geometry='geometry', crs = "epsg:2587")

    return homes_gdf.sjoin(raipop, how = "left", op = "intersects")





def makePlot(join, x_val = "raster_val", y_val = "caid",  clr = "#72aaa1"):
    fig, ax = plt.subplots()
    m, b, r_value, p_value, std_err = scipy.stats.linregress(join[y_val], join[x_val])
    join.plot(x=x_val, y=y_val, kind='scatter', ax = ax, color = clr, s = 50, alpha = 0.5)
    line_coords = make_center_line_coords(join,x_val,y_val)
    ax.plot(line_coords[0],line_coords[1], color = 'k',linestyle='-', linewidth=1)
    ax.annotate('R^2: ' + str("{:.2f}".format(r_value**2)) + ", Corr.: " +  str("{:.2f}".format(r_value)), 
                xy=(0.6, 0.95), xycoords='axes fraction', 
                # increase font size
                size = 12
                )
    # set x label
    ax.set_xlabel("Population (million)")
    ax.set_ylabel("Number of devices residing (million)")
    
    # make graph prettier
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # set major grid only, in light grey
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

    return ax
    
    
if __name__ == "__main__":
    
    
    ## load homes data
    homes_gdf = getHomesData()

    #--------------------------------------------
    #### WorldPop data
    #--------------------------------------------

    join = prepWorldPopData()


    #### raion-level plot
    join = join.groupby(['name_2']).agg({'raster_val' : 'mean', 'caid' : "count"}).reset_index()


    join[['caid', 'raster_val']] = join[['caid', 'raster_val']] / 1000000

    makePlot(clr = "#009392", join = join)
    plt.savefig(os.path.join(dir_main, 'output', 'figs', 'worldpop_vs_caid_raion.pdf'), dpi = 600, bbox_inches = 'tight')


    #### oblast level plot
    raion = gpd.read_file(os.path.join(dir_main, 'data', 'raw', 'geo', 'ukr_admin_level2'))
    raion = raion[['name_2', 'name_1']]

    join = join.merge(raion, left_on='name_2', right_on='name_2')

    joinoblast = join.groupby(['name_1']).agg({'raster_val' : 'sum', 'caid' : "sum"}).reset_index()

    makePlot(join = joinoblast, clr = "#e5b9ad")
    plt.savefig(os.path.join(dir_main, 'output', 'figs', 'worldpop_vs_caid_oblast.pdf'), dpi = 600, bbox_inches = 'tight')



    #--------------------------------------------
    #### COD-PS data
    #--------------------------------------------

    # get oblast "shapefile" (this is at admin level 4 but we can aggregate it)
    oblast = gpd.read_file(os.path.join(dir_main, 'data', 'raw', 'geo', 'ukr_admbnda_sspe_20230201_SHP'))
    oblast = oblast[['ADM1_PCODE', 'geometry']]

    # get oblast population
    oblastpop = pd.read_csv(os.path.join(dir_main, 'data', 'raw', 'geo', 'ukr_admpop_adm1_2022.csv'))

    oblastpop = oblastpop[['ADM1_PCODE', 'T_TL']]

    oblast = oblast.merge(oblastpop, left_on='ADM1_PCODE', right_on='ADM1_PCODE')
    oblast = oblast.to_crs("epsg:2587")

    join = homes_gdf.sjoin(oblast, how = "inner", op = "intersects")

    #### oblast level plot
    join = join.groupby(['ADM1_PCODE']).agg({'T_TL' : 'mean', 'caid' : "count"}).reset_index()
    
    join[['caid', 'T_TL']] = join[['caid', 'T_TL']] / 1000000

    
    makePlot(x_val = "T_TL", clr = "#d0587e", join = join)
    plt.savefig(os.path.join(dir_main, 'output', 'figs', 'codps_vs_caid_oblast.pdf'), dpi = 600, bbox_inches = 'tight')



