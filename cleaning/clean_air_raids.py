import os 

import geopandas as gpd
import pandas as pd
from bs4 import BeautifulSoup
import pytz
import folium as folium
import io
import json
import codecs




def storeSoupValuesAndText(soup):
    # store all values and text in a dictionary
    # remove "${indent}" from text
    d = {}
    for i in soup.find_all('option'):
        for val in i.get('value').split(", "):
            d[val] = i.text.replace("${indent}", "")
    return d

def concatDfByColumn(df1, df2, col):
    # concat two dataframes by column
    # df1 and df2 must have the same number of rows
    # col is the column to concat by
    # returns a dataframe
    df1 = df1.reset_index(drop = True)
    df2 = df2.reset_index(drop = True)
    return pd.concat([df1, df2[col]], axis = 1)

# convert to Kyiv datetime
def convertUTC(col): 
    ukr = pytz.timezone("Europe/Kiev")
    return pd.to_datetime(airlong[col], unit = "s", utc = True).dt.tz_convert(ukr).dt.tz_localize(None)



# get shapefile
def getAndMergeShapefile(fn, airlong, name_col, names_col = "VARNAME_1"):
    # fn is the filename of the shapefile
    # returns a geodataframe
    path = os.path.join(dir_data, 'raw', 'geo', fn)
    gdf = gpd.read_file(path)

    # split VARNAME_1 into separate columns by splitting on the vertical bar
    temp = gdf[names_col].str.split("|", expand = True)
    I = temp.shape[1]
    gdf = pd.concat([gdf, temp], axis = 1)

    airlong = airlong.merge(gdf[[name_col]], left_on = 'region_english', right_on = name_col, how = 'left')
    for i in range(0, I):
        # merge airlong with gdf on column i of VARNAME_1
        # if there is a match, assign ID_1 value to airlong
        temp = airlong.drop(columns = name_col).merge(gdf[[name_col, i]], left_on = 'region_english', right_on = i, how = 'left')
        temp = temp[temp[name_col].notnull()]
        temp = temp[[name_col]]
        airlong.update(temp)
    return airlong
    

def readGdf(fn):
    path = os.path.join(dir_data, 'raw', 'geo', fn)
    return gpd.read_file(path)


# check whether the pyspark dataframe column 'minute' is between the columns start_dt and end_dt of the pandas dataframe airlong
def checkMinute(minute, start_dt, end_dt):
    # minute is a pyspark dataframe column
    # start_dt and end_dt are pandas dataframe columns
    # returns a pyspark dataframe column
    return (minute >= start_dt) & (minute <= end_dt)


def plotHeatMap(airlong):
    """ plot heat map with days on the x axis, regions on the y axis, and the number of alerts as the color scale """
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    
    # duplicate rows for each alert for every day between start_dt and end_dt
    temp = airlong[['region_english', 'start_dt', 'end_dt']].copy()
    temp['start_dt'] = pd.to_datetime(temp['start_dt'])
    temp['end_dt'] = pd.to_datetime(temp['end_dt'])
    temp['days'] = temp.apply(lambda x: pd.date_range(x['start_dt'], x['end_dt']), axis = 1)
    temp = temp.explode('days')
    temp.drop(columns = ['start_dt', 'end_dt'], inplace = True)
    temp['days'] = temp['days'].dt.date
    temp = temp.groupby(['region_english', 'days']).size().reset_index()
    temp.columns = ['region_english', 'days', 'alerts']
    
    max_alerts = temp.alerts.max()

    
    # pivot table
    temp = temp.pivot(index = 'region_english', columns = 'days', values = 'alerts') 

    temp = temp.fillna(0)
    # sort by number of alerts
    temp = temp.reindex(temp.sum(axis = 1).sort_values(ascending = False).index)
    # plot heatmap
    fig, ax = plt.subplots(figsize = (20, 20))


    # create a list of the colors in the preceding comments
    myColors = ['#002E34', '#003E4B', '#004B61', '#015576', '#045D8B', '#07629F', '#2A7CAD', '#4D95BB', '#70ADC9', '#93C3D6', '#B6D8E4', '#D9ECF1', '#FDFEFE']
    myColors.reverse()
    #myColors = ["#f7feae","#b7e6a5","#7ccba2","#46aea0","#089099","#00718b","#045275"]
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, max_alerts)

    # sort temp by average number of alerts
    temp = temp.reindex(temp.mean(axis = 1).sort_values(ascending = False).index)
    #temp = temp.reindex(sorted(temp.index), axis = 0)
    
    sns.heatmap(temp, cmap = cmap, ax = ax, linecolor = "lightgray", linewidths = 0.05, 
                cbar_kws = {'label': 'Number of alerts', "shrink" : 0.5})
    ax.set_title("")
    
    # drop axis labels
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # increase font size
    ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 14)
    
    # increase legend font size
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 14)
    
    # add a tick for each number
    cbar.set_ticks([i for i in range(0, max_alerts + 1)])
    


    # increase legend title size
    cbar.ax.set_ylabel(cbar.ax.get_ylabel(), fontsize = 16)

    
    # set background white
    ax.set_facecolor("white")
    
    # tight layout
    plt.tight_layout()

    plt.savefig(os.path.join(dir_main, "output", "figs", "raid_heatmap.pdf"), dpi = 600)
    plt.show()
    
    
## produce shapefile of cities
dir_main = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..'))
dir_data = os.path.join(dir_main, 'data')
    




if __name__ == '__main__':
    
    fn = os.path.join(dir_data, "raw", "geo", "ukr_admbnda_sspe_20230201_SHP", "ukr_admbnda_adm3_sspe_20230201.shp")
    gdf = gpd.read_file(fn)
    
    # export city polygons for referee 1 question 2
    path = os.path.join(dir_data, 'raw', 'geo', 'ukr_admin_level2')
    gdf = gpd.read_file(path)
    cities = gdf[gdf.type_2.isin(["Misto", "Mis'ka Rada", None])]
    cities = cities[['geometry', 'name_2']].rename(columns = {'name_2': 'city'})
    
    # add in Kyiv and Sebastopol, which are classified as level 1 
    path = os.path.join(dir_data, 'raw', 'geo', 'ukr_admin_level1')
    gdf = gpd.read_file(path)
    gdf = gdf[gdf.NAME_1.isin(["Kiev City", "Sevastopol'"])][['geometry', 'NAME_1']].rename(columns = {'NAME_1': 'city'})
    cities = pd.concat([cities.to_crs('epsg:4326'), gdf.to_crs('epsg:4326')], axis = 0)
    
    cities.to_file(os.path.join(dir_data, 'raw', 'geo', 'cities'))
    
    # get cities with pop >30,000 from OSM and check which ones aren't in the shapefile
    # result: only 4 actual cities aren't, and all of them have a population <250,000
    fn = 'ukr_cities'
    path = os.path.join(dir_data, 'raw', 'geo', fn)
    gdf = gpd.read_file(path)
    gdf = gdf.sjoin(cities.to_crs('epsg:4326'), how = "left", op = "intersects")
    



    #--------
    # load json and convert it to the csv format Austin originally converted it to
    # get air raids data
    fn = os.path.join(dir_data, 'raw', 'airraid_alerts', 'sirens@106.json')
    with codecs.open(fn, 'r', 'utf-8-sig') as f:
        test = json.load(f)
    # convert json to dataframe
    df = pd.DataFrame(test)
    df = df.explode('locations')
    df = df.reset_index()

    # create a column that is an increasing index for each index value
    df['index2'] = df.groupby('index').cumcount()

    # convert df to wide format with each column a location of a region
    df = (df.pivot(index = "index2", columns = 'index', values = 'locations')
            .reset_index(drop = True).add_prefix("locations__"))

    air = df.copy()
    air = air.rename_axis(None, axis = 1) # remove index name
    del(df)


    test = pd.read_csv(os.path.join(dir_data, 'raw', 'airraid_alerts', 'sirens@99.csv'))

    # parse html list mapping ukrainian region names to english region names
    with open(os.path.join(dir_data, 'raw', 'ukr_to_eng.html')) as f:
        soup = BeautifulSoup(f, "html.parser")


    d = storeSoupValuesAndText(soup)

    # code even time values as start of raid and uneven as end of raid
    temp = pd.melt(air, value_vars = air.columns)
    temp = temp[temp.value.notnull()]
    start = temp.iloc[::2].rename(columns = {"value": "start"})
    end = temp.iloc[1::2].rename(columns = {"value": "end"})


    airlong = concatDfByColumn(start, end, "end")


    airlong['start_dt'] = convertUTC('start')
    airlong['end_dt'] = convertUTC('end')

    col = "start"
    ukr = pytz.timezone("Europe/Kiev")
    pd.to_datetime(airlong[col], unit = "s", utc = True).dt.tz_convert(ukr).dt.tz_localize(None)

    # add english names
    eng_cols = [d.get(i) for i in airlong.variable.str.replace("locations__", "").tolist()]
    airlong['region_english'] = eng_cols
    airlong.rename(columns = {"variable": "region_ukrainian"}, inplace = True)

    airlong[airlong['region_english'].notnull()].shape

    # write to file
    airlong.to_csv(os.path.join(dir_data, 'final', 'air_raid_alerts_times.csv'), index = False)

    # rename Kiev to Kiev City and Kiev Oblast to Kiev
    airlong['region_english'] = airlong['region_english'].str.replace("Kyiv oblast", "Kiev")
    airlong['region_english'] = airlong['region_english'].str.replace("Kyiv", "Kiev City")



    # assign ID_1 row value to airlong row if region_english matches any of the VARNAME_1 columns
    airlong = airlong[airlong['region_english'].notnull()]
    airlong['region_english'] = airlong['region_english'].str.replace(' oblast', '').str.replace(' raion', '')

    airlong.reset_index(inplace = True)


    # merge level 1 boundaries
    airlong = getAndMergeShapefile('ukr_admin_level1', airlong, 'NAME_1')

    # merge level 2 boundaries
    airlong = getAndMergeShapefile('ukr_admin_level2', airlong, 'name_2', names_col = "varname_2")



    # fix some unresolved matches
    fix = {'Bakhmut' :  "Bakhmats'kyi", # same as below
    # 'Bakhmut hromada', -- district abolished
    'Bila Tserkva' : "Bilia?vs'kyi", # -- city but shapefile has district
    'Dobropillia hromada' : "Dobropil's'kyi", # -- not sure
    "Izium" : "Iziums'ka", # -- city, but not sure which hromada it belongs to (Iziums'ka or Iziums'kyi? and whether these cover more than the city)
    'Ivano-Frankivsk' : "Ivano-Frankivs'k",
    # 'Kropyvnytskyi',
    'Kryvyi Rih' : "KryvyiRig",
    'Lubny' : "Lubens'ka", # -- city, checked the name difference, this is correct
    'Nikopol' : "Nikopol's'ka", # city
    'Pershotravensk': "Pershotravnevyi", # same as below
        'Pokrovsk' : "Pokrovs'kyi", # same as below
    'Polohiv' : "Polohivs'kyi", # same as below
    "Starokostiantyniv" : "Starokostiantynivs'kyi", # this is a city, but the shapefile has a district
    # 'Toretsk hromada', -- Torezka hromada?
    'Uman hromada' : "Umans'kyi",
    'Vasylkiv' : "Vasylivs'kyi",
    'Zaporizhzhia' : "Zaporiz'kyi"
    }

    # match airlong region_english with fix dictionary
    airlong['name_3'] = airlong['region_english'].map(fix)

    # if name_1 missing, assign name_2 to name_1, if name_2 and name_1 missing, assign name_3 to name_1
    airlong['NAME_1'] = airlong['NAME_1'].fillna(airlong['name_2'])
    airlong['NAME_1'] = airlong['NAME_1'].fillna(airlong['name_3'])

    print("Missing {} rows out of {}".format(airlong[airlong['NAME_1'].isnull()].shape[0], airlong.shape[0]))

    airlong.drop(columns = ['name_2', 'name_3'], inplace = True)

    gdf1 = readGdf('ukr_admin_level1')
    gdf = readGdf('ukr_admin_level2').to_crs('epsg:4326').append(gdf1.to_crs('epsg:4326'))
    gdf['NAME_1'] = gdf['NAME_1'].fillna(gdf['name_2'])

    airlong = airlong.merge(gdf[['NAME_1', 'geometry']], on = 'NAME_1', how = 'left')
    
    # check alert duration for matched regions
    airlong['duration'] = (airlong['end_dt'] - airlong['start_dt']).dt.total_seconds() / 60
    
    
    matched = pd.read_csv(os.path.join(dir_data, 'final', 'estimates_filtered', 'coeffs_distance_sum_.csv'))
    matched['region_english'] = matched['id'].str.slice(19)
    matched['start_date'] = matched['id'].str.slice(0, 19)
    matched.start_date = pd.to_datetime(matched.start_date)
    # round start_date down to minute 
    airlong['start_date'] = airlong.start_dt.dt.floor('min')
    # left antijoin airlong with matched
    check = airlong.merge(matched[['region_english', 'start_date']], on = ['region_english', 'start_date'], how = 'inner', indicator = True)
    check.duration.describe()
    
    # 90th percentile of duration
    check.duration.quantile(.99)

    airlong.drop(columns = ['geometry', 'start_date']).to_csv(os.path.join(dir_data, 'final', 'air_raid_alerts.csv'), index = False)

    # retain only unique region_english and geometry rows
    temp = airlong[['region_english', 'geometry']].drop_duplicates(subset = ['region_english', 'geometry'])
    temp = gpd.GeoDataFrame(temp, geometry = temp.geometry, crs = {'init': 'epsg:4326'})


    # plot raid heat map
    plotHeatMap(airlong)

    # inspect region-name match on map
    m = temp[['geometry', 'region_english']][temp.geometry != None].explore(color = "orange")

    path = os.path.join(dir_data, 'raw', 'geo', "ukr_admin_level1")
    gdf = gpd.read_file(path)

    # overlay gdf.boundary onto folium map m
    def style(feature):
        return {
            'color': 'dimgray'
        }
    folium.GeoJson(gdf.boundary, name = 'geojson', style_function = style).add_to(m)
    folium.GeoJson(temp.boundary, name = 'geojson', style_function = style).add_to(m)

    img_data = m._to_png(1)
    img = Image.open(io.BytesIO(img_data))
    img = img.resize((int(img.size[0]), int(img.size[1])), Image.ANTIALIAS)
    # font = ImageFont.truetype(r'/usr/share/fonts/truetype/lyx/cmex10.ttf', 50) 
    # draw = ImageDraw.Draw(img)
    # draw.text((30, 30), mon, font = font, align ="left") 
    img.save(os.path.join(dir_main, "output", "figs", "raid_regions.png"))



    temp.to_file(os.path.join(dir_data, 'final', 'air_geos/air_geos.shp'))
        
        



        

    

