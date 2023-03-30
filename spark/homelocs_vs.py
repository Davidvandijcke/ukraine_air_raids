"""
Likely Home Locations
"""

import logging
import sys

import pandas as pd
import numpy as np


from collections import namedtuple
from datetime import timedelta, date, datetime
from statistics import *
import re


from pyspark import SparkContext
from pyspark.sql import SQLContext, SparkSession
from pyspark.sql.types import (
    StructType,
    StructField,
    StringType,
    IntegerType,
    DoubleType,
    TimestampType,
)
from pyspark.sql.functions import (
    from_utc_timestamp,
    to_utc_timestamp,
    dayofyear,
    col,
    regexp_replace,
    unix_timestamp,
    monotonically_increasing_id,
    pandas_udf,
    PandasUDFType,
    col,
    asc,
    lit,
    countDistinct,
)


from sklearn.cluster import dbscan

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger("dwell_analysis.homelocations")

Config = namedtuple("Config", ["timezone"])

def floor_date(day):
    ''' Returns start of week (Sunday before) for compatibility with R's lubridate::floor_date
    @day: datetime object
    '''
    return day - timedelta(days=day.weekday() + 1)



def get_dwell_cluster(data: pd.DataFrame) -> pd.DataFrame:
    """Find a cluster if it exists and return the center and classification. For
    cases where there are more than one cluster, 
    
    all clusters are returned along with
    the counts of number of members (overnight dwells)

    Returns:
      pandas.DataFrame: results

    """
    caid = data.iloc[0].caid
    temp = data.copy()
    coords = np.radians(temp[["latitude", "longitude"]].values)
    # choose appropriate epsilon values
    # here we use ~50 meters
    kms_per_radian = 6371
    epsilon = 0.05 / kms_per_radian

    # calculate clusters
    # use haversine metric for calculating approximate distances on earth's
    #  surface (crow fly)
    _, cluster_labels = dbscan(
        coords, eps=epsilon, min_samples=3, algorithm="ball_tree", metric="haversine",
    )
    temp["cluster_labels"] = [f"l{l}" for l in cluster_labels]
    temp = temp[temp.cluster_labels != "l-1"]
    

    cluster_centers = (
        temp.groupby("cluster_labels")
        .agg({"latitude": "mean", "longitude": "mean", "caid": "count"})
        .reset_index()
        .rename(columns={"caid": "num_members"})
        .sort_values("num_members", ascending=False)
    )
    cluster_centers["caid"] = caid
    # select the top one (ordered by num_members)

#     if (cluster_centers.empty) or (len(cluster_centers.columns) != 3):
#         # change the schema as needed
#         cluster_centers = pd.DataFrame({'caid': pd.Series([], dtype='str'),
#                             'latitude': pd.Series([], dtype='float64'),
#                             'longitude': pd.Series([], dtype='float64')})

    return cluster_centers[["caid", "latitude", "longitude"]].iloc[0:1]




if __name__ == "__main__":

    # read commandline args
    #s3://ipsos-dvd/ukr/data/dwells_2022-01-01_2022-01-03 s3://ipsos-dvd/ukr/data/homes_2022-01-01_2022-01-03 2022-01-01 2022-01-03 GMT+3 1082
    bucket = sys.argv[1] # year for which to parse data
    dest = sys.argv[2]
    START_DATE = datetime.strptime(sys.argv[3], "%Y-%m-%d") # enter start date for considering home location in given format
    END_DATE = datetime.strptime(sys.argv[4], "%Y-%m-%d") 
    TIMEZONE = "EET"
    N_PARTITIONS = int(sys.argv[6])
    


    config = Config(TIMEZONE)
    # setup spark
    spark = SparkSession.builder.appName(
        f"run_home_locations"
    ).getOrCreate()


    cluster_schema = StructType(
        [
            StructField("caid", StringType(), False),
            StructField("latitude", DoubleType(), False),
            StructField("longitude", DoubleType(), False),
            StructField("utc_timestamp_start", DoubleType(), False),
            StructField("stay_length_sec", IntegerType(), False),
        ]
    )

    if config.timezone != "": # TODO
            overnight_dwells_full = (
        spark.read.parquet(bucket)
        .where("stay_length_sec > 7200")
        .withColumn( "date", 
        from_utc_timestamp(
                    col("utc_timestamp_start").cast(dataType=TimestampType()),
                    config.timezone,
                )
        )
        .where( col("date").between(START_DATE, END_DATE )) # only include dates in window
        .withColumn(
            "doy_start",
            dayofyear(
                col("date")
            ),
        )
        .withColumn(
            "doy_end",
            dayofyear(
                from_utc_timestamp(
                    (col("utc_timestamp_start") + col("stay_length_sec")).cast(
                        dataType=TimestampType()
                    ),
                    config.timezone,
                )
            ),
        )
        .where("doy_start != doy_end")
        .withColumn("rowid", monotonically_increasing_id())
        .withColumn("num_days", col("doy_end") - col("doy_start"))
        .select("caid", "latitude", "longitude", "num_days", "rowid", "date")
    )

    else:
        # Load data for period of interest, computing number of days spent overnight at each dwell location
        overnight_dwells_full = (
            spark.read.parquet(bucket)
            .where("stay_length_sec > 7200")
            .withColumn('TZID', regexp_replace(col("TZID"), "[\[\'\]]", "")) # split out frame lists
            .withColumn( "date", 
            from_utc_timestamp(
                        col("utc_timestamp_start").cast(dataType=TimestampType()),
                        col('TZID'),
                    )
            )
            .where( col("date").between(START_DATE, END_DATE )) # only include dates in window
            .withColumn(
                "doy_start",
                dayofyear(
                    col("date")
                ),
            )
            .withColumn(
                "doy_end",
                dayofyear(
                    from_utc_timestamp(
                        (col("utc_timestamp_start") + col("stay_length_sec")).cast(
                            dataType=TimestampType()
                        ),
                        col('TZID'),
                    )
                ),
            )
            .where("doy_start != doy_end")
            .withColumn("rowid", monotonically_increasing_id())
            .withColumn("num_days", col("doy_end") - col("doy_start"))
            .select("caid", "latitude", "longitude", "num_days", "rowid", "date")
        )

    ## loop over 5-week windows in selected period
    numdays = (END_DATE - START_DATE).days
    date_list = [START_DATE + timedelta(days = x) for x in range(1,numdays,30)]


    for date in date_list: # should make this a pandas function

        overnight_dwells = overnight_dwells_full.where(col("date").between(date - timedelta(days = 21), 
        date + timedelta(days = 21))) # get 5-week window around date

        # explode dwells into multiple overnight stays
        out_schema = StructType(
            [
                StructField("caid", StringType(), False),
                StructField("latitude", DoubleType(), False),
                StructField("longitude", DoubleType(), False),
            ]
        )

        def dwell_exploder(row):
            """Wrapper for dstar_clusters"""
            temp = row.iloc[0]
            df = pd.DataFrame(
                [
                    {
                        "caid": temp.caid,
                        "latitude": temp.latitude,
                        "longitude": temp.longitude,
                    }
                    for idx in range(max(1,temp.num_days)) # this is just defensive programming to fix a bug, num_days should always be >= 1
                ] 
            )


            return df

        exploded = (
            overnight_dwells.groupby("rowid")
            .applyInPandas(dwell_exploder, schema = out_schema)
            .select("caid", "latitude", "longitude")
            .cache()
        )

        # caids with sufficient night dwells
        keepers = (
            exploded.groupby("caid")
            .agg({"latitude": "count"})
            .where(col("count(latitude)") >= lit(4))
            .select("caid")
        )
        final = (
            exploded.join(keepers, exploded.caid == keepers.caid, how="inner")
            .drop(keepers.caid)
            .repartition(N_PARTITIONS, "caid")
        )

        # Cluster exploded overnight dwells 
        cluster_schema = StructType(
            [
                StructField("caid", StringType(), False),
                StructField("latitude", DoubleType(), False),
                StructField("longitude", DoubleType(), False),
            ]
        )
        def cluster_getter(data):
            """Wrapper for get_dwell_cluster"""
            return get_dwell_cluster(data)

        home_locations = final.groupby("caid").applyInPandas(cluster_getter, schema = cluster_schema)

        home_locations = home_locations.withColumn('start_date', lit(datetime.strftime(date, "%Y-%m-%d")))
        # home_locations = home_locations.toPandas()

        if date == date_list[0]:
            home_df = home_locations
        else:
            home_df = home_df.union(home_locations) # pd.concat([home_df, home_locations], ignore_index = True,axis = 0)


    home_df.write.mode("overwrite").parquet(dest)