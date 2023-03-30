"""Stay location and duration finding using the D-Star algorithm


aws emr add-steps --cluster-id <Your EMR cluster id> --steps Type=spark,Name=TestJob,Args=[--deploy-mode,cluster,--master,yarn,--conf,spark.yarn.submit.waitAppCompletion=true,s3a://your-source-bucket/code/pythonjob.py,s3a://your-source-bucket/data/data.csv,s3a://your-destination-bucket/test-output/],ActionOnFailure=CONTINUE
"""

from collections import namedtuple
import logging
import sys

from geopy.distance import great_circle
import pandas as pd

from datetime import timedelta, date, datetime
from statistics import *

from pyspark import SparkContext
from pyspark.sql import SQLContext, SparkSession
from pyspark.sql.types import (
    StructType,
    LongType,
    StructField,
    StringType,
    DoubleType,
    IntegerType,
    TimestampType,
)
from pyspark.sql.functions import (
    from_utc_timestamp,
    to_utc_timestamp,
    dayofyear,
    col,
    unix_timestamp,
    monotonically_increasing_id,
    pandas_udf,
    PandasUDFType,
    col,
    asc,
    lit,
    countDistinct,
)

# Range of IDs to look at when considering a cluster candidate
WindowRange = namedtuple("WindowRange", ["min", "max"])

# Time range of candidate cluster
TimeRange = namedtuple("Range", ["start", "end", "diff"])



def get_window_range(current_index, window_size, index_min, index_max) -> WindowRange:
    """Give the trajectory window of points to include
    For cases at the beginning or end, the values are clamped and the upper or
    lower range is restricted.

    Normal use case:
    |   ...x... |

    Lower bound adjustment:
    |..x...     |

    Upper bound adjustment:
    |      ...x.|
    """
    min_range = max(current_index - window_size, index_min)
    max_range = min(current_index + window_size, index_max)
    # TODO: should this error if we get np.nan?

    return WindowRange(min=min_range, max=max_range)

def get_stay_period(
    data: pd.DataFrame, cluster_indexes: set, self_id: int
) -> TimeRange:
    """
    Get the length of time in a cluster.

    Assumptions:
    1. Data is already time ordered
    2. Index starts at zero and increases monotonically to n-1 where the dataset
       is of size (n, 4)
    3. If there are no neighbors, self should be start and end and diff == 0

    Args:
      data (pandas.DataFrame): trajectory data
      cluster_indexes (set of ints): Indexes of points that are part of stay period
        cluster
      self_id (int): Index of current focus point

    Returns:
      TimeRange: stay period of the indexes that correspond to a stay cluster
    """
    idx_min = self_id if len(cluster_indexes) == 0 else min(*cluster_indexes, self_id)
    idx_max = self_id if len(cluster_indexes) == 0 else max(*cluster_indexes, self_id)

    max_val = data.at[idx_max, "utc_timestamp"]
    min_val = data.at[idx_min, "utc_timestamp"]

    return TimeRange(start=min_val, end=max_val, diff=(max_val - min_val),)


def combine_stay_periods(period1: TimeRange, period2: TimeRange) -> TimeRange:
    """Given two stay periods, combine them to give an overall time range

    Args:
      period1 (TimeRange): duration to be combined
      period2 (TimeRange): duration to be combined

    Returns
      TimeRange: duration that encompasses period1 and period2
    """
    start_val = min(period1.start, period2.start)
    end_val = max(period1.end, period2.end)
    return TimeRange(start=start_val, end=end_val, diff=(end_val - start_val))


def duration_joinable_ids(target_start_time: int, stay_clusters: dict) -> set:
    """
    Enumerate all core points that have overlapping time ranges for clusters

    If a previous cluster ends after the target cluster starts, then they are
    duration joinable

    Args:
      target_start_time (int): Start time of current stay cluster
      stay_clusters (dict): Dict of cluster indexes, their neighbors, and their
        durations as TimeRanges

    Returns:
      set: Indexes (arbitrary) of other clusters that overlap with target cluster

    TODO: shouldn't only the most recent cluster be used for consideration here?
          One way to test is to see whether other clusters are marked as joinable
    """
    joinable_ids = set()

    for core_id in stay_clusters:
        if stay_clusters[core_id]["timerange"].end >= target_start_time:
            joinable_ids.add(core_id)
    return joinable_ids


def dstar(
    trajectory: pd.DataFrame,
    window_size: int = 10,
    dist_threshold: float = 25,
    min_time_threshold: float = 60,
) -> dict:
    """
    NOTE:
      all dataframes need to have an index that monotonically increases and starts
      at zero
    Args:
      trajectory is a pd.DataFrame with schema:
        * latitude (float64)
        * longitude (float64)
        * utc_timestamp (int64)
      window_size (int): number of time adjacent points to consider in
        calculating dwell location. Window size is added before and after the focus
        point. In cases where the focus point is closer to the beginning or end than
        the size of this window size, the lower/upper difference is the lower/upper
        window size instead.
      dist_threshold (float): distance in meters for points to be considered as
        part of a cluster
      min_time_threshold (float): minimum amount of time a cluster must endure
        to be considered a dwell event

    `trajectory` is assumed to be time-ordered ascending, and to represent a
      single user (by 'caid')
    """
    index_min, index_max = trajectory.index[0], trajectory.index[-1]
    # TODO: make this into a class with accessor methods
    # candidate_clusters and stay_clusters are both structured as:
    #  {'<point_id>': {
    #      'timerange': TimeRange(start=..., end=..., diff=...),
    #      'neighbors': set(...)
    #    },
    #   '<point_id2>': ... }
    #
    candidate_clusters = dict()
    stay_clusters = dict()

    for point in trajectory.itertuples():

        window_range = get_window_range(point.Index, window_size, index_min, index_max)
        current_window = trajectory.loc[window_range.min : window_range.max]

        # add neighbor entry if not exists
        if point.Index not in candidate_clusters:
            candidate_clusters[point.Index] = dict(neighbors=set())

        # compare neighbors to core point
        for point_j in current_window.itertuples():
            # skip self
            if point_j.Index == point.Index:
                continue
            # NOTE: this function is likely called redundantly. The results could be
            #       stored for lookup. if not in dict, calculate and add to dict,
            #       otherwise use. Is it possible to do a symmetric dict lookup?
            dist_apart = great_circle(
                (point_j.latitude, point_j.longitude),
                (point.latitude, point.longitude),
            ).meters

            if dist_apart <= dist_threshold:
                # If j is close to i, add it to j's and i's neighbor set
                candidate_clusters[point.Index]["neighbors"].add(point_j.Index)
                if point_j.Index not in candidate_clusters:
                    candidate_clusters[point_j.Index] = dict(neighbors=set())
                candidate_clusters[point_j.Index]["neighbors"].add(point.Index)

        # set stay period for current core point
        candidate_clusters[point.Index]["timerange"] = get_stay_period(
            trajectory, candidate_clusters[point.Index]["neighbors"], point.Index
        )

        # Avoid combining clusters if they don't have a lagged comparison point
        # This is the case when point index - window_size < 0
        if point.Index < window_size:
            continue

        # get stay duration of candidate cluster
        lagged_point_duration = get_stay_period(
            trajectory,
            candidate_clusters[window_range.min]["neighbors"],
            window_range.min,
        )

        # If stay is long enough, keep the cluster
        if lagged_point_duration.diff >= min_time_threshold:
            joinable_ids = duration_joinable_ids(
                lagged_point_duration.start, stay_clusters
            )
            if len(joinable_ids) > 0:
                # combine clusters since they overlap
                # NOTE: we only pop off one of the IDs, so we are ignore others if
                #       they exist
                if len(joinable_ids) > 2:
                    logging.warning(
                        "More than one joinable IDs were returned: "
                        f"{str(joinable_ids)} for `{window_range.min}`"
                    )
                current_cluster = joinable_ids.pop()
                # Update existing cluster with new
                stay_clusters[current_cluster]["neighbors"].update(
                    candidate_clusters[window_range.min]["neighbors"]
                )
                stay_clusters[current_cluster]["timerange"] = combine_stay_periods(
                    candidate_clusters[window_range.min]["timerange"],
                    stay_clusters[current_cluster]["timerange"],
                )
            else:
                # add to passed cluster set
                stay_clusters[window_range.min] = {
                    "timerange": candidate_clusters[window_range.min]["timerange"],
                    "neighbors": candidate_clusters[window_range.min]["neighbors"],
                }

    return stay_clusters


def get_trajectory_cluster(neighbors: dict, trajectory: pd.DataFrame) -> pd.DataFrame:
    """Given points making up a cluster, find the center lat/lng,
    stay duration, etc.

    Args:
        neighbors (dict): A dictionary of point IDs and the neighbors of the
          point
        trajectory (pandas.DataFrame): DataFrame of the entire trajectory
    """
    centers = []
    caid = trajectory.at[0, "caid"]
    for n_id in neighbors:
        temp = trajectory.loc[list(neighbors[n_id]["neighbors"])].agg(
            {"latitude": "mean", "longitude": "mean", "utc_timestamp": ["max", "min"],}
        )
        centers.append(
            {
                "caid": caid,
                "latitude": temp.loc["mean"].latitude,
                "longitude": temp.loc["mean"].longitude,
                "utc_timestamp_start": temp.loc["min"].utc_timestamp,
                "stay_length_sec": (
                    temp.loc["max"].utc_timestamp - temp.loc["min"].utc_timestamp
                ),
            }
        )
    out_cols = [
        "caid",
        "latitude",
        "longitude",
        "utc_timestamp_start",
        "stay_length_sec",
    ]
    if len(centers) > 0:
        return pd.DataFrame(centers)[out_cols]
    else:
        return pd.DataFrame(columns=out_cols)


def dstar_clusters(
    data: pd.DataFrame,
    window_size: int = 18,
    dist_threshold: float = 25,
    min_time_threshold: float = 60,
) -> pd.DataFrame:
    """Calculate clusters, prepare them as output dataframe"""
    return get_trajectory_cluster(
        dstar(data, window_size, dist_threshold, min_time_threshold), data
    )

def dstar_agg(data):
    """Wrapper for dstar_clusters"""
    return dstar_clusters(data)




def parse_dates(x):
    if "/" in x:
        start_date = datetime.strptime(x.split('/')[0], "%Y-%m-%d")
        end_date = datetime.strptime(x.split('/')[1], "%Y-%m-%d")
        delta = end_date - start_date
        date_list = []
        for i in range(delta.days + 1):
            date = start_date + timedelta(days = i)
            date_list.append(date.strftime("%Y/%m/%d"))
        return(date_list)
    else: 
        return([x.replace("-", "/")])
    


if __name__ == "__main__":

    ## set up spark



    # import os
    # import sys
    # spark_path = r"/home/antonvocalis/spark/spark-3.2.1-bin-hadoop3.2" # spark installed folder
    # os.environ['SPARK_HOME'] = spark_path
    # sys.path.insert(0, spark_path + "/bin")
    # sys.path.insert(0, spark_path + "/python/pyspark/")
    # sys.path.insert(0, spark_path + "/python/lib/pyspark.zip")
    # sys.path.insert(0, spark_path + "/python/lib/py4j-0.10.9.3-src.zip")
    # os.environ['PYSPARK_PYTHON'] = "python3"

    ## LOAD PINGS, FILTER ON COUNTRY AND DATES, AND WRITE TO FILE

    # get arguments
    dates = [sys.argv[1]]
    country = str(sys.argv[2])
    dest_pings = str(sys.argv[3])
    dest_dwells = str(sys.argv[4])
    N_PARTITIONS = int(sys.argv[6])

    # launch spark instance
    spark = SparkSession.builder.appName(f"parse_{country}_veraset").getOrCreate()

    # set up date ranges in path form
    datelist = []
    #args = ["2020-03-18/2020-03-20", "2020-04-01"]
    for arg in dates:
        temp = parse_dates(arg)
        datelist.extend(temp)

    # create list of s3 folders for dates to get
    base_bucket = sys.argv[5] # f"s3://external-veraset-data-us-west-2/us/" # f"s3://external-veraset-data-us-west-2/movement/"
    #base_bucket = f"/home/antonvocalis/ipsos/data/in/veraset/"
    pathlist = [base_bucket + date for date in datelist]


    north, east, south, west = 39.850721,-95.163574,47.376035,-74.234619

    # load data
    essential_fields = [
            StructField("utc_timestamp",LongType(),False),
            StructField("caid",StringType(),False),
            StructField("latitude",DoubleType(),False),
            StructField("longitude",DoubleType(),False),
            StructField("altitude",DoubleType(),False),
    ]
    raw_schema = StructType(
        essential_fields + [
            StructField("id_type",StringType(),False),
            StructField("geo_hash",StringType(),False),
            StructField("horizontal_accuracy",DoubleType(),False),
            StructField("ip_address",StringType(),False),
            #StructField("altitude",DoubleType(),False),
            StructField("iso_country_code",StringType(),False)]
    )
    df = (spark.read.schema(raw_schema).parquet(*pathlist)
            .select("caid", "latitude", "longitude", "utc_timestamp", "horizontal_accuracy", "altitude",)
            .filter(col("horizontal_accuracy") <= 150) # 150m or less horizontal accuracy (similar to safegraph)
            .drop(col('horizontal_accuracy'))
            #.filter(col('longitude').between(east,west) & col('latitude').between(north, south))
            )

    # filter jumpy pings: to do 
    #over(Window.partitionBy('caid').orderBy('utc_timestamp') 
    # calculate time difference
    # take gps pings within less than x seconds from each other and calculate the speed
    # drop the ones with extremely high speed


    df = df.filter( col('iso_country_code') == country) # filter on country
    df.write.mode("overwrite").parquet(dest_pings)
    df.unpersist()


    ## CALCULATE DWELLS AND WRITE TO FILE
    if dest_dwells != 'none':
        cluster_schema = StructType(
            [
                StructField("caid", StringType(), False),
                StructField("latitude", DoubleType(), False),
                StructField("longitude", DoubleType(), False),
                StructField("utc_timestamp_start", LongType(), False),
                StructField("stay_length_sec", IntegerType(), False),
            ]
        )
        essential_schema = StructType(essential_fields)
        clusters = (spark.read.schema(essential_schema).parquet(dest_pings)
                        .drop("altitude")
                        .repartition(N_PARTITIONS, "caid")
                        .sortWithinPartitions('utc_timestamp', ascending = True)
                        .groupby("caid").applyInPandas(dstar_agg, schema = cluster_schema))

        clusters.write.mode("overwrite").parquet(
            dest_dwells
        )
