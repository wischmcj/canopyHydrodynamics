
import datetime
from datetime import time 
import matplotlib.cm as cm

import os
import sys
import glob
from itertools import product

# dependencies
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#adding the path to the sys.path to allow for 
# easier relative referencing
sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())

from src.canhydro.global_vars import log

#would be better if we read them all in
# seperately
dataframes = pd.read_excel(
    "../data/validation_data/all_validation_data.xlsx",
    sheet_name=None,
)


#Dates were added to the data log in diff formats 
# i.e. mm-dd-yy vs mm-dd-yy
#   this dict mapps the incorrect dates to the correct dates 
date_fix_dict = {
    datetime.datetime(2020, 1, 8): datetime.datetime(2020, 8, 1),
    datetime.datetime(2020, 3, 8): datetime.datetime(2020, 8, 3),
    datetime.datetime(2020, 4, 8): datetime.datetime(2020, 8, 4),
    datetime.datetime(2020, 5, 8): datetime.datetime(2020, 8, 5),
    datetime.datetime(2020, 6, 8): datetime.datetime(2020, 9, 6),
    datetime.datetime(2020, 7, 9): datetime.datetime(2020, 9, 7),
    datetime.datetime(2022, 1, 8): datetime.datetime(2020, 8, 1),
    datetime.datetime(2022, 3, 8): datetime.datetime(2020, 8, 3),
    datetime.datetime(2022, 4, 8): datetime.datetime(2020, 8, 4),
    datetime.datetime(2022, 5, 8): datetime.datetime(2020, 8, 5),
    datetime.datetime(2022, 6, 8): datetime.datetime(2020, 9, 6),
    datetime.datetime(2022, 7, 9): datetime.datetime(2020, 9, 7),
    datetime.datetime(2020, 1, 10):datetime.datetime(2020, 10, 1),
    datetime.datetime(2020, 4, 10):datetime.datetime(2020, 10, 4),
    datetime.datetime(2022, 2, 8): datetime.datetime(2022, 8, 2),
    datetime.datetime(2022, 11, 8):datetime.datetime(2022, 8, 11),
    datetime.datetime(2022, 9, 10):datetime.datetime(2022, 9, 10),
    datetime.datetime(2022, 2,11): datetime.datetime(2022, 11, 2),
    datetime.datetime(2022, 1, 10):datetime.datetime(2020, 10, 1),
    datetime.datetime(2022, 4, 10):datetime.datetime(2020, 10, 4),
    datetime.datetime(2020, 2, 8): datetime.datetime(2022, 8, 2),
    datetime.datetime(2020, 11, 8):datetime.datetime(2022, 8, 11),
    datetime.datetime(2020, 9, 10):datetime.datetime(2022, 9, 10),
    datetime.datetime(2020, 2,11): datetime.datetime(2022, 11, 2),
}

#Hydromet data columns
sub_categories = ['N','SE','SW','C'] # to be ignored
rain_columns = ['Rain','10-08','16-14','18-13','02-26','24-7','26-3','27-5','04-12','32-1','32-6','32-07','T1'] #hold total rain amt

def get_event(rain_measurement: pd.DataFrame, 
              start_end_tuples: pd.Series = None):
    """
        Determines if a rain measurement was taken 
            during an already labeled event
    """
    if not start_end_tuples:
        start_end_tuples = get_rain_events()
    event_id =None
    for event, dates in start_end_tuples.iterrows():
        print(f'Event: {event} Start: {dates["start_datetime"]} End: {dates["end_datetime"]}')
        if rain_measurement >= dates['start_datetime'] and rain_measurement <= dates['end_datetime']:
            if event_id:
                print(f'Error: {rain_measurement["TIMESTAMP"]} is in multiple events')
                raise Exception
            else:
                event_id = event
    return event_id 

def get_rain_events(rain_events):
    start_dates = rain_events[['start_datetime','Event']]
    end_dates = rain_events[['end_datetime', 'Event']]
    start_end_tuples = start_dates.set_index('Event').join(end_dates.set_index('Event'), on='Event',rsuffix='_')
    return start_end_tuples

def get_hydromet_data(get_existing = True):
    if get_existing:
        try: 
            hydro_met_data_agg = pd.read_csv('../data/output/aggregated/hydro_met_data_agg.csv')
            return hydro_met_data_agg
        except FileNotFoundError as e:
            log.info('Could not find existing hydromet data, creating new data...')
    else:
        #measurements from a hydrometer, every few seconds durring storm
        hydro_met_data_with_seperators = dataframes['hydromet']
        
        #rain events, grouped based on periods of consistent
        #  rain, as measured by hydrometer
        rain_events = dataframes['hg_event_summaries']

        # filtering out subcategories that break down the total rain
        # measurement by cardinal direciton
        needed_columns = [x for x in hydro_met_data_with_seperators.columns if not any([sub in x for sub in sub_categories])]    
        
        # Finding start/end time, duration
        hydro_met_data_with_seperators['is_date'] = hydro_met_data_with_seperators['Date'].map(lambda x: isinstance(x, datetime.datetime))
        hydro_met_data = hydro_met_data_with_seperators[hydro_met_data_with_seperators['is_date']][needed_columns]
        
        #getting mins and maxs by date (needs fixing as 
        #   some storms) occured over night 
        date_groupings = hydro_met_data.groupby('Date')
        date_group_indexes = date_groupings['TIMESTAMP'].transform('min')
        hydro_met_data = hydro_met_data.join(date_group_indexes, rsuffix = '_min')
        hydro_met_data = hydro_met_data.join(date_group_indexes, rsuffix = '_max')

        #finding how soon after the start of a day's rain
        # a measurement was taken
        time_start_to_curr= hydro_met_data_agg['TIMESTAMP'] - hydro_met_data_agg['TIMESTAMP_min']
        hydro_met_data_agg['sec_since_start']=time_start_to_curr.apply(lambda x: x.total_seconds())

        # getting a list of pre-labeled rain events 
        event_start_end_dtms = get_rain_events()

        # Assign event id if measurement was in an event 
        # for those not in event, go off date for now
        #   new events to be created to account for these
        hydro_met_data_agg['event_id'] = hydro_met_data_agg['TIMESTAMP'].apply(lambda x: get_event(x, event_start_end_dtms))
        hydro_met_data_agg['group'] = hydro_met_data_agg.apply(lambda x: x['event_id'] if x['event_id'] else x['Date'])
        

        # Sum of rain for each event/ non-event date
        event_groupings = hydro_met_data.groupby('group')
        sum_group_indexes = event_groupings[rain_columns].transform('sum')
        cumsum_group_indexes = event_groupings[rain_columns].cumsum()
        hydro_met_data = hydro_met_data.join(cumsum_group_indexes, rsuffix = '_cum')
        hydro_met_data_agg= hydro_met_data.join(sum_group_indexes, rsuffix = '_sum')

        hydro_met_data_agg.to_csv('../data/output/aggregated/hydro_met_data_agg.csv')
    return hydro_met_data_agg


# using hydrometdata
def combine_rain_and_tree_data():
    #place holder    

    # rain_data_by_tree = hydro_met_data[needed_columns]
    pass

def consolidate_rain_data(rain_data:pd.DataFrame = None, 
                            volume:pd.DataFrame = None, 
                            depth:pd.DataFrame = None):
    #Here we start setting indicies and joining the data together
    if not rain_data:
        rain_data = get_hydromet_data() 
    if not volume:
        volume = dataframes['Volume']
    if not depth:    
        depth = dataframes['Depth']

    # Set indexes to easily access rows/join data sets

    # tree_mappings.set_index("Code", inplace=True)
    try:
        depth.set_index("Date", inplace=True)
    except KeyError as e:
        print(f'Key already set {e}')
    try:
        volume.set_index("Date", inplace=True)
    except KeyError as e:
        print(f'Key already set {e}')
    try:
        rain_data.set_index("Date", inplace=True)
    except KeyError as e:
        print(f'Key already set {e}')


    # date_fix_dict = {v:k for k,v in date_fix_dict.items()}

    volume.index = volume.index.map(lambda x: date_fix_dict.get(x, x))

    depth.index = depth.index.map(lambda x: date_fix_dict.get(x, x))


    volume_w_intensity = volume.join(rain_data, how='left',on='Date', rsuffix = '_hm')
    depth_w_intensity = depth.join(rain_data, how='left',on='Date', rsuffix = '_hm')
    volume_w_intensity['Duration'] = volume_w_intensity['TIMESTAMP_min'] - volume_w_intensity['EndTime']
    volume_w_intensity['Duration'] = volume_w_intensity['Duration'].map(lambda x: x.total_seconds())

    # tree_metrics = dataframes['Mapping'].loc[[tree_name]]
    rain_sorted_vol = volume_w_intensity.sort_values(by=["Rain"]).query("Rain>0")
    # rain_sorted_vol   = volume.sort_values(by=['Rain']).replace(0, np.nan)

    # rain_sorted_vol.replace(0, np.nan, inplace=True)

    rain_sorted_depth = depth_w_intensity.sort_values(by=["Rain"]).query("Rain>0")

    return rain_sorted_vol, rain_sorted_depth


def plot_rain_events(to_plot:pd.DataFrame =None, bin=10):
    if to_plot:
        hydro_met_data_agg = to_plot
    else:
        hydro_met_data_agg = get_hydromet_data()
    date_groupings = hydro_met_data_agg.groupby('event_id')
    event_ids = hydro_met_data_agg['event_id'].unique()

    # date_groupings = hydro_met_data_agg.groupby('event_id')
    # date_group_indexes = date_groupings[rain_columns].agg('sum')

    # date_group_indexes = date_groupings[rain_columns].agg('sum')
    # print(date_group_indexes)

    # rain_data_by_tree = hydro_met_data[needed_columns]
    # event_dates = hydro_met_data_agg['Date'].unique()
    
    # for idx, date in enumerate(event_dates):
    for event_id in event_ids:
        fig, ax = plt.subplots()
        # event = hydro_met_data_agg[hydro_met_data_agg['Date'] ==date].copy()
        event = hydro_met_data_agg[hydro_met_data_agg['event_id'] ==event_id].copy()
        dates = event['Date'].unique()
        # event  = event[~event['event_id'].isna()]
        # periods_of_significant_rain = retain_quantile(event, 'Rain_sum', .25)
        # event[event['Rain']>]
        colors = iter(cm.rainbow(np.linspace(0, 1, len(rain_columns)*2)))
        tree_cols = [col for col in rain_columns if not col=='Rain']
        ys = []
        for date in dates:
            print(f'running data for {date}')
            x_base = event[event['Date'] ==date]
            ax1 = ax.twinx()
            for column in tree_cols:
                print(event[column].max())
                # for idy, ver in enumerate([column, column +'_cum']):
                color = next(colors)
                # print(color)
                cut = (event[column]>0)
                y = x_base[column]
                x = x_base['sec_since_start']/60
                # print(column)
                # ax[idy].set_title(f'{ver} {date}')
                # ax[idy].plot(x,y,c=color, label= ver )
                ax.set_title(f'{event_id}-{date}')
                ax.plot(x,y,c=color, label= column ) 
            ax1.plot (x,event['Rain'], c = 'red')      
            # ax[0].legend()
            # ax[1].legend()
        plt.show()


def get_all_data(get_existing = True):
    # Try to get existing, if not there get new
    #read in all sheets of data from
    # combined validation data excel
    volume          = dataframes["Volume"]
    depth            = dataframes["Depth"]
    tree_mappings = dataframes["Mapping"]
    rain_event_summary = dataframes['hg_event_summaries']

    #Reading diff sheets into diff dataframes
    # data to be cleaned and joined
    get_hydromet_data(get_existing = False)

if __name__ == "__main__":
    get_all_data(get_existing = True)
    # get_event_ids(rain_events