# build in libs
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

#mostly for documentation but these
# are use for dynamically created 
# variable names 
abbrvs = { 'projected_area':'pa',
          'surface_area':'sa',
          'volume':'vol'}



# evently spaced angles to run for 
#  each tree in the sensitivity analysis
angles_to_run = [ 2.04,1.96,1.86, 1.8,1.72,1.66, 1.64, 1.58, 1.5, 1.42, 1.34,
            1.26, 1.18, 1.1, 1.02, 0.96, 0.88,
            0.8, 0.72, 0.64, 0.56, 0.48, 0.4, 0.32, 0.24, 0.16,
            0.08, 0, -0.08, -0.16, -0.24, -0.32, -0.4, -0.48,
            -0.56, -0.64, -0.72, -0.8, -0.88, -0.96, -1.02,
            -1.1, -1.18, -1.26, -1.34, -1.42, -1.5, -1.58
            ,-1.66,-1.64,-1.72,-1.8,-1.88,-1.96,-2.04]

# Columns output by CylinderCollection's
#   statistics and flow statistics funcs.
flow_cols = ["num_cylinders","projected_area","surface_area","angle_sum","volume","sa_to_vol","drip_node_id","drip_node_loc",'timestamp']
stats_cols = ["total_psa","psa_w_overlap","stem_psa","stem_psa_w_overlap","tot_surface_area","stem_surface_area","tot_hull_area","tot_hull_boundary","stem_hull_area","stem_hull_boundary","num_drip_points","max_bo","topQuarterTotPsa","topHalfTotPsa","topThreeQuarterTotPsa","TotalShade","top_quarter_shade","top_half_shade","top_three_quarter_shade","DBH","volume","X_max","Y_max","Z_max","X_min","Y_min","Z_min","Order_zero_angle_avg","Order_zero_angle_std","Order_one_angle_avg","Order_one_angle_std","Order_two_angle_avg","Order_two_angle_std","Order_three_angle_avg","Order_three_angle_std","order_gr_four_angle_avg","order_gr_four_angle_std","file_name","timestamp"]

#ratio metrics to calculate
#   (a,b,c) => a = b/c
volume_ratios = [('vol_over_tps','volume','total_psa'),
                        ('vol_over_sps','volume','stem_psa'),    
                        ('vol_over_tsa','volume','tot_surface_area'),
                        ('vol_over_ssa','volume','stem_surface_area'),
                        ('vol_over_tha','volume','tot_hull_area'),
                        ('vol_over_sha','volume','stem_hull_area'),
                        ('vol_over_dbh','volume','DBH'),
                        ('vol_over_oza','volume','Order_zero_angle_avg'),
                        ('vol_over_ota','volume','Order_two_angle_avg')]
stem_ratios = [('stem_vol_over_tps','volume','total_psa')]

# Flow Stats Config

# Calculate for drip nodes (dn) which
#   are defined as the entire drip shed (0th %ile)
# as well as drip points (pd), the 98th %ile
scopes = [('dn',0),('dp',.98)]

# List of metrics to calculate for the drip shed
deciding_drip_metrics       = ['projected_area',
                                'surface_area',
                                'volume', 
                                'sa_to_vol']
drip_point_summary_funcs = [  np.max,
                                np.sum,
                                np.mean]
stem_metrics        = ['projected_area',
                                'surface_area',
                                'volume', 
                                'sa_to_vol',
                                'num_cylinders']
stem_summary_metrics = [np.sum]



# Specific to paper 
paper_cases = [0.04,-0.14, -0.34]
paper_trees = 'Secrest27-05_000000','Secrest32-06_000000'
# paper_trees[['dp_count', 'dp_max_pa', 'dp_pa', 'dp_avg_pa', 'dp_sa', 'dp_vol',
#        'dp_sum_sa_to_vol', 'dp_avg_sa_to_vol', 'dp_agg_sa_to_vol', 'dn_count',
#        'dn_max_pa', 'dn_pa', 'dn_avg_pa', 'dn_avg_sa']]



def get_files(tree_codes:list[str]=None, 
              dir = 'output', 
              types = ['statistics','flows']):
    files = []
    for f_type in ['statistics','flows']:
        if f_type in types:
            files.append(glob.glob(f'./data/{dir}/{f_type}/*.csv'))
    if tree_codes:
        breakpoint()
        files=[f for f in files if any([code in f for code in tree_codes])]
    return files


#Throws alot of warnings since the last column, time stamp gets dropped off since it doesnt have a col name
def clean_sa_file_names(file_list,remove_zeros, flows):
    file_names = [(f,os.path.basename(f).replace('_1_','').replace(
            'Secrest07-32_000000_Secrest07-32_000000','Secrest07-32_000000').replace(
            'Secrest08-24c_000000_Secrest08-24c_000000','Secrest08-24c_000000').replace(
            'Secrest10-02_000000_Secrest10-02_000000','Secrest10-02_000000')) for f in file_list]
    
    to_remove_from_name = ['_flows','_statistics','.csv','_flows']
    for rem_string in to_remove_from_name:
        file_names = [(x,y.replace(rem_string,'')) for x,y in file_names]

    if remove_zeros: 
        if flows:
            name_parts = [(name.split('_000000_')[0], name.split('_000000_')[1], ) for f, name in file_names if '_000000_' in name]
        else:
            name_parts = [(name.split('_000000_')[0], name.split('_000000_')[1], pd.read_csv(f,index_col=False,names=stats_cols,skiprows=1  )) for f, name in file_names if '_000000_' in name]
    else:
        if flows:
            name_parts = [(name.split('_000000_')[0] + '_000000', name.split('_000000_')[1], pd.read_csv(f,index_col=False,names=flow_cols,skiprows=1)) for f, name in file_names if '_000000_' in name]
        else:
            name_parts = [(name.split('_000000_')[0] + '_000000', name.split('_000000_')[1], pd.read_csv(f,index_col=False,names=stats_cols,skiprows=1)) for f, name in file_names if '_000000_' in name]

    # print(name_parts)
    return name_parts


# Read each

def read_and_get_cases(tree_codes:list[str]=None):
    # read by default 1st sheet of an excel file

    # Get a list of all the csv files in the target directory
    stats_files = get_files(types = ['statistics'])
    flow_files = get_files(types = ['flows'])
    
    all_tree_stats = pd.DataFrame()
    dfs_and_cases = clean_sa_file_names(stats_files, remove_zeros = False, flows = False)

    file_names = []
    dfs=[]
    for tup in dfs_and_cases:
        file_name, case_name, df = tup
        df['file_name'] = file_name
        df['case_name'] = case_name
        df['join_field']=f'{file_name},{str(case_name)}'
        file_names.append(file_name)
        dfs.append(df)

    all_tree_stats = pd.concat(dfs)

    # print(all_tree_stats.size)
    # if debug: if debug: if debug: if debug: print(len(all_tree_stats['file_name']))
    # print(len(dfs_and_cases))

    flow_dfs_and_cases = clean_sa_file_names(flow_files, remove_zeros = False, flows = True)
    flow_dfs = []

    for idx, tup in enumerate(flow_dfs_and_cases):
        # print(tup)
        file_name, case_name, df = tup
        df['file_name'] = file_name
        df['case_name'] = case_name
        df['join_field'] = f'{file_name},{str(case_name)}'
        is_stem = [cyl_id==0 for cyl_id in df['drip_node_id']]
        df['is_stem'] = is_stem
        flow_dfs.append(df)


    all_tree_flows = pd.concat(flow_dfs)

    return all_tree_stats, all_tree_flows
    # test_read= pd.read_csv('../data/output/aggregated/flow_stats.csv')


def retain_quantile(df, field, percentile):
    """
        Used to identify drip points (e.g. the top X percentile
            of drip nodes)
    """
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val]

def return_quantile(df, field, percentile):
    quantile = retain_quantile(df, field, percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return quantile[field]

def add_ratio_fields(df:pd.DataFrame, 
                        ratio_metrics_list:list[tuple],
                        force = False):
    for metric in ratio_metrics_list:
        name, numerator, denominator = metric
        if name in df.columns:
            log.warning(f'''Column {name} already exists in dataframe,
                                if still desired, override with force = True''')
        df[name] = df[numerator].div(df[denominator])


#extract the tree codes from the file names
def summarize_flow_chars(all_tree_flows:pd.DataFrame, all_tree_stats:pd.DataFrame):
    stem_flows = all_tree_flows[all_tree_flows['is_stem']]
    drip_flows = all_tree_flows[all_tree_flows['is_stem']==False]
    drip_flows = drip_flows[drip_flows['projected_area']>0.001]

    drip_grouped = drip_flows.groupby(['file_name','case_name'])
    drip_stats = {}

    # all_metrics = product(scopes, drip_point_summary_funcs, deciding_drip_metrics)
    # get_name = lambda a,b,c: f'{a[0]}_{abbrvs.get(c,c)}_{b.__name__}'

    # # get_scope = lambda (a,b,c): return_quantile(x[1], metric, pct)

    # # i.e. [(('dn',0,),'projected_area', max),(('dp',.98),'surface_area', sum), ...]
    # criteria_w_names = [(get_name(tup), tup[1], tup[2], tup[0][1], ) for tup in all_metrics]
    # # [('dn_pa_max', max, 'projected_area', 0), ('dp_sa_sum', sum, 'surface_area', .98), ... ]

    # get_dp_grps = lambda x: {name:func(return_quantile(x[1], metric, pct)) for name, func, metric, pct in criteria_w_names}
    # drip_point_groups = map(get_dp_grps,drip_grouped)
    # # drip_grouped = [((<file_name>,<case_name>), <df_of_fields>), ...]
    # # dpg = [((<file_name>,<case_name>)):{'dn_pa_max': max(return_quantile(df, 'projected_area', 0)), ...}, ...] 

    # calculating stats for different definitions of drip points
    drip_stats = []
    for group in drip_grouped:
        drip_stat = {}
        for scope in scopes:
            for metric in deciding_drip_metrics:
                drip_points = return_quantile(group[1], metric, scope[1])
                for func in drip_point_summary_funcs:
                    metric_name = f'{scope[0]}_{abbrvs.get(metric, metric)}_{func.__name__}'
                    value = func(drip_points)
                    drip_stats[metric_name] = value
        drip_stat['join_field'] = f'{group[0][0]},{str(group[0][1])}'
        drip_stats.append(drip_stat)
    
    #Calculating stem flow stats, for which there 
    #  is only one value per (file,case) pair
    stem_stats = []
    stem_grouped = stem_flows.groupby(['file_name','case_name'])
    for group in stem_grouped:
        stem_stat = {}
        for metric in deciding_drip_metrics:
            for func in drip_point_summary_funcs:
                metric_name = f'stem_{abbrvs.get(metric, metric)}_{func.__name__}'
                value = func(group[1][metric])
                drip_stats[metric_name] = value
        stem_stat['join_field'] = f'{group[0][0]},{str(group[0][1])}'
        stem_stats.append(stem_stat)

    # preparing to join the two dataframes
        
        # the old way of doing the above (with less metrics )
        # drip_stats = [{'join_field': f'{group[0][0]},{str(group[0][1])}',
        #             'drip_avg_num_cyls': np.mean(group[1]['num_cylinders']), 
        #             'drip_max_num_cyls': np.max(group[1]['num_cylinders']), 
        #             'dp_prct': group[1]['projected_area'].quantile(0.98),
        #             'dp_count' :     len((top_by_pa)),
        #             'dp_max_pa' : np.max (top_by_pa),
        #             'dp_pa'     : np.sum (top_by_pa),
        #             'dp_avg_pa' : np.mean(top_by_pa),
        #             'dp_sa'     : np.sum(top_by_sa),
        #             'dp_vol'     : np.sum(top_by_volume),
        #             'dp_sum_sa_to_vol'     : np.sum(return_quantile(group[1],'sa_to_vol',0.98)),
        #             'dp_avg_sa_to_vol'     : np.mean(return_quantile(group[1],'sa_to_vol',0.98)),
        #             'dp_agg_sa_to_vol': np.sum(top_by_sa)/np.sum(top_by_volume),

        #             'dn_count' : len(group[1]['projected_area']),
        #             'dn_max_pa': np.max(group[1]['projected_area']),
        #             'dn_pa': np.sum(group[1]['projected_area']),
        #             'dn_avg_pa': np.mean(group[1]['projected_area']),
        #             'dn_avg_sa'     : np.sum(top_by_sa)
        #             }   for group in drip_grouped] 

    # stem_grouped = stem_flows.groupby(['file_name','case_name'])
    # stem_stats = [
    #                 { 
    #                     'join_field':  f'{group[0][0]},{str(group[0][1])}',
    #                     'stem_avg_num_cyls' :np.mean(group[1]['num_cylinders']), 
    #                     'stem_pa' :np.sum(group[1]['projected_area']), 
    #                     'stem_vol' :np.sum(group[1]['volume']),
    #                     'stem_sa' :np.sum(group[1]['surface_area']),
    #                     'stem_avg_sa_vol' :np.mean(group[1]['sa_to_vol'])
    #                 } for group in stem_grouped]



    drip_df = pd.DataFrame(drip_stats)
    stem_df = pd.DataFrame(stem_stats)
    stem_df['stem_agg_sa_to_vol'] = stem_df['stem_sa']/stem_df['stem_vol']
    return stem_df, drip_df

def get_tree_code(df_row):
    code = df_row['file_name'].replace('Secrest','').split('_')[0]
    return code
    
def prep_sensitivity_data(tree_codes:list[str]):
    # combine all flow data and stats data into respective dataframes
    all_tree_stats, all_tree_flows = read_and_get_cases(tree_codes)

    # identify drop point
    # s and aggregate flow data
    stem_chars, drip_chars = summarize_flow_chars()

    #extract tree codes (mostly for plot labels)
    all_tree_stats['tree_codes'] = all_tree_stats.apply(get_tree_code,  axis = 1)
    
    #sorting and joining
    all_tree_stats = all_tree_stats.dropna(axis=1, how='all')
    all_tree_stats = all_tree_stats.sort_values(by=['file_name', 'case_name'])
    all_tree_flows = all_tree_flows.sort_values(by=['file_name', 'case_name'])

    # w_flow_stats = all_tree_stats.join(drip_df, how='left', on= 'join_field', lsuffix='stats',rsuffix='drip')
    w_drip_stats = all_tree_stats.set_index('join_field').join(drip_chars.set_index('join_field'))
    w_flow_stats = w_drip_stats.join(stem_chars.set_index('join_field'))
    

    # Case names read from file names are in strings 
    # this portion converts them to decimals
    all_tree_stats['case_name'] = pd.to_numeric(all_tree_stats['case_name'], errors='coerce')
    w_flow_stats['case_name'] = pd.to_numeric(w_flow_stats['case_name'], errors='coerce')

    # Sort the dataframe by 'case_name' to 
    # make analytics easier 
    all_tree_stats = all_tree_stats.sort_values(by='case_name')
    w_flow_stats = w_flow_stats.sort_values(by='case_name')

    w_flow_stats = w_flow_stats.dropna(axis=1, how='all')

    # calculating remaining statistics;
    # Those that need both flow and stat data
    w_flow_stats['stem_over_dn_pa'] = w_flow_stats['stem_psa']/w_flow_stats['dn_pa']
    w_flow_stats['stem_over_dp_pa'] = w_flow_stats['stem_psa']/w_flow_stats['dp_pa']

    w_flow_stats['stem_over_else'] = w_flow_stats['stem_psa']/(w_flow_stats['total_psa'] - w_flow_stats['stem_psa'])


    w_flow_stats = w_flow_stats.dropna(axis=1, how='all')


    tree_traits = w_flow_stats
    try:
        tree_traits.set_index("tree_codes", inplace=True)
    except KeyError as e:
        print(f'Key already set {e}')
    return w_flow_stats


def evaluate_cases():
    """
        Determines what cases have been run and 
        which cases still need to be run 
    """
#getting cases that we have already run a sensitivity analysis on
    arr = [(x,1.3) if y =='-1.3.' else (x,y) for x,y in all_tree_stats[['file_name', 'case_name']].values ]

    run_cases = set([ (x,float(y)) for x,y in arr ])

    files = set([case[0] for case in run_cases])

    all_runs = [case[1] for case in run_cases]
    run_angles = list(set(all_runs))


    all_cases = product(files, angles_to_run)
    needed_cases = [case for case in all_cases if case not in run_cases]

    needed_files = set([case[0] for case in needed_cases])

    needed_all_runs = [case[1] for case in needed_cases]
    needed_run_angles = list(set(all_runs))

    return needed_files, needed_run_angles

def plot(run_angles, all_runs):
    # plot num run by angle 
    colors = ['' for _ in run_angles]
    cnts = [0 for _ in run_angles]
    for idx, angle in enumerate(run_angles):
        if angle in angles_to_run:
            colors[idx] = 'b'
        else:
            colors[idx] = 'r'
        cnts[idx] = all_runs.count(angle)
    print(run_angles)
    print(cnts)
    _, ax = plt.subplots()

    ax.scatter(run_angles, cnts, color=colors)
    plt.show()

if __name__ == "__main__":
    prep_sensitivity_data()
    evaluate_cases()