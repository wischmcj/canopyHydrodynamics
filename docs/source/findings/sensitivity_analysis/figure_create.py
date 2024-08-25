from __future__ import annotations

import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from geopandas import GeoSeries
from sensitivity_analysis import run_test_cases

from src.canopyhydro.configuration import log


def retain_quantile(df, field, percentile):
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val]


def return_quantile(df, field, percentile):
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val][field]


def get_drip_point_coords(collection, percentile):
    """## Returns the x, y, z coordinates of drip points
     in the percentile indicated by projected area (found in'find_flow_components')

    ### Args:
        - `collection (_type_)`: _description_
        - `percentile (_type_)`: _description_

    ### Returns:
        - `_type_`: _description_
    """
    scale = 1
    flows = pd.DataFrame([flow.__dict__ for flow in collection.flows])
    flows = flows[flows["drip_node_id"] != 0]
    drip_points = retain_quantile(flows, "projected_area", percentile)
    drip_point_locs_x = [pt[0] * scale for pt in drip_points["drip_node_loc"]]
    drip_point_locs_y = [pt[1] * scale for pt in drip_points["drip_node_loc"]]
    drip_point_size = [pt * 5 for pt in drip_points["projected_area"]]
    return drip_point_locs_x, drip_point_locs_y, drip_point_size


if __name__ == "__main__":
    run_test_cases(
        [
            ("5_SmallTree", -0.16666),
            ("5_SmallTree", -0.16666),
            ("5_SmallTree", -0.16666),
        ],
        stats=True,
        fig=False,
        from_pickle=False,
    )

    # load_from_pickle('Secrest32-06_000000', 'stats', -0.36)
    # load_from_pickle('Secrest32-06_000000', 'stats', 0.04)
    # load_from_pickle('Secrest32-06_000000_1', 'stats', -0.14)

    db_files = [
        open(
            "/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.36",
            "rb",
        ),
        # open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_0.04', 'rb'),
        # open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.14', 'rb')
    ]

    file = "32-06_low_drip_points_36"
    # file = '32-06_high_drip_points_04'
    # file = '32-06_mid_drip_points_14'

    collections = [pickle.load(db_file) for db_file in db_files]

    drip_info = [plot_drip_points(col, 0.98) for col in collections]
    max_drip = 0
    for x, y, drip_point_size in drip_info:
        potential_max_drip = np.max(drip_point_size)
        if potential_max_drip > max_drip:
            max_drip = potential_max_drip

    for x, y, drip_point_size in drip_info:
        # ensuring each drip array has the same max drip_size
        # allows us to standarize the size
        x.append(-10)
        y.append(-10)
        # if we want to specify the max drip size from
        # a separate run of this code to ensure they are on
        # the same scale
        # drip_point_size.append(6.603)

    # polyss = [[cyl.projected_data['XY']['polygon'] for cyl in col.cylinders] for col in collection]

    # geohull = GeoSeries(just_right.stem_hull)
    # geohull = GeoSeries(high_end.stem_hull)

    # ax.scatter(x, y+8, marker = 'o')
    # fig, ax = plt.subplots()
    for idx, graph_data in enumerate(drip_info):
        if "32" in file:
            x_trans = 9
            y_trans = 4
            x_lim = 14
            y_lim = 13
        if "27" in file:
            x_trans = 2
            y_trans = 8
            x_lim = 12
            y_lim = 11

        x, y, drip_point_size = graph_data
        col = collections[idx]
        polys = [cyl.projected_data["XY"]["polygon"] for cyl in col.cylinders]
        stem_hull = col.stem_hull
        geopolys = GeoSeries(polys)
        geohull = GeoSeries(stem_hull)
        tot_geoHull = GeoSeries(col.hulls["XY"])

        fig, ax = plt.subplots()
        ext = "stem_highlight_xz"
        try:
            polys = [cyl.projected_data["XZ"]["polygon"] for cyl in col.cylinders]
            geopolys = GeoSeries(polys)
            colors = ["Blue" if cyl.is_stem else "Grey" for cyl in col.cylinders]
            geopolys_trans = geopolys.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
            geopolys_trans.plot(ax=ax, color=colors)
            ax.set_xlim(x_lim)
            ax.set_ylim(y_lim)
            plt.show()
            # plt.savefig(f'/code/code/canopyHydrodynamics/data/output/draw/{file}_{ext}.svg')
            # plt.savefig(
            #     f"/media/penguaman/Healthy/BranchHighlight/branchHighlight/{file}_{ext}.svg"
            # )
        except KeyError as e:
            log.error(f"error getting XZ projected data {e}")
            #   trying to run project_cylinders {e}')
            col.project_cylinders("XZ")
            polys = [cyl.projected_data["XZ"]["polygon"] for cyl in col.cylinders]

        fig, ax = plt.subplots()
        ext = "stem_highlight_xy"
        polys = [cyl.projected_data["XY"]["polygon"] for cyl in col.cylinders]
        geopolys = GeoSeries(polys)
        colors = ["Blue" if cyl.is_stem else "Grey" for cyl in col.cylinders]
        geopolys_trans = geopolys.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        geopolys_trans.plot(ax=ax, color=colors)
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        # plt.savefig(
        #     f"/media/penguaman/Healthy/BranchHighlight/branchHighlight/{file}_{ext}.svg"
        # )

        plt.show()

        fig, ax = plt.subplots()
        ext = "ghost_tree"
        geopolys_trans = geopolys.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        geohull_trans = geohull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)

        geopolys_trans.plot(ax=ax, color="lightgrey", alpha=0.7)
        geohull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)
        ax.scatter(
            [x_val + x_trans for x_val in x],
            [y_val + y_trans for y_val in y],
            facecolors="none",
            edgecolors="slategray",
            s=[x * 5 for x in drip_point_size],
        )

        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        # plt.savefig(f"./images/dripsWithTrees/{file}_{ext}.svg")
        plt.show()

        fig, ax = plt.subplots()
        ext = "two_hulls"
        tot_geoHull = GeoSeries(col.hulls["XY"])
        tot_geoHull_trans = tot_geoHull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        geohull_trans = geohull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        tot_geoHull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)
        geohull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)

        ax.scatter(
            [x_val + x_trans for x_val in x],
            [y_val + y_trans for y_val in y],
            facecolors="none",
            edgecolors="slategray",
            s=[x * 5 for x in drip_point_size],
        )

        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        plt.savefig(f"./images/dripsWithTrees/{file}_{ext}.svg")
        plt.show()

        fig, ax = plt.subplots()
        ext = "stem_hull"

        geohull_trans = geohull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        geohull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)
        ax.scatter(
            [x_val + x_trans for x_val in x],
            [y_val + y_trans for y_val in y],
            facecolors="none",
            edgecolors="slategray",
            s=[x * 5 for x in drip_point_size],
        )

        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        # plt.savefig(f"./images/dripsWithHulls/{file}_{ext}.svg")
        plt.show()

        fig, ax = plt.subplots()
        ext = "only_hulls"
        tot_geoHull = GeoSeries(col.hulls["XY"])
        tot_geoHull_trans = tot_geoHull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        geohull_trans = geohull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        tot_geoHull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)
        geohull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)

        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        # plt.savefig(f"./images//hulls/{file}_{ext}.svg")
        plt.show()


# As class functions
# def percentile_by_field(to_filter:list[dict], field, percentile):
#     to_filter_list= to_filter[field]
#     percentile_val = np.percentile(to_filter_list, percentile)
#     filter_arr = np.array(to_filter_list) >= percentile_val
#     return to_filter[filter_arr]


# def get_drip_point_coords(self, percentile):
#     """## Returns the x, y, z coordinates of drip points
#     in the percentile indicated by projected area (found in'find_flow_components')

#     ### Args:
#         - `collection (_type_)`: _description_
#         - `percentile (_type_)`: _description_

#     ### Returns:
#         - `_type_`: _description_
#     """
#     scale = 1
#     flows = np.array([flow.__dict__ for flow in self.flows])
#     flows = flows[flows["drip_node_id"] != 0]
#     drip_points = self.percentile_by_field(flows, "projected_area", percentile)
#     drip_point_locs_x = [pt[0] * scale for pt in drip_points["drip_node_loc"]]
#     drip_point_locs_y = [pt[1] * scale for pt in drip_points["drip_node_loc"]]
#     drip_point_size = [pt * 5 for pt in drip_points["projected_area"]]
#     return drip_point_locs_x, drip_point_locs_y, drip_point_size
