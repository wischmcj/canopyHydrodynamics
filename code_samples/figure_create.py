from __future__ import annotations

import pickle

import matplotlib.pyplot as plt
import pandas as pd
from geopandas import GeoSeries

from src.canopyhydro.global_vars import log


def retain_quantile(df, field, percentile):
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val]


def return_quantile(df, field, percentile):
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val][field]


def plot_drip_points(collection, percentile):
    scale = 1
    flows = pd.DataFrame([flow.__dict__ for flow in collection.flows])
    flows = flows[flows["drip_node_id"] != 0]
    drip_points = retain_quantile(flows, "projected_area", percentile)
    drip_point_locs_x = [pt[0] * scale for pt in drip_points["drip_node_loc"]]
    drip_point_locs_y = [pt[1] * scale for pt in drip_points["drip_node_loc"]]
    drip_point_size = [pt * 5 for pt in drip_points["projected_area"]]
    return drip_point_locs_x, drip_point_locs_y, drip_point_size


if __name__ == "__main__":
    # data/output/pickle/5_SmallTree_pickle__prep_-0.1

    # load_from_pickle('Secrest32-06_000000', 'stats', 0.3666)
    # load_from_pickle('Secrest32-06_000000', 'stats', 0.3666)
    # load_from_pickle('Secrest32-06_000000_1', 'stats', -1.5)
    # load_from_pickle('5_SmallTree_1', 'stats', 0.36666)
    #    sensitivity_analysis()
    #
    # run_test_cases([('5_SmallTree',-0.16666),
    #                 ('5_SmallTree',-0.16666),
    #                 ('5_SmallTree',-0.16666)],
    #                 stats = True, fig = False, from_pickle = False)
    # run_test_cases([('Secrest32-06_000000',1.2),
    #                 ('Secrest32-06_000000',1.4),
    #                 ('Secrest32-06_000000',1.5),c
    #                 ('Secrest27-05_000000',1.2),
    #                 ('Secrest27-05_000000',1.4),
    #                 ('Secrest27-05_000000',1.5)],
    #                 stats = True, from_pickle = False)

    db_files = [
        #         open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.36', 'rb')]
        # file = '32-06_low_drip_points_36'
        #         open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_0.04', 'rb')]
        # file = '32-06_high_drip_points_04'
        #         open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.14', 'rb')]
        # file = '32-06_mid_drip_points_14
        open(
            "/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_-0.34",
            "rb",
        )
    ]
    file = "27-05_low_end_drip_points_34"
    # pick_file = 'Secrest27-05_000000'
    # pik_designation = 'stats_-0.34'
    #   open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_0.04.pickle', 'rb')]
    # file = '27-05_high_end_drip_points_04'
    #     open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_-0.14.pickle', 'rb')]
    # file = '27-05_mi_drip_points_14'

    collections = [pickle.load(db_file) for db_file in db_files]
    # collections[0].project_cylinders('XZ')
    # pickle_collection(collections[0],designation=pik_designation)

    drip_info = [plot_drip_points(col, 0.98) for col in collections]
    # max_drip = 0
    # for x,y,drip_point_size in drip_info:
    #     potential_max_drip = np.max(drip_point_size)
    #     if potential_max_drip>max_drip:
    #         max_drip = potential_max_drip\
    # 32-06
    # 10.613
    # 19.25299317980633
    # 28.45471360468393

    # 27-05
    # 6.603
    # 6.603
    # 6.603

    for x, y, drip_point_size in drip_info:
        # ensuring each drip array has the same max drip_size
        # allows us to standarize the size
        x.append(-10)
        y.append(-10)
        drip_point_size.append(28.4547)
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
        tot_geoHull = GeoSeries(col.hull)

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
            plt.savefig(
                f"/media/penguaman/Healthy/BranchHighlight/branchHighlight/{file}_{ext}.svg"
            )
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
        plt.savefig(
            f"/media/penguaman/Healthy/BranchHighlight/branchHighlight/{file}_{ext}.svg"
        )

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
        plt.savefig(
            f"/media/penguaman/Healthy/BranchHighlight/dripsWithTrees/{file}_{ext}.svg"
        )
        plt.show()

        fig, ax = plt.subplots()
        ext = "two_hulls"
        tot_geoHull = GeoSeries(col.hull)
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
        plt.savefig(
            f"/media/penguaman/Healthy/BranchHighlight/dripsWithHulls/{file}_{ext}.svg"
        )
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
        plt.savefig(
            f"/media/penguaman/Healthy/BranchHighlight/dripsWithHulls/{file}_{ext}.svg"
        )
        plt.show()

        fig, ax = plt.subplots()
        ext = "only_hulls"
        tot_geoHull = GeoSeries(col.hull)
        tot_geoHull_trans = tot_geoHull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        geohull_trans = geohull.translate(xoff=x_trans, yoff=y_trans, zoff=0.0)
        tot_geoHull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)
        geohull_trans.plot(ax=ax, color="darkgrey", alpha=0.3)

        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        plt.savefig(f"/media/penguaman/Healthy/BranchHighlight/hulls/{file}_{ext}.svg")
        plt.show()
