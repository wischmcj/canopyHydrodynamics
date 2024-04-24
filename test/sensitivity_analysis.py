import pickle
import math
import matplotlib.pyplot as plt
import numpy as np
from geopandas import GeoSeries
import pandas as pd
import matplotlib.colors as colors

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
    flows = flows[flows['drip_node_id']!=0]
    drip_points  = retain_quantile(flows, 'projected_area', percentile)
    drip_point_locs_x = [pt[0] * scale for pt in drip_points['drip_node_loc']]
    drip_point_locs_y = [pt[1] * scale for pt in drip_points['drip_node_loc']]
    drip_point_size = [pt for pt in drip_points['projected_area']]
    drip_node = [pt for pt in drip_points['drip_node_id']]
    return drip_point_locs_x,drip_point_locs_y,drip_point_size

def draw_for_paper():
    import matplotlib.pyplot as plt
    class nlcmap(object):
        def __init__(self, cmap, levels):
            self.cmap = cmapr
            self.N = cmap.N
            self.monochrome = self.cmap.monochrome
            self.levels = np.asarray(levels, dtype='float64')
            self._x = self.levels
            self.levmax = self.levels.max()
            self.transformed_levels = np.linspace(0.0, self.levmax,
                len(self.levels))

        def __call__(self, xi, alpha=1.0, **kw):
            yi = np.interp(xi, self._x, self.transformed_levels)
            return self.cmap(yi / self.levmax, alpha)
    
    
    # import numpy as np
    # import matplotlib.pyplot as plt

    # x = y = np.linspace(1, 10, 10)

    # t1mean, t2mean = 2, 9
    # sigma1, sigma2 = .3, .01
    # t1 = np.random.normal(t1mean, sigma1, 10)
    # t2 = np.random.normal(t2mean, sigma2, 10)

    class nlcmap(object):
        def __init__(self, cmap, levels):
            self.cmap = cmap
            self.N = cmap.N
            self.monochrome = self.cmap.monochrome
            self.levels = np.asarray(levels, dtype='float64')
            self._x = self.levels
            self.levmax = self.levels.max()
            self.transformed_levels = np.linspace(0.0, self.levmax,
                len(self.levels))

        def __call__(self, xi, alpha=1.0, **kw):
            yi = np.interp(xi, self._x, self.transformed_levels)
            return self.cmap(yi / self.levmax, alpha)

    # tmax = max(t1.max(), t2.max())
    # #the choice of the levels depends on the data:
    # levels = np.concatenate((
    #     [0, tmax],
    #     np.linspace(t1mean - 4 * sigma1, t1mean + 4 * sigma1, 5),
    #     np.linspace(t2mean - 4 * sigma2, t2mean + 4 * sigma2, 5),
    #     ))

    # levels = levels[levels <= tmax]
    # levels.sort()
    # # breakpoint()

    # cmap_nonlin = nlcmap(plt.cm.jet, levels)

    # fig, (ax1, ax2) = plt.subplots(1, 2)

    # ax1.scatter(x, y, edgecolors=cmap_nonlin(t1), s=15, linewidths=4)
    # ax2.scatter(x, y, edgecolors=cmap_nonlin(t2), s=15, linewidths=4)

    # fig.subplots_adjust(left=.25)
    # cbar_ax = fig.add_axes([0.10, 0.15, 0.05, 0.7])

    # #for the colorbar we map the original colormap, not the nonlinear one:
    # sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, 
    #                 norm=plt.Normalize(vmin=0, vmax=tmax))
    # sm._A = []
    # plt.show()

if __name__ == "__main__":
    # data/output/pickle/5_SmallTree_pickle__prep_-0.1

    # load_from_pickle('Secrest32-06_000000', 'stats', 0.3666)
    # load_from_pickle('Secrest32-06_000000', 'stats', 0.3666)
    # load_from_pickle('Secrest32-06_000000_1', 'stats', -1.5)
    # load_from_pickle('5_SmallTree_1', 'stats', 0.36666)
    sensitivity_analysis()       
    # draw_for_paper()
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
