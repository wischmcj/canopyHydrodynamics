# How To

This section contains a set of example workflows that demonstrate different capabilities of `canoPyHydro`.

```{toctree}
:maxdepth: 1
creating_figures.ipynb
```

You can get the tutorials by cloning the canoPyHydro repo and navigating to the "examples" folder:

```bash
git clone https://github.com/wischmcj/canopyHydrodynamics.git
cd canopyHydrodynamics/code_samples/examples
```

## Tutorials

  The below code can be run at the first breakpoint in the test_collection_integration.py file

### Displaying, Filtering and Highlighting
```{python}
    flexible_collection.draw(plane = 'XZ')
    flexible_collection.draw(plane = 'XZ', a_lambda = lambda: cyl_id>100)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>50)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>75)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>75, highlight_lambda = lambda:branch_order==2)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100, highlight_lambda = lambda:branch_order==2)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100, highlight_lambda = lambda:is_stem)
```
### Draw all projections
```{python}
    import geopandas as geo  # only import what we need
    import matplotlib.pyplot as plt
    happy_path_projection.project_cylinders('XY')
    happy_path_projection.project_cylinders('XZ')
    happy_path_projection.project_cylinders('YZ')
    xz_poly = [cyl.projected_data['XZ']['polygon'] for cyl in happy_path_projection.cylinders[1:20]]
    xy_poly = [cyl.projected_data['XY']['polygon'] for cyl in happy_path_projection.cylinders[1:20]]
    yz_poly = [cyl.projected_data['YZ']['polygon'] for cyl in happy_path_projection.cylinders[1:20]]
    geoPolys_xy = geo.GeoSeries(xy_poly)
    geoPolys_xz = geo.GeoSeries(xz_poly)
    geoPolys_yz = geo.GeoSeries(yz_poly)
    fig, ax = plt.subplots(3)
    geoPolys_xy.plot(ax=ax[0,0])
    geoPolys_xy.plot(ax=ax[0])
    geoPolys_xz.plot(ax=ax[1])
    geoPolys_yz.plot(ax=ax[2])
```