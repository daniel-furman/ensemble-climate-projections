from __future__ import print_function
import rasterio
import numpy as np
import os
import math
import logging
#from sklearn import metrics
#from sklearn import model_selection
logger = logging.getLogger('pyimpute')


def load_training_vector(response_shapes, explanatory_rasters, response_field, metric='mean'):
    """
    Parameters
    ----------
    response_shapes : Source of vector features for raster_stats;
                      can be OGR file path or iterable of geojson-like features
    response_field : Field name containing the known response category (must be numeric)
    explanatory_rasters : List of Paths to GDAL rasters containing explanatory variables
    metric : Statistic to aggregate explanatory data across line and polygon vector features
             Defaults to 'mean' (optional)

    Returns
    -------
    train_xs : Array of explanatory variables
    train_y : 1xN array of known responses
    """
    from rasterstats import zonal_stats
    all_means = []
    all_zones = None

    for i, raster in enumerate(explanatory_rasters):
        logger.debug("Rasters stats on %s" % raster)

        stats = zonal_stats(response_shapes, raster, stats=metric, prefix="pyimpute_", geojson_out=True)

        zones = [x['properties'][response_field] for x in stats]
        if all_zones:
            assert zones == all_zones
        else:
            all_zones = zones

        means = [x['properties']['pyimpute_' + metric] for x in stats]
        all_means.append(means)

    train_y = np.array(all_zones)
    train_xs = np.array(all_means).T

    return train_xs, train_y


def load_targets(explanatory_rasters):
    """
    Parameters
    ----------
    explanatory_rasters : List of Paths to GDAL rasters containing explanatory variables

    Returns
    -------
    expl : Array of explanatory variables
    raster_info : dict of raster info
    """

    explanatory_raster_arrays = []
    aff = None
    shape = None
    crs = None

    for raster in explanatory_rasters:
        logger.debug(raster)
        with rasterio.open(raster) as src:
            ar = src.read(1)  # TODO band num? 

            # Save or check the shape
            if not shape:
                shape = ar.shape
            else:
                assert shape == ar.shape

            # Save or check the geotransform
            if not crs:
                crs = src.crs
            else:
                assert crs == src.crs

        # Flatten in one dimension
        arf = ar.flatten()
        explanatory_raster_arrays.append(arf)

    expl = np.array(explanatory_raster_arrays).T

    raster_info = {
        'shape': shape,
        'crs': crs
    }
    return expl, raster_info


def impute(target_xs, clf, raster_info, outdir="output", linechunk=1000, class_prob=True, certainty=True):
    """
    Parameters
    ----------
    target_xs: Array of explanatory variables for which to predict responses
    clf: instance of a scikit-learn Classifier
    raster_info: dictionary of raster attributes with key 'gt', 'shape' and 'srs'

    Options
    -------
    outdir : output directory
    linechunk : number of lines to process per pass; reduce only if memory is constrained
    class_prob : Boolean. Should we create a probability raster for each class?
    certainty : Boolean. Should we produce a raster of overall classification certainty?
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    shape = raster_info['shape']

    profile = {
        'blockxsize': shape[1],
        'height': shape[0],
        'blockysize': 1,
        'count': 1,
        'crs': raster_info['crs'],
        'driver': u'GTiff',
        'dtype': 'int16',
        'nodata': -32768,
        'tiled': False,
        'width': shape[1]}

    try:
        response_path = os.path.join(outdir, "responses.tif")
        response_ds = rasterio.open(response_path, 'w', **profile)

        profile['dtype'] = 'float32'
        if certainty:
            certainty_path = os.path.join(outdir, "certainty.tif")
            certainty_ds = rasterio.open(certainty_path, 'w', **profile)

        class_dss = []
        if class_prob:
            classes = list(clf.classes_)
            class_paths = []
            for i, c in enumerate(classes):
                ods = os.path.join(outdir, "probability_%s.tif" % c)
                class_paths.append(ods)
            for p in class_paths:
                class_dss.append(rasterio.open(p, 'w', **profile))

        # Chunky logic
        if not linechunk:
            linechunk = shape[0]
        chunks = int(math.ceil(shape[0] / float(linechunk)))

        for chunk in range(chunks):
            logger.debug("Writing chunk %d of %d" % (chunk+1, chunks))
            row = chunk * linechunk
            if row + linechunk > shape[0]:
                linechunk = shape[0] - row
            # in 1D space
            start = shape[1] * row
            end = start + shape[1] * linechunk
            line = target_xs[start:end, :]

            window = ((row, row + linechunk), (0, shape[1]))

            # Predict
            responses = clf.predict(line)
            responses2D = responses.reshape((linechunk, shape[1])).astype('int16')
            response_ds.write_band(1, responses2D, window=window)

            if certainty or class_prob:
                proba = clf.predict_proba(line)

            # Certainty
            if certainty:
                certaintymax = proba.max(axis=1)
                certainty2D = certaintymax.reshape((linechunk, shape[1])).astype('float32')
                certainty_ds.write_band(1, certainty2D, window=window)

            # write out probabilities for each class as a separate raster
            for i, class_ds in enumerate(class_dss):
                proba_class = proba[:, i]
                classcert2D = proba_class.reshape((linechunk, shape[1])).astype('float32')
                class_ds.write_band(1, classcert2D, window=window)

    finally:
        response_ds.close()
        if certainty:
            certainty_ds.close()
        for class_ds in class_dss:
            class_ds.close()