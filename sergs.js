// ****************************************************************************************************************** //
// ***************************** Module to implement Sequential Linear Regression Slopes (SeRGS) ******************** //
// ****************************************************************************************************************** //

/**
 * Compute The Sequential Linear Regression Slopes (SeRGS) of a satellite time series.
 * The algorithm has been tested with Sentinel-2, Landsat and MODIS data.
 * For a deeper understanding of how arrays in GEE work, check out: https://medium.com/google-earth/runs-with-arrays-400de937510a
 * The SeRGS approach is taken from: Abel, C., Horion, S., Tagesson, T., Brandt, M., & Fensholt, R. (2019).
 *                                   Towards improved remote sensing based monitoring of dryland ecosystem functioning using
 *                                   sequential linear regression slopes (SeRGS). Remote Sensing of Environment, 224, 317-332.
 * @param {ImageCollection} ts: A satellite time series in the form of an EE Image Collection.
 * @param {Array} independent_var: An array of independent variables band names.
 *                                 This can be time, precipitation, soil moisture, or any other abiotic factor.
 * @param {Array} dependent_var: An array of dependent variables band names.
 *                               This can be NDVI, NDMI, EVI, or any vegetation index.
 * @param {Dictionary} options: The options consist of the following parameters:
 *                              - @param {Number} tw_size: Size of the temporal window to apply.
 *                              - @param {Number} sw_size: Size of the spatial window to apply.
 * @returns {Image}
 * @ignore
 */
exports.sergs = function(ts, independent_var, dependent_var, options){

  // Define temporal and spatial window sizes for Sequential Linear Regression Slopes
  var tw_size = options.tw_size || 5;
  var sw_size = options.sw_size || 5;

  // Wrapper function to allow passing global parameters in a mappable function
  var wrap = function(tw_idx){
    var ts_arr = ts.select([independent_var].concat([dependent_var])).toArray();
    var ts_list = ts.select([independent_var].concat([dependent_var])).toList(ts.size());
    tw_idx = ee.Number(tw_idx);

    // Subset the time series to match the temporal window size
    var temp_subset = ts_list.slice(tw_idx.subtract(tw_size/2).add(1).int(), tw_idx.add(tw_size/2).int());
    // Initialize the first spatially-reduced image corresponding to the above temporal subset
    var init_arr = ee.Image(temp_subset.get(0)).neighborhoodToArray(ee.Kernel.square(sw_size/2));
    init_arr = init_arr.arrayReshape(ee.Image(sw_size*sw_size).toArray(), 1);

    // Loop over the rest of the temporal subsets to generate the other spatially-reduced images
    var spatial_subset = ee.Image(ee.ImageCollection(temp_subset.slice(1)).iterate(function(img, prev){
      var win_arr = img.neighborhoodToArray(ee.Kernel.square(sw_size/2));
      win_arr = win_arr.arrayReshape(ee.Image(sw_size*sw_size).toArray(), 1);

      return ee.Image(prev).arrayCat(win_arr, 0)
    }, init_arr));

    var bandNames = [independent_var, dependent_var];
    // Run the linear regression on the stack of reduced images
    var sergs_img = spatial_subset.toArray(1).toArray(2).arrayReduce(ee.Reducer.linearFit(), [0, 2], 1)
                                  .arrayProject([1])
                                  .arrayFlatten([['scale', 'offset']]);

    // Retrieve original centre image of the temporal window and its associated metadata
    var timestamp_img = ts.toList(ts.size()).get(tw_idx);
    return ee.Image(timestamp_img).addBands(sergs_img);
  };

  return wrap

};
