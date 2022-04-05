// ****************************************************************************************************************** //
// ************************************** SeRGS example for Senegal ************************************************* //
// ****************************************************************************************************************** //

// Import external dependencies
var palettes = require('users/gena/packages:palettes');
var S2Masks = require('users/soilwatch/soilErosionApp:s2_masks.js');
var composites = require('users/soilwatch/soilErosionApp:composites.js');
var drawingTools = require('users/soilwatch/soilErosionApp:drawing_tools.js');

// Import the Dynamic Time Warping script
var sergs = require('users/soilwatch/functions:sergs.js');
var sergs_palette = palettes.colorbrewer.PiYG[11];

// Define years for which to apply SeRGS
var START_YEAR = '2000-01-01';
var END_YEAR = '2021-12-31';
var TW_SIZE = 4;
var SW_SIZE = 7;

// Define Geography for which to apply SeRGS
var country = ee.FeatureCollection("FAO/GAUL/2015/level0")
              .filter(ee.Filter.or(ee.Filter.eq('ADM0_NAME', 'Senegal'),
                                   ee.Filter.eq('ADM0_NAME', 'Gambia')));

Map.addLayer(country, {}, 'Senegal + Gambia')

// Load independent variables to use: Mean Soil Moisture and Log of Annual Precipitation
var soil_moisture = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate(START_YEAR, END_YEAR).select('soil');
soil_moisture = soil_moisture.map(function(img){return img.divide(10).unmask(0).copyProperties(img, ['system:time_start'])});
var precipitation = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD").filterDate(START_YEAR, END_YEAR).select('precipitation');
var precipitation_log = precipitation.map(function(img){return img.log().unmask(0).copyProperties(img, ['system:time_start'])});

// Landsat harmonization coefficients
var coefficients = {
  itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172])
             .multiply(10000),
  slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};

// Function to get and rename bands of interest from OLI.
function renameOli(img) {
  return img.select(
      ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa'],
      ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'pixel_qa']);
}

// Function to get and rename bands of interest from ETM+.
function renameEtm(img) {
  return img.select(
      ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'],
      ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'pixel_qa']);
}

// Function to adapt ETM+ to OLI using above-specified coefficients.
function etmToOli(img) {
  return img.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'])
      .multiply(coefficients.slopes)
      .add(coefficients.itcps)
      .round()
      .toShort()
      .addBands(img.select('pixel_qa'));
}

// Cloud Masking procedure with Landsat
function fmask(img) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = img.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask)
                 .eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return img.updateMask(mask);
}

// Function calculate NDVI and rescale to [-2000, 10000]
function calcNDVI(img) {
  return img.normalizedDifference(['B8', 'B4']).multiply(10000).rename('NDVI');
}

// Define function to prepare OLI images.
function prepOli(img) {
  var orig = img;
  img = renameOli(img);
  img = fmask(img);
  var img_ndvi = calcNDVI(img);
  return ee.Image(img.addBands(img_ndvi).copyProperties(orig, orig.propertyNames()));
}

// Define function to prepare ETM+ images.
function prepEtm(img) {
  var orig = img;
  img = renameEtm(img);
  img = fmask(img);
  img = etmToOli(img);
  var img_ndvi = calcNDVI(img);
  return ee.Image(img.addBands(img_ndvi).copyProperties(orig, orig.propertyNames()));
}

// Join time series images by year
var harmonize_ts = function(collection){
  collection = collection.map(function(img) {
    return img.set('year', ee.String(img.date().get('year'))
                  //.cat(ee.String(img.date().get('month')))
                  );
  });

  var distinctYearCol = collection.distinct('year');

  var filter = ee.Filter.equals({leftField: 'year', rightField: 'year'});
  var join = ee.Join.saveAll('year_matches');

  var joinCol = ee.ImageCollection(join.apply(distinctYearCol, collection, filter));

  return joinCol
}

// This function adds a time band to the image.
var createTimeBand = function(image) {
  // Scale milliseconds by a large constant to avoid very small slopes
  // in the linear regression output.
  return image.addBands(image.metadata('system:time_start').divide(3.154e10).add(1970));
};

// Look-up table for kendall coefficients
function get_kendall_coef(n, level){
    // The minus 4 is because the indexing below for a sample size of 4
    n = n.subtract(4);
    var coefs = {90: ee.List([4, 6, 7, 9, 10, 12, 15, 17, 18, 22, 23, 27, 28, 32, 35, 37, 40, 42,
                  45, 49, 52, 56, 59, 61, 66, 68, 73, 75, 80, 84, 87, 91, 94, 98, 103,
                  107, 110, 114, 119, 123, 128, 132, 135, 141, 144, 150, 153, 159,
                  162, 168, 173, 177, 182, 186, 191, 197, 202]),
               95: ee.List([4, 6, 9, 11, 14, 16, 19, 21, 24, 26, 31, 33, 36, 40, 43, 47, 50, 54,
                    59, 63, 66, 70, 75, 79, 84, 88, 93, 97, 102, 106, 111, 115, 120,
                    126, 131, 137, 142, 146, 151, 157, 162, 168, 173, 179, 186, 190,
                    197, 203, 208, 214, 221, 227, 232, 240, 245, 251, 258]),
               99: ee.List([6, 8, 11, 18, 22, 25, 29, 34, 38, 41, 47, 50, 56, 61, 65, 70, 76, 81,
                    87, 92, 98, 105, 111, 116, 124, 129, 135, 142, 150, 155, 163, 170,
                    176, 183, 191, 198, 206, 213, 221, 228, 236, 245, 253, 260, 268,
                    277, 285, 294, 302, 311, 319, 328, 336, 345, 355, 364])}
    return coefs[level].get(n);
}

// Derive Significance masks from the kendall coefficients
function signif_mask(img, mk_trend, period){
  // Define Kendall parameter values for a significance of 0.05
  //var period = end_year.get('year').subtract(start_year.get('year')).add(1);
  var kendall90 = ee.Number(get_kendall_coef(period, 90));
  var kendall95 = ee.Number(get_kendall_coef(period, 95));
  var kendall99 = ee.Number(get_kendall_coef(period, 99));
  // Create final productivity trajectory output layer. Positive values are
  // significant increase, negative values are significant decrease.
  return ee.Image(-32768)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall90)), 1)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall95)), 2)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall99)), 3)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall90)), -1)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall95)), -2)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall99)), -3)
        .where(mk_trend.abs().lte(kendall90), 0);
        //.where(img.abs().lte(0.001), 0).toFloat();
}

// Calculate Mann-Kendall S statistic
function mann_kendall(imageCollection){
  /***Calculate Mann Kendall's S statistic.
  //This function returns the Mann Kendall's S statistic, assuming that n is
  //less than 40. The significance of a calculated S statistic is found in
  //table A.30 of Nonparametric Statistical Methods, second edition by
  //Hollander & Wolfe.
  //Args:
  //    imageCollection: A Google Earth Engine image collection.
  //Returns:
  //    A Google Earth Engine image collection with Mann Kendall statistic for
  //        each pixel.
  */

  var afterFilter = ee.Filter.lessThan({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

  var joined = ee.ImageCollection(ee.Join.saveAll('after').apply({
    primary: imageCollection,
    secondary: imageCollection,
    condition: afterFilter
  }));

  var sign = function(i, j) { // i and j are images
    //return ee.Image(j).neq(i) // Zero case
    //    .multiply(ee.Image(j).subtract(i).clamp(-1, 1)).int();
      var concordant = ee.Image(i).lt(j).rename('concordant');
      var discordant = ee.Image(i).gt(j).rename('discordant');
      return concordant.addBands(discordant);
  };

  var mk = ee.ImageCollection(joined.map(function(current) {
    var afterCollection = ee.ImageCollection.fromImages(current.get('after'));
    return afterCollection.map(function(image) {
      // The unmask is to prevent accumulation of masked pixels that
      // result from the undefined case of when either current or image
      // is masked.  It won't affect the sum, since it's unmasked to zero.
      return ee.Image(sign(current, image)).unmask(0);
    });
    // Set parallelScale to avoid User memory limit exceeded.
  }).flatten()).reduce('sum', 4);

  //var palette = ['red', 'white', 'green'];
  //Stretch this as necessary.
  //Map.addLayer(kendall, {palette: palette}, 'kendall');

  /*
  var mk_list = imageCollection.toList(imageCollection.size().subtract(1)).map(function(CurrentImage){
    var new_col = imageCollection.filter(ee.Filter.gt('system:time_start', ee.Image(CurrentImage).get('system:time_start')));
    return new_col.toList(new_col.size()).map(function(nextImage){
      var concordant = ee.Image(CurrentImage).lt(nextImage).rename('concordant');
      var discordant = ee.Image(CurrentImage).gt(nextImage).rename('discordant');
      return concordant.addBands(discordant);
    })
  }).flatten();
  */

  //var mk_stat = ee.ImageCollection(mk_list).reduce('sum', 4);
  var mk_stat = mk.select('concordant_sum').subtract(mk.select('discordant_sum'));
  return mk_stat.toFloat()

}

// General time intervals based on a time period and interval size (in days)
var extractTimeRanges = function(start, end, agg_interval){
    /*
    Extract the time range data from the received time range and aggregation interval e.g.,
    input time interval: time_interval = ['2019-01-01','2020-01-01'], agg_interval: 60 days
    generate the following time intervals:
    time_range = [("2019-01-01T00:00:00Z", "2019-03-01T00:00:00Z"),
                ("2019-03-01T00:00:00Z", "2019-05-01T00:00:00Z"),
                ("2019-05-01T00:00:00Z", "2019-07-01T00:00:00Z"),
                ("2019-07-01T00:00:00Z", "2019-09-01T00:00:00Z"),
                ("2019-09-01T00:00:00Z", "2019-11-01T00:00:00Z"),
                ("2019-11-01T00:00:00Z", "2020-01-01T00:00:00Z")
    */

    var start_date = ee.Date(start);
    var end_date = ee.Date(end);

    // Number of intervals in the given "time_range" based on the specified "agg_interval" period
    var interval_no = ee.Date(end).difference(ee.Date(start), 'day').divide(agg_interval).round();
    var month_check = ee.Number(30.4375 / agg_interval).ceil(); // The number of aggregation intervals within a month

    // Compute the relative date delta (in months) to add to each preceding period to compute the new one
    var rel_delta = ee.Number(end_date.difference(start_date, 'day'))
                    .divide(ee.Number(30.4375).multiply(interval_no)).ceil(); // 30.4375 days = average month length

    // Compute the first time interval end date by adding the relative date delta (in months) to the start date
    end_date = start_date.advance(start_date.advance(rel_delta, 'month')
                                  .difference(start_date, 'day')
                                  .divide(month_check), 'day')
                                  .advance(-1, 'second');

    var time_intervals = ee.List([ee.List([start_date, end_date])]);
    time_intervals = ee.List(ee.List.sequence(1, interval_no.subtract(1)).iterate(function(x,previous){
        start_date = ee.Date(ee.List(ee.List(previous).reverse().get(0)).get(1))
                     .advance(1, 'second'); //end_date of last element
        end_date = start_date
                   .advance(start_date.advance(rel_delta, 'month')
                   .difference(start_date, 'day')
                   .divide(month_check), 'day')
                   .advance(-1, 'second');

        return ee.List(previous).add(ee.List([start_date, end_date]));
    }, time_intervals));

    return time_intervals;
}

// Merge landsat time series from different sensors (TM, ETM+, OLI)
function merge_landsat(geom){
  var oliCol = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate(START_YEAR, END_YEAR);
  var etmCol = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterDate(START_YEAR, END_YEAR);
  var tmCol = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterDate(START_YEAR, END_YEAR);

  var colFilter = ee.Filter.and(
      ee.Filter.bounds(geom),
      ee.Filter.lt('CLOUD_COVER', 50), ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
      ee.Filter.or(
          ee.Filter.eq('IMAGE_QUALITY', 9),
          ee.Filter.eq('IMAGE_QUALITY_OLI', 9))
      );

  // Filter collections and prepare them for merging.
  oliCol = oliCol.filter(colFilter).map(prepOli);
  etmCol = etmCol.filter(colFilter).map(prepEtm);
  tmCol = tmCol.filter(colFilter).map(prepEtm);

  // Merge the collections.
  var col = oliCol.merge(etmCol).merge(tmCol);

  return col
}

// Apply SeRGS by specifying an independent and dependent variable, as well as a temporal and spatial window size
// Unfortunately, multivariate linear regressions using ee.LinearRegression() is not possible using EE arrays.
// Only ee.LinearFit() is applied with a single independent and dependent variable.
function apply_sergs(independent_col, dependent_col, independent_var, dependent_var, tw_size, sw_size){

  // Aggregate annual composites for the dependent variable, in this case NDVI time series of Landsat or MODIS
  var dependentComp = harmonize_ts(dependent_col).map(function(img) {
    var yearCol = ee.ImageCollection.fromImages(img.get('year_matches'));

    // Calculation of the NDVI small integral as per Eklundh et al., 2004.
    // Percentiles are preferred to actual min and max to avoid noisy/outlier pixels.
    // Moreover, the mean of the values corresponding to the upper half of the NDVI amplitude is taken,
    // as a simplification to the actual integral
    var min_col = yearCol.select(['NDVI']).reduce(ee.Reducer.percentile([5]), 4)
                  .set('system:time_start', ee.Date(img.date().update({year: img.date().get('year'),
                                                                       month: 7,
                                                                       day: 1, hour:0, minute:0, second:0})).millis());

    var max_col = yearCol.reduce(ee.Reducer.percentile([95]), 4)
                        .set('system:time_start', ee.Date(img.date().update({year: img.date().get('year'),
                                                                             month: 7,
                                                                             day: 1, hour:0, minute:0, second:0})).millis());

    var amplitude = max_col.select(['NDVI_p95']).subtract(min_col);

    var mean_col = yearCol.select(['NDVI']).map(function(img){return img.updateMask(img.gte(amplitude.divide(2)))}).reduce(ee.Reducer.mean(), 4)
    //var mean_col = yearCol.reduce(ee.Reducer.mean(), 4)
                   .set('system:time_start', ee.Date(img.date().update({year: img.date().get('year'),
                                                                       month: 7,
                                                                       day: 1, hour:0, minute:0, second:0})).millis());

    return mean_col.rename(['NDVI_integ'])
  })
  .map(createTimeBand);

  // Generate a list of time intervals for which to generate a harmonized time series
  var time_intervals = extractTimeRanges(START_YEAR, END_YEAR, 364);

  // Generate harmonized annual time series of the independent variable
  var independentComp = composites.harmonizedTS(independent_col, [independent_var], time_intervals, {agg_type: 'sum'})
                          .map(function(img){return img.set('system:time_start', ee.Date(img.date().update({year: img.date().get('year'),
                                                                                 month: 7,
                                                                                 day: 1, hour:0, minute:0, second:0})).millis());
                                                                                 });
  // Join the dependent and independent time series using 'system:time_start'
  dependentComp = ee.ImageCollection(ee.Join.saveFirst(independent_var).apply({
      primary: dependentComp,
      secondary: independentComp,
      condition: ee.Filter.equals({
          leftField: 'system:time_start',
          rightField: 'system:time_start'})
  }))
  .map(function(img){return img.addBands(ee.Image(img.get(independent_var)));})
  .sort('system:time_start');

  /*
  // Apply the RESTREND methodology as per Wessel et al., 2007
  // Compute robust linear regression coefficients.
  var vpd_ols = dependentComp.select(['precipitation_sum','NDVI_integ'])
                .reduce(ee.Reducer.linearFit(), 4);

  var dependentComp_restrend = dependentComp.map(function(img){
    var pred_ndvi = img.select('precipitation_sum')
                    .multiply(vpd_ols.select('scale'))
                    .add(vpd_ols.select('offset'));
    var obs_ndvi = img.select('NDVI_integ');

    return img.select('system:time_start', 'precipitation_sum', 'NDVI_integ')
              .addBands(obs_ndvi.subtract(pred_ndvi).rename('NDVI_residual'))
  });
  */

  // Apply SeRGS with the pre-processed and joined time series containing the independent and dependent variable
  var dependentComp_sergs = ee.ImageCollection(ee.List.sequence(Math.ceil(tw_size/2-1), dependentComp.size().subtract(parseInt(tw_size/2+1)))
             .map(sergs.sergs(dependentComp,
                             independent_var,
                             dependent_var,
                             {tw_size: tw_size, sw_size: sw_size})
                  )
    )

  return dependentComp_sergs
}

// Apply the sequence of functions to get to the mann-kendall significance masked SeRGS
var landsat_ts = merge_landsat(country.geometry());
var prec_landsat_sergs = apply_sergs(precipitation_log, landsat_ts, 'precipitation', 'NDVI_integ', TW_SIZE, SW_SIZE);
var prec_landsat_mk = mann_kendall(prec_landsat_sergs.select('scale')).rename('mannkendall');
prec_landsat_mk = ee.Image('users/SoilWatch/Senegal/landsat_prec_ndvi_sergs_2000_2021').select('mannkendall');
var prec_landsat_ols =  prec_landsat_sergs.select(['system:time_start', 'scale'])
                                  .reduce(ee.Reducer.linearFit(), 4).select('scale');
prec_landsat_ols = ee.Image('users/SoilWatch/Senegal/landsat_prec_ndvi_sergs_2000_2021').select('scale');
var prec_landsat_signif = signif_mask(prec_landsat_ols.select('scale'), prec_landsat_mk, prec_landsat_sergs.size()).rename('significance');
prec_landsat_signif = ee.Image('users/SoilWatch/Senegal/landsat_prec_ndvi_sergs_2000_2021').select('significance');

Map.centerObject(country);
Map.setOptions("SATELLITE");

Map.addLayer(prec_landsat_signif.updateMask(prec_landsat_signif.neq(0)).clip(country.geometry()),
             {min: -3, max:3, palette: sergs_palette},
             'SeRGS of Log-Rainfall-NDVI 2000-2021, Landsat')

// Apply the sequence of functions to get to the mann-kendall significance masked SeRGS
var modis_ts = ee.ImageCollection("MODIS/006/MOD13Q1").filterDate(START_YEAR, END_YEAR).filterBounds(country.geometry());
var prec_modis_sergs = apply_sergs(precipitation_log, modis_ts, 'precipitation', 'NDVI_integ', TW_SIZE, SW_SIZE);
var prec_modis_mk = mann_kendall(prec_modis_sergs.select('scale')).rename('mannkendall');
prec_modis_mk = ee.Image('users/SoilWatch/Senegal/modis_prec_ndvi_sergs_2000_2021').select('mannkendall');
var prec_modis_ols =  prec_modis_sergs.select(['system:time_start', 'scale'])
                                  .reduce(ee.Reducer.linearFit(), 4).select('scale');
prec_modis_ols = ee.Image('users/SoilWatch/Senegal/modis_prec_ndvi_sergs_2000_2021').select('scale');
var prec_modis_signif = signif_mask(prec_modis_ols.select('scale'), prec_modis_mk, prec_modis_sergs.size()).rename('significance');
prec_modis_signif = ee.Image('users/SoilWatch/Senegal/modis_prec_ndvi_sergs_2000_2021').select('significance');

Map.addLayer(prec_modis_signif.updateMask(prec_modis_signif.neq(0)).clip(country.geometry()),
             {min: -3, max:3, palette: sergs_palette},
             'SeRGS of Log-Rainfall-NDVI 2000-2021, MODIS')

// Apply the sequence of functions to get to the mann-kendall significance masked SeRGS
var ssm_landsat_sergs = apply_sergs(soil_moisture, landsat_ts, 'soil', 'NDVI_integ', TW_SIZE, SW_SIZE);
var ssm_landsat_mk = mann_kendall(ssm_landsat_sergs.select('scale')).rename('mannkendall');
ssm_landsat_mk = ee.Image('users/SoilWatch/Senegal/landsat_ssm_ndvi_sergs_2000_2021').select('mannkendall');
var ssm_landsat_ols =  ssm_landsat_sergs.select(['system:time_start', 'scale'])
                                  .reduce(ee.Reducer.linearFit(), 4).select('scale');
ssm_landsat_ols = ee.Image('users/SoilWatch/Senegal/landsat_ssm_ndvi_sergs_2000_2021').select('scale');
var ssm_landsat_signif = signif_mask(ssm_landsat_ols.select('scale'), ssm_landsat_mk, ssm_landsat_sergs.size()).rename('significance');
ssm_landsat_signif = ee.Image('users/SoilWatch/Senegal/landsat_ssm_ndvi_sergs_2000_2021').select('significance');

Map.addLayer(ssm_landsat_signif.updateMask(ssm_landsat_signif.neq(0)).clip(country.geometry()),
             {min: -3, max:3, palette: sergs_palette},
             'SeRGS of Soil Moisture-NDVI 2000-2021, Landsat')

// Apply the sequence of functions to get to the mann-kendall significance masked SeRGS
var ssm_modis_sergs = apply_sergs(soil_moisture, modis_ts, 'soil', 'NDVI_integ', TW_SIZE, SW_SIZE);
var ssm_modis_mk = mann_kendall(ssm_modis_sergs.select('scale')).rename('mannkendall');
ssm_modis_mk = ee.Image('users/SoilWatch/Senegal/modis_ssm_ndvi_sergs_2000_2021').select('mannkendall');
var ssm_modis_ols =  ssm_modis_sergs.select(['system:time_start', 'scale'])
                                  .reduce(ee.Reducer.linearFit(), 4).select('scale');
ssm_modis_ols = ee.Image('users/SoilWatch/Senegal/modis_ssm_ndvi_sergs_2000_2021').select('scale');
var ssm_modis_signif = signif_mask(ssm_modis_ols.select('scale'), ssm_modis_mk, ssm_modis_sergs.size()).rename('significance');
ssm_modis_signif = ee.Image('users/SoilWatch/Senegal/modis_ssm_ndvi_sergs_2000_2021').select('significance');

Map.addLayer(ssm_modis_signif.updateMask(ssm_modis_signif.neq(0)).clip(country.geometry()),
             {min: -3, max:3, palette: sergs_palette},
             'SeRGS of Soil Moisture-NDVI 2000-2021, MODIS')

// Define arguments for animation function parameters.
var video_args = {
  'dimensions': 2500,
  'region': country.geometry(),
  'framesPerSecond': 1,
  'crs': 'EPSG:4326',
  'min': -3,
  'max': 3,
  'palette': sergs_palette
}

// Function to export GIF
var generateGIF = function(img, band_name, args){
  print(img.select(band_name).getVideoThumbURL(args));
}

// Export GIF of the DTW classification + DTW classification distinguishing between abandoned/active cropland
print('Log-Rainfall-NDVI (Landsat/MODIS) SeRGS GIF:');
generateGIF(ee.ImageCollection(ee.List([prec_landsat_signif.updateMask(prec_landsat_signif.neq(0))])
               .add(prec_modis_signif.updateMask(prec_modis_signif.neq(0)))), 'significance', video_args);

// Export GIF of the DTW classification + DTW classification distinguishing between abandoned/active cropland
print('Soil Moisture-NDVI (Landsat/MODIS) SeRGS GIF:');
generateGIF(ee.ImageCollection(ee.List([ssm_landsat_signif.updateMask(ssm_landsat_signif.neq(0))])
               .add(ssm_modis_signif.updateMask(ssm_modis_signif.neq(0)))), 'significance', video_args);

// Export Video of the temporal NDVI composites for each year produced
Export.video.toDrive({
  collection: ee.ImageCollection(prec_landsat_sergs.map(function(img){return ee.Image(img).select('NDVI_integ').divide(100)
                                                                           .addBands(ee.Image(img).select('precipitation').divide(10))
                                                                           .addBands(ee.Image(img).select('scale').add(10).multiply(10))
                                                                           .toByte()
  })).select(["scale", "NDVI_integ", 'precipitation']),
  description:'prec_ndvi_rgb_landsat',
  region: country.geometry(),
  scale: 250,
  crs:'EPSG:4326',
  maxPixels: 1e13
});

// Export Video of the temporal NDVI composites for each year produced
Export.video.toDrive({
  collection: ee.ImageCollection(prec_modis_sergs.map(function(img){return ee.Image(img).select('NDVI_integ').divide(100)
                                                                           .addBands(ee.Image(img).select('precipitation').divide(10))
                                                                           .addBands(ee.Image(img).select('scale').add(10).multiply(10))
                                                                           .toByte()
  })).select(["scale", "NDVI_integ", 'precipitation']),
  description:'prec_ndvi_rgb_modis',
  region: country.geometry(),
  scale: 250,
  crs:'EPSG:4326',
  maxPixels: 1e13
});

// Export Video of the temporal NDVI composites for each year produced
Export.video.toDrive({
  collection: ee.ImageCollection(ssm_landsat_sergs.map(function(img){return ee.Image(img).select('NDVI_integ').divide(100)
                                                                           .addBands(ee.Image(img).select('soil'))
                                                                           .addBands(ee.Image(img).select('scale').add(10).multiply(10))
                                                                           .toByte()
  })).select(["scale", "NDVI_integ", 'soil']),
  description:'ssm_ndvi_rgb_landsat',
  region: country.geometry(),
  scale: 250,
  crs:'EPSG:4326',
  maxPixels: 1e13
});

// Export Video of the temporal NDVI composites for each year produced
Export.video.toDrive({
  collection: ee.ImageCollection(ssm_modis_sergs.map(function(img){return ee.Image(img).select('NDVI_integ').divide(100)
                                                                           .addBands(ee.Image(img).select('soil'))
                                                                           .addBands(ee.Image(img).select('scale').add(10).multiply(10))
                                                                           .toByte()
  })).select(["scale", "NDVI_integ", 'soil']),
  description:'ssm_ndvi_rgb_modis',
  region: country.geometry(),
  scale: 250,
  crs:'EPSG:4326',
  maxPixels: 1e13
});

// Export classified data. This is recommended to get the full extent of the data generated and saved,
// so it can be explored and consulted seamlessly.
Export.image.toAsset({
  image: ssm_landsat_ols.toFloat().addBands(ssm_landsat_mk).addBands(ssm_landsat_signif).clip(country.geometry()),
  description:'landsat_ssm_ndvi_sergs_2000_2021',
  assetId: 'users/SoilWatch/Senegal/landsat_ssm_ndvi_sergs_2000_2021',
  region: country.geometry(),
  crs: 'EPSG:4326',
  scale: 30,
  maxPixels:1e13
});

// Export classified data. This is recommended to get the full extent of the data generated and saved,
// so it can be explored and consulted seamlessly.
Export.image.toAsset({
  image: ssm_modis_ols.toFloat().addBands(ssm_modis_mk).addBands(ssm_modis_signif).clip(country.geometry()),
  description:'modis_ssm_ndvi_sergs_2000_2021',
  assetId: 'users/SoilWatch/Senegal/modis_ssm_ndvi_sergs_2000_2021',
  region: country.geometry(),
  crs: 'EPSG:4326',
  scale: 250,
  maxPixels:1e13
});

// Export classified data. This is recommended to get the full extent of the data generated and saved,
// so it can be explored and consulted seamlessly.
Export.image.toAsset({
  image: prec_landsat_ols.toFloat().addBands(prec_landsat_mk).addBands(prec_landsat_signif).clip(country.geometry()),
  description:'landsat_prec_ndvi_sergs_2000_2021',
  assetId: 'users/SoilWatch/Senegal/landsat_prec_ndvi_sergs_2000_2021',
  region: country.geometry(),
  crs: 'EPSG:4326',
  scale: 30,
  maxPixels:1e13
});

// Export classified data. This is recommended to get the full extent of the data generated and saved,
// so it can be explored and consulted seamlessly.
Export.image.toAsset({
  image: prec_modis_ols.toFloat().addBands(prec_modis_mk).addBands(prec_modis_signif).clip(country.geometry()),
  description:'modis_prec_ndvi_sergs_2000_2021',
  assetId: 'users/SoilWatch/Senegal/modis_prec_ndvi_sergs_2000_2021',
  region: country.geometry(),
  crs: 'EPSG:4326',
  scale: 250,
  maxPixels:1e13
});

// Initiate drawing tools functionalities
var drawing_elements =  drawingTools.initializeDrawingTools();
var drawing_tools = drawing_elements[0];
var control_panel = drawing_elements[1];

Map.add(control_panel);

drawing_tools.onDraw(ui.util.debounce(chartCustomTimeSeries, 500));

function chartCustomTimeSeries(){
  // Get the drawn geometry; it will define the reduction region.
  var aoi = drawing_tools.layers().get(0).getEeObject();
  var short_date_range = ee.Dictionary({'start': START_YEAR+'-01-01', 'end': END_YEAR+'-12-31'});

    // Set the drawing mode back to null; turns drawing off.
  drawing_tools.setShape(null);

  // Reduction scale is based on map scale to avoid memory/timeout errors.
  var map_scale = Map.getScale();
  var scale = map_scale > 250 ? 250 : 90;

    // Generate plot title
  var pt_title =  ' - centroid coordinates (lon/lat): '
                  + ee.Number(aoi.centroid(0.001).coordinates().get(0)).multiply(1e6).round().divide(1e6).getInfo() + ', '
                  + ee.Number(aoi.centroid(0.001).coordinates().get(1)).multiply(1e6).round().divide(1e6).getInfo();

  // Define LandTrendr parameters.
  var ltParams = {
    maxSegments: 6,
    spikeThreshold: 0.5,
    vertexCountOvershoot: 3,
    preventOneYearRecovery: false,
    recoveryThreshold: 0.25,
    pvalThreshold: 0.05,
    bestModelProportion: 0.75,
    minObservationsNeeded: 6
  };

  ltParams.timeSeries = prec_landsat_sergs.select(['scale']);

  var ltArrImg = ee.Algorithms.TemporalSegmentation.LandTrendr(ltParams);

  // Reduce the LandTrendr result by the aoi. Use ee.Reducer.first() to select
  // the pixel that intersects the point.
  var ltPoints = ltArrImg.reduceRegions({
    collection: aoi,
    reducer: ee.Reducer.first(),
    scale: scale,
    tileScale: 4
  });

  // Get the LandTrendr segmentation results from point: ID 0.
  var ltPoint = ee.Array(
    ltPoints.first().get('LandTrendr')
  );

  // Slice out data to plot in time series chart.
  var yearValues = ltPoint.slice(0, 0, 1).transpose();
  var yValues = ltPoint.slice(0, 1, 3).transpose();

  // Make a time series chart.
  var LTchart_sergs = ui.Chart.array.values(yValues, 0, yearValues)
    //.setChartType('LineChart')
    .setSeriesNames(['NDVI_sergs', 'NDVI_fitted'])
    .setOptions({
      title: 'Log-Rainfall-Landsat NDVI SeRGS Trend for ' + START_YEAR.slice(0,4) + '-' + END_YEAR.slice(0,4) + ', ' + pt_title,
      trendlines: {0:{type: 'linear', color: 'black', visibleInLegend: true}
      },
      colors: ['#33b3a6', 'red'],
      hAxis: {title: 'Date'},
      vAxis: {title: 'NDVI SeRGS'},
      lineWidth: 2,
      pointSize: 0,
    });

  print(LTchart_sergs);

  ltParams.timeSeries = prec_modis_sergs.select(['scale']);

  ltArrImg = ee.Algorithms.TemporalSegmentation.LandTrendr(ltParams);

  // Reduce the LandTrendr result by the aoi. Use ee.Reducer.first() to select
  // the pixel that intersects the point.
  ltPoints = ltArrImg.reduceRegions({
    collection: aoi,
    reducer: ee.Reducer.first(),
    scale: scale,
    tileScale: 4
  });

  // Get the LandTrendr segmentation results from point: ID 0.
  ltPoint = ee.Array(
    ltPoints.first().get('LandTrendr')
  );

  // Slice out data to plot in time series chart.
  yearValues = ltPoint.slice(0, 0, 1).transpose();
  yValues = ltPoint.slice(0, 1, 3).transpose();

  // Make a time series chart.
  LTchart_sergs = ui.Chart.array.values(yValues, 0, yearValues)
    //.setChartType('LineChart')
    .setSeriesNames(['NDVI_sergs', 'NDVI_fitted'])
    .setOptions({
      title: 'Log-Rainfall-MODIS NDVI SeRGS Trend for ' + START_YEAR.slice(0,4) + '-' + END_YEAR.slice(0,4) + ', ' + pt_title,
      trendlines: {0:{type: 'linear', color: 'black', visibleInLegend: true}
      },
      colors: ['#33b3a6', 'red'],
      hAxis: {title: 'Date'},
      vAxis: {title: 'NDVI SeRGS'},
      lineWidth: 2,
      pointSize: 0,
    });

  print(LTchart_sergs);

  ltParams.timeSeries = ssm_landsat_sergs.select(['scale']);

  ltArrImg = ee.Algorithms.TemporalSegmentation.LandTrendr(ltParams);

  // Reduce the LandTrendr result by the aoi. Use ee.Reducer.first() to select
  // the pixel that intersects the point.
  ltPoints = ltArrImg.reduceRegions({
    collection: aoi,
    reducer: ee.Reducer.first(),
    scale: scale,
    tileScale: 4
  });

  // Get the LandTrendr segmentation results from point: ID 0.
  ltPoint = ee.Array(
    ltPoints.first().get('LandTrendr')
  );

  // Slice out data to plot in time series chart.
  yearValues = ltPoint.slice(0, 0, 1).transpose();
  yValues = ltPoint.slice(0, 1, 3).transpose();

  // Make a time series chart.
  LTchart_sergs = ui.Chart.array.values(yValues, 0, yearValues)
    //.setChartType('LineChart')
    .setSeriesNames(['NDVI_sergs', 'NDVI_fitted'])
    .setOptions({
      title: 'Soil Moisture-Landsat NDVI SeRGS Trend for ' + START_YEAR.slice(0,4) + '-' + END_YEAR.slice(0,4) + ', ' + pt_title,
      trendlines: {0:{type: 'linear', color: 'black', visibleInLegend: true}
      },
      colors: ['#9b7653', 'red'],
      hAxis: {title: 'Date'},
      vAxis: {title: 'NDVI SeRGS'},
      lineWidth: 2,
      pointSize: 0,
    });

  print(LTchart_sergs);

  ltParams.timeSeries = ssm_modis_sergs.select(['scale']);

  ltArrImg = ee.Algorithms.TemporalSegmentation.LandTrendr(ltParams);

  // Reduce the LandTrendr result by the aoi. Use ee.Reducer.first() to select
  // the pixel that intersects the point.
  ltPoints = ltArrImg.reduceRegions({
    collection: aoi,
    reducer: ee.Reducer.first(),
    scale: scale,
    tileScale: 4
  });

  // Get the LandTrendr segmentation results from point: ID 0.
  ltPoint = ee.Array(
    ltPoints.first().get('LandTrendr')
  );

  // Slice out data to plot in time series chart.
  yearValues = ltPoint.slice(0, 0, 1).transpose();
  yValues = ltPoint.slice(0, 1, 3).transpose();

  // Make a time series chart.
  LTchart_sergs = ui.Chart.array.values(yValues, 0, yearValues)
    //.setChartType('LineChart')
    .setSeriesNames(['NDVI_sergs', 'NDVI_fitted'])
    .setOptions({
      title: 'Soil Moisture-MODIS NDVI SeRGS Trend for ' + START_YEAR.slice(0,4) + '-' + END_YEAR.slice(0,4) + ', ' + pt_title,
      trendlines: {0:{type: 'linear', color: 'black', visibleInLegend: true}
      },
      colors: ['#9b7653', 'red'],
      hAxis: {title: 'Date'},
      vAxis: {title: 'NDVI SeRGS'},
      lineWidth: 2,
      pointSize: 0,
    });

  print(LTchart_sergs);

}

var makeLegend = function(title, palette, class_names, class_length){
  // Create a legend for the different crop types
  // set position of panel
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '12px 15px'
    }
  });

  // Create legend title
  var legendTitle = ui.Label({
    value: title,
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
      }
  });

  // Add the title to the panel
  legend.add(legendTitle);

  // Creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
        // Create the label that is actually the colored box.
        var colorBox = ui.Label({
          style: {
            backgroundColor: color,
            // Use padding to give the box height and width.
            padding: '8px',
            fontSize: '12px',
            margin: '0 0 4px 0'
          }
        });

        // Create the label filled with the description text.
        var description = ui.Label({
          value: name,
          style: {margin: '0 0 4px 6px'}
        });

        // return the panel
        return ui.Panel({
          widgets: [colorBox, description],
          layout: ui.Panel.Layout.Flow('horizontal')
        });
  };

  // Add color and and names
  for (var i = 0; i <= class_length; i++) {
    legend.add(makeRow(palette[i], class_names[i]));
    }

  return legend
}

var class_names = ['decrease (p <= 0.1)', 'decrease (p <= 0.05)', 'decrease (p <= 0.1)', 'not significant (p > 0.1)',
                   'increase (p <= 0.1)', 'increase (p <= 0.05)', 'increase (p <= 0.01)']
var palette = ['c51b7d', 'e9a3c9', 'fde0ef', 'f7f7f7', 'e6f5d0', 'a1d76a', '4d9221'];

var sergs_legend = makeLegend('Linear trend in SeRGS', palette, class_names, 6);
Map.add(sergs_legend);
