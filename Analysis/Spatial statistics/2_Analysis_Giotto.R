
# A- Configuring the Giotto Environment
# Ensure Giotto Suite is installed.
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}

library(Giotto)
library(ggplot2)
library(future)

#removeGiottoEnvironment()
#installGiottoEnvironment()

# B- Package Accessibility

# Creating Giotto Instructions without specifying a Python path will make
# reticulate activate the default Giotto environment.
#my_python_path = "/Users/joaoluizsfilho/Library/r-miniconda/envs/giotto_env/bin/python" #if using the r-miniconda env, it does not have 3 packages
#my_python_path = "/Users/joaoluizsfilho/opt/anaconda3/envs/giotto_env/bin/pythonw" #if using the conda env
default_instrs <- createGiottoInstructions()

# Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}

# Extract python path information
default_python_path <- default_instrs$python_path

# Make reticulate iteratively check for the packages
pkg_check <- function(){
  py_pkgs = c('pandas','networkx', 'igraph', 'leidenalg','community','sklearn','python.app','scanpy')
  py_pkg_error = character()
  test_availability = TRUE
  
  for (i in py_pkgs){
    if(i == 'python.app' & Sys.info()[['sysname']] != "Darwin"){
      # If the machine OS is not OSX (Mac), break out of the loop
      # Otherwise, also check for python.app
      break
    }
    test_availability <- reticulate::py_module_available(i)
    if(!test_availability) {py_pkg_error <- c(py_pkg_error,i)}
  }
  
  if(test_availability){
    cat('All Python packages for Giotto are accessible at environment:\n', default_python_path)
  }else{
    for (x in py_pkg_error) cat(x,'was not found within environment:\n',default_python_path,'\n\n')
  }
  
  return(py_pkg_error)
}

pkg_check()

#library(reticulate)
# indicate that we want to use a specific virtualenv
#use_miniconda(condaenv = "/Users/joaoluizsfilho/Library/r-miniconda/envs/giotto_env/bin/pythonw", required = NULL)
# install SciPy
#conda_install("/Users/joaoluizsfilho/Library/r-miniconda/envs/giotto_env/bin/pythonw", "scanpy")
#py_install("scanpy")

# 1- Setup

# Custom color palettes from rcartocolor
# pal10 = rcartocolor::carto_pal(n = 10, name = 'Pastel')
pal10 = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#FF0000",   #87C55F
          "#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B3B3B3", "black")
# viv10 = rcartocolor::carto_pal(n = 10, name = 'Vivid')
viv10 = c("#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0",
          "#24796C","#DAA51B","#2F8AC4","#764E9F","#A5AA99")

# set working directory
results_folder = "./Plots/"

# Optional: Specify a path to a Python executable within a conda or miniconda environment. 
# If set to NULL (default), the Python executable within the previously installed Giotto environment will be used.
my_python_path = NULL # alternatively, "/local/python/path/python" if desired.

## Set object behavior
# by directly saving plots, but not rendering them you will save a lot of time
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  return_plot = FALSE,
                                  python_path = my_python_path)


# 15. Saving and Loading the giotto object
# Giotto uses many objects that include pointers to information that live on disk instead of loading everything into memory. 
# This includes both giotto image objects (giottoImage, giottoLargeImage) and also subcellular information (giottoPoints, giottoPolygon).
# When saving the project as a .RDS or .Rdata, these pointers are broken and can produce errors when loaded again.
# saveGiotto() is a function that can save Giotto Suite projects into a contained structured directory that can then be properly loaded again later using loadGiotto().

saveGiotto(gobject = COVID_gobject,
           dir = './COVID_Giotto_object', overwrite = TRUE)


COVID_gobject <- loadGiotto('./COVID_Giotto_object/saveGiottoDir/')


# 1.1- CosMx Project loading function

 #Convenience function for loading in the CosMx data. 
# It loads subcellular transcript information and polygons and generates a giotto object with giottoPoints objects for both ‘rna’ and ‘neg_probe’ nested in the gobject feat_info slot, 
# and a giottoPolygon object for the ‘cell’ spatial unit in the spatial_info slot.
# This function performs the manual object creation steps described below.
# To skip those steps and preliminary data exploration, go to Section 5.
# Additionally, a comparison of the count matrix produced through the convenience function ‘subcellular’ workflow and Nanostring’s provided matrix can be found at Section 6.4.

## provide path to nanostring folder
data_path = './Files/'

## create giotto cosmx object
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular', # only subcellular
                                   FOVs = c(01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25), #looks like it can load several fovs at once
                                   instructions = instrs)

# [Expected Directory] This function generates a giotto object when given a link to a cosmx output directory. 
# It expects the following items within the directory where the bolded portions are what this function matches against:
# CellComposite (folder of images)
# CellLabels (folder of images)
# CellOverlay (folder of images)
# CompartmentLabels (folder of images)
# experimentname_exprMat_file.csv (file)
# experimentname_fov_positions_file.csv (file)
# experimentname_metadata_file.csv (file)
# experimentname_tx_file.csv (file)

# ‘all’ - loads and requires subcellular information from tx_file and fov_positions_file and also the existing aggregated information (expression, spatial locations, and metadata) from exprMat_file and metadata_file.
# ‘subcellular’ - loads and requires subcellular information from tx_file and fov_positions_file only.
# ‘aggregate’ - loads and requires the existing aggregate information (expression, spatial locations, and metadata) from exprMat_file and metadata_file.

# [Images] Images in the default CellComposite, CellLabels, CompartmentLabels, and CellOverlay folders will be loaded as giotto largeImage objects in all workflows as long as they are available. 
# Additionally, CellComposite images will be converted to giotto image objects, making plotting with these image objects more responsive when accessing them from a server. 
# showGiottoImageNames can be used to see the available images.


showGiottoFeatInfo(fov_join)
showGiottoSpatialInfo(fov_join)

# 1.2- Creating Giotto object from Anndata

# Reference: https://giottosuite.readthedocs.io/en/latest/subsections/datasets/interoperability_04122023.html
# To convert an AnnData Object back into a Giotto object, it must first be saved as a .h5ad file. 
# The name of said file may then be provided to anndataToGiotto() for conversion.
# If a nearest neighbor network or spatial netowkr was created using the key_added argument, they may be provided to arguments n_key_added and/or spatial_n_key_added, respectively.

fov_integrated_gobject <- anndataToGiotto(anndata_path = "./adata_COVID_giotto",
                                  python_path = my_python_path, deluanay_spat_net = FALSE)

#aggregate_rna_gobject <- anndataToGiotto(anndata_path = "./giotto_anndata_conversion/aggregate_rna_converted_gobject.h5ad",
 #                                        python_path = my_python_path,
  #                                       n_key_added = list("sNN.pca","new_network"),
   #                                      spatial_n_key_added = "aggregate_rna_spatial_network_keys_added.txt")


saveGiotto(gobject = fov_integrated_gobject,
           dir = './COSMX_Giotto/from_Anndata', overwrite = TRUE)

fov_integrated <- loadGiotto('./COSMX_Giotto/from_Anndata')



################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# Part I - Single-cell Analysis

# 2- Data exploration and loading

# 2.1- Subcellular detections (points info)

# tx_file.csv contains the subcellular detections information. It contains information on each of the individual feature detections within the sample.
# - fov which FOV the detection happened in
# - cell_ID the ID of the cell the detection happened in
# - x_global_px the global spatial x location in pixels
# - y_global_px the global spatial y location in pixels
# - x_local_px the spatial x location in pixels within the FOV
# - y_local_px the spatial y location in pixels within the FOV
# - z the z plane the detection was called in (-1 to 16)
# - target the feature the probe is targeted against
# - CellComp Cellular compartment the detection happened in (0, Cytoplasm, Membrane, Nuclear)

## provide path to nanostring folder
data_path = './Files/'

# load transcript coordinates
library(data.table)
tx_coord_all = fread(paste0(data_path, 'S1_TMA_tx_file.csv'))

colnames(tx_coord_all)
cat('\n')
# z planes
tx_coord_all[, table(z)]
cat('\n')
# Cell compartment
tx_coord_all[, table(CellComp)]

# 2.2- Split detections by features vs negative probes

# tx_file.csv contains information on both actual features (960 targeted gene probes in this dataset) and negative probes (20) 
# that are targeted to alien sequences defined by the External RNA Controls Consortium (ERCC) that do not exist in human tissue.
# These two types of detections will be treated as separate feature types (feat_type) and placed in separate expression matrices.

all_IDs = tx_coord_all[, unique(target)]
# negative probe IDs
neg_IDs = all_IDs[grepl(pattern = 'NegPrb', all_IDs)]
cat('Negative Probe IDs\n')
neg_IDs
cat('\nFeature IDs\n')
feat_IDs = all_IDs[!all_IDs %in% neg_IDs]
length(feat_IDs)

# split detections
feat_coords_all = tx_coord_all[target %in% feat_IDs]
neg_coords_all = tx_coord_all[target %in% neg_IDs]

cat('\nFeatures: ', feat_coords_all[, .N], '\n',
    'NegProbes: ', neg_coords_all[, .N])

feat_IDs

# 2.2.1- Preview negative probes (optional)

# Previewing the probe information can be done by converting to giottoPoints and then using plot().
# Here we show a preview of the negative probes.
# Note: if previewing the rna expression information, it is highly recommended to set a subset of features using the feats param. 
# The default is to plot all points, which can be very slow for large data.

neg_points = createGiottoPoints(
  x = neg_coords_all[, .(target, x_global_px, y_global_px)]
)
plot(neg_points, point_size = 0.5, feats = neg_IDs) # plot weird if using fovs from different samples!

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 2.3 FOV shifts - IF ANALYSING MORE THAN 1 FOV AT A TIME

# fov_positions_file.csv contains information on the x and y shifts needed in order to put the FOVs tiles together into a cohesive whole. 
# This information is needed during the gobject join process.

#  load field of vision (fov) positions
fov_offset_file = fread(paste0(data_path, 'S1_TMA_fov_positions_file.csv'))

# 2.4 Choose field of view for analysis - IF ANALYSING MORE THAN 1 FOV AT A TIME

# CosMx data is large and Giotto loads in the subcellular information by FOV. 
# This dataset includes 28 FOVs which can be difficult for most computers to handle at once.
# This tutorial will use FOVs ‘02’, ‘03’, and ‘04’ which correspond to the 3 FOVs visible on the bottom right in the negative probe preview above.

gobjects_list = list()

id_set = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25')

# 3. Create a Giotto Object for each FOV - IF ANALYSING MORE THAN 1 FOV AT A TIME

# 3.1- Option 1

image1 = paste0(data_path, 'CellComposite/CellComposite_F001','.jpg')


# cell centroids are now used to provide the spatial locations
fov_integrated_gobject = addSpatialCentroidLocations(fov_integrated_gobject,
                                        poly_info = 'cell')

# create and add Giotto images
fov001_image = createGiottoLargeImage(raster_object = image1,
                                   negative_y = FALSE,
                                   name = 'fov001-image')

fov_integrated_gobject = addGiottoImage(gobject = fov_integrated_gobject,
                           largeImages = list(fov001_image))


fov_integrated_gobject = convertGiottoLargeImageToMG(giottoLargeImage = fov001_image,
                                        #mg_name = 'composite',
                                        gobject = fov_integrated_gobject,
                                        return_gobject = TRUE)

# 3.2- Option 2
for(fov_i in 1:length(id_set)) {
  
  fov_id = id_set[fov_i]
  
  # 1. original composite image as png
  original_composite_image = paste0(data_path, 'CellComposite/CellComposite_F0', fov_id,'.jpg')
  
  # 2. input cell segmentation as mask file
  segmentation_mask = paste0(data_path, 'CellLabels/CellLabels_F0', fov_id, '.tif')
  
  # 3. input features coordinates + offset
  tx_coord = tx_coord_all[fov == as.numeric(fov_id)]
  tx_coord = tx_coord[,.(x_local_px, y_local_px, z, target)]
  colnames(tx_coord) = c('x', 'y', 'z', 'gene_id')
  tx_coord = tx_coord[,.(x, y, gene_id)]
  
  fovsubset = createGiottoObjectSubcellular(
    gpoints = list('rna' = tx_coord),
    gpolygons = list('cell' = segmentation_mask),
    polygon_mask_list_params = list(
      mask_method = 'guess',
      flip_vertical = TRUE,
      flip_horizontal = FALSE,
      shift_horizontal_step = FALSE
    ),
    instructions = instrs
  )
  
  
  # cell centroids are now used to provide the spatial locations
  fovsubset = addSpatialCentroidLocations(fovsubset,
                                          poly_info = 'cell')
  
  # create and add Giotto images
  composite = createGiottoLargeImage(raster_object = original_composite_image,
                                     negative_y = FALSE,
                                     name = 'composite')
  
  fovsubset = addGiottoImage(gobject = fovsubset,
                             largeImages = list(composite))
  
  
  fovsubset = convertGiottoLargeImageToMG(giottoLargeImage = composite,
                                          #mg_name = 'composite',
                                          gobject = fovsubset,
                                          return_gobject = TRUE)
  
  gobjects_list[[fov_i]] = fovsubset
  
  
}

# 4. Join FOV Giotto Objects - IF ANALYSING MORE THAN 1 FOV AT A TIME
new_names = paste0("fov0", id_set)

id_match = match(as.numeric(id_set), fov_offset_file$fov)
x_shifts = fov_offset_file[id_match]$x_global_px
y_shifts = fov_offset_file[id_match]$y_global_px

# Create Giotto object that includes all selected FOVs
fov_join = joinGiottoObjects(gobject_list = gobjects_list,
                             gobject_names = new_names,
                             join_method = 'shift',
                             x_shift = x_shifts,
                             y_shift = y_shifts)

# 5- Visualize Cells and Genes of Interest

# When plotting subcellular data, Giotto uses the spatInSituPlot functions. 
# Spatial plots showing the feature points and polygons are plotted using spatInSituPlotPoints().

showGiottoImageNames(fov_join)

# Set up vector of image names
id_set = c('01')
new_names = paste0("fov0", id_set)
image_names = paste0(new_names, '-image')

spatInSituPlotPoints(fov_join,
                     show_image = TRUE,
                     image_name = image_names,
                     feats = list('rna' = c('MMP2', 'VEGFA', 'IGF1R',
                                            'MKI67', 'EPCAM', "KRT18", "IGF2", "MIF", "CD74")),
                     feats_color_code = pal10,
                     spat_unit = 'cell',
                     point_size = 0.05,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.03,
                     save_param = list(base_height = 3,
                                       save_name = '1_inSituFeats'))

# 5.1- Visualize Cell Centroids

#The standard spatPlot2D() function can also be used, but this works off only the aggregated information that is assembled based on the subcellular information. 
# Plotting information based on cell centroids can be done through this function.

spatPlot2D(gobject = fov_join,
           show_image = TRUE,
           point_shape = 'no_border',
           point_size = 0.01,
           point_alpha = 0.5,
           coord_fix_ratio = 1,
           save_param = list(base_height = 2,
                             save_name = '2_spatCentroids'))


# 6. Aggregate subcellular features

# Giotto supports working directly with the subcellular features in order to generate cell by feature matrices. 
# The data generated this way is then given the spatial unit 'cell'. 
# This workflow is recommended over loading the provided cell by feature (aggregated expression) matrix
# and then including the subcellular information as secondary data.
# When both the raw subcellular information and the pre-made expression matrix are loaded in at the same time, 
# the subcellular data and all data generated from it should be given the spatial unit 'cell' 
# and the pre-generated aggregated information should be given a different spatial unit such as 'cell_agg' to differentiate between the two sources of information.

# In this step, we will be aggregating the feature points of 'rna' and 'neg_probe' into the 'cell' spatial unit.
# Find the feature points overlapped by polygons. This overlap information is then
# returned to the relevant giottoPolygon object's overlaps slot.
fov_join = calculateOverlapRaster(fov_join, feat_info = 'rna')
fov_join = calculateOverlapRaster(fov_join, feat_info = 'neg_probe')

# Convert the overlap information into a cell by feature expression matrix which
# is then stored in the Giotto object's expression slot
fov_join = overlapToMatrix(fov_join, feat_info = 'rna')
fov_join = overlapToMatrix(fov_join, feat_info = 'neg_probe')

showGiottoExpression(fov_join)

# 6.1 Plot histograms of total counts per cell

filterDistributions(fov_join,
                    plot_type = 'hist',
                    detection = 'cells',
                    method = 'sum',
                    feat_type = 'rna',
                    nr_bins = 100,
                    save_param = list(base_height = 3,
                                      save_name = '3.1_totalexpr'))

filterDistributions(fov_join,
                    plot_type = 'hist',
                    detection = 'cells',
                    method = 'sum',
                    feat_type = 'neg_probe',
                    nr_bins = 25,
                    save_param = list(base_height = 3,
                                      save_name = '3.2_totalnegprbe'))
# 6.2- 2D Density Plots

# Density-based representations may sometimes be preferred instead of viewing the raw points information, 
# especially when points are dense enough that there is overplotting. 
# After overlaps information has been calculated, spatInSituPlotDensity() can be used in order to get a general idea of how much expression there is of a feature.
spatInSituPlotDensity(gobject = fov_join,
                      feats = c("MMP2", "VEGFA", "IGF1R",
                                'MKI67', 'EPCAM', 'KRT8'),
                      cow_n_col = 2,
                      save_param = list(base_height = 4,
                                        save_name = '4_inSituDens'))

# 6.3- Extract Data from Giotto Object

# combine cell data
morphometa = combineCellData(fov_join,
                             feat_type = 'rna')

# combine feature data
featmeta = combineFeatureData(fov_join,
                              feat_type = c('rna'))

# combine overlapping feature data
featoverlapmeta = combineFeatureOverlapData(fov_join,
                                            feat_type = c('rna'))

# 6.4- Comparison of Giotto aggregated and Nanostring provided matrices

# Comparison of Giotto’s aggregated matrix results and those provided by Nanostring.
# Only FOV1 will be used in this comparison. 
# Matrices are expected to be similar when the same sets of cell polygons/masks are used for both.
# Load and prepare data

nanoDT = data.table::fread(paste0(data_path, 'Primary_Fov01_exprMat_file.csv')) # remove column before fov or cell_ID if there is in the original file
test1 = nanoDT[fov == 1]
# Set up cell_IDs
test1[, cell_ID := paste0('cell_', cell_ID)]
test1[, cell_ID := paste0('f', fov, '-', cell_ID)]
test1[, fov := NULL]

test1mat = Giotto:::t_flex(Giotto:::dt_to_matrix(test1)) 
testnano_f2 = test1mat
# Remove cell_0 (all tx counts that do not fall within a polygon)
testnano_f2 = testnano_f2[, -1]
# Remove negative probe counts
testnano_f2 = testnano_f2[!grepl('NegPrb', rownames(testnano_f2)),]

# giotto matrix
testg = fov_join@expression$cell$rna$raw[]
testg_f2 = testg[, grepl('', colnames(testg))]
sorted_rownames = sort(rownames(testg_f2))
testg_f2 = testg_f2[sorted_rownames, ]

# Prepare matrix comparison
# Summarise sparse matrices (i and j are matrix indices, x is value)
testg_f2_DT = data.table::as.data.table(Matrix::summary(testg_f2))
testg_f2_DT[, method := 'giotto']

testnano_f2_DT = data.table::as.data.table(Matrix::summary(testnano_f2))
testnano_f2_DT[, method := 'nanostring']

testDT = data.table::rbindlist(list(testg_f2_DT, testnano_f2_DT))
# Combine sparse matrix indices
testDT[, combo := paste0(i,'-',j)]

# Plot results
library(ggplot2)

# matrix index similarity
pl_n = ggplot()
pl_n = pl_n + geom_tile(data = testnano_f2_DT, aes(x = i, y = j, fill = log(x+1)))
pl_n = pl_n + ggtitle('Nanostring Sparse Matrix')
pl_n = pl_n + scale_fill_gradient(low = 'blue', high = 'red')
pl_n = pl_n + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "black"))

pl_g = ggplot()
pl_g = pl_g + geom_tile(data = testg_f2_DT, aes(x = i, y = j, fill = log(x+1)))
pl_g = pl_g + ggtitle('Giotto Sparse Matrix')
pl_g = pl_g + scale_fill_gradient(low = 'blue', high = 'red')
pl_g = pl_g + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "black"))


combplot = cowplot::plot_grid(pl_n, pl_g,
                              nrow = 2,
                              labels = 'AUTO')
print(combplot)

# directly compare differences in matrix values (counts assigned)
vartestDT = testDT[, list(var = var(x), diff = diff(x), mean = mean(x)), by = .(i,j)]
data.table::setorder(vartestDT, var)

# check arbitrary index values
testDT[i == '812' & j == '2']
testDT[i == '667' & j == '1072']
testDT[i == '667' & j == '2880']

# plot difference in values
pl = ggplot()
pl = pl + geom_bar(data = vartestDT, aes(x = diff))
pl = pl + theme_bw()
pl = pl + labs(x = 'difference nanostring - Giotto')
pl

testDT[order(x)]

testDT[, .N, by = 'method']

testDT[, method, by = combo][, sum(duplicated(combo))]

# Overall, the nanostring matrix has 207672 - 207496 = 176 more non-zero values than giotto’s matrix for FOV1. 
# Within the 204397 shared entries that were called by both methods (common i and j indices), there appears to be no major bias in terms of counts/values assigned. 
# Moreover, the vast majority of these shared entries have the same values (difference of 0).

# 7. Filtering and normalization

# After the expression matrix is generated from the subcellular information, analysis proceeds through data filtering and normalization.
# For the normalization step, we will employ two types.

# A- standard normalization method: library size normalization and log normalization. 
# This method will produce both normalized and scaled values that are be returned as the ‘normalized’ and ‘scaled’ expression matrices respectively. 
# In this tutorial, the normalized values will be used for generating expression statistics and plotting expression values. 
# The scaled values will be ignored. We will also generate normalized values for the negative probes for visualization purposes during which the library normalization step will be skipped.

# B- pearson residuals: A normalization that uses the method described in Lause/Kobak et al. 2021. 
# This produces a set of values that are most similar in utility to a scaled matrix and offer improvements to both HVF detection and PCA generation. 
# These values should not be used for statistics, plotting of expression values, or differential expression analysis.

# filter (feat_type = 'rna' by default)
fov_join <- filterGiotto(gobject = fov_join,
                         feat_type = 'rna',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5)

# normalize
# standard method of normalization (log normalization based)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            norm_methods = 'standard',
                            verbose = TRUE)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'neg_probe',
                            norm_methods = 'standard',
                            library_size_norm = FALSE,
                            verbose = TRUE)

# new normalization method based on pearson correlations (Lause/Kobak et al. 2021)
# this normalized matrix is given the name 'pearson' using the update_slot param
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            scalefactor = 5000,
                            verbose = TRUE,
                            norm_methods = 'pearson_resid',
                            update_slot = 'pearson')

showGiottoExpression(fov_join)

# add statistics based on log normalized values for features rna and negative probes
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'rna')
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'neg_probe')

# View cellular data (default is feat = 'rna')
showGiottoCellMetadata(fov_join)
# View feature data
showGiottoFeatMetadata(fov_join)

# 8. View Transcript Total Expression Distribution

# 8.1 Histogram of log normalized data
filterDistributions(fov_join,
                    detection = 'cells',
                    feat_type = 'rna',
                    expression_values = 'normalized',
                    method = 'sum',
                    nr_bins = 100,
                    save_param = list(base_height = 3,
                                      save_name = '5.1_rna_norm_total_hist'))

filterDistributions(fov_join,
                    detection = 'cell',
                    feat_type = 'neg_probe',
                    expression_values = 'normalized',
                    method = 'sum',
                    nr_bins = 20,
                    save_param = list(base_height = 3,
                                      save_name = '5.2_neg_norm_total_hist'))

# 8.2 Plot spatially as centroids
spatPlot2D(gobject = fov_join,
           cell_color = 'total_expr',
           color_as_factor = FALSE,
           show_image = FALSE,
           image_name = "image",
           point_size = 0.9,
           point_alpha = 0.75,
           save_param = list(base_height = 2,
                             save_name = '5.3_color_centroids'))

# 8.3 Plot spatially as color-scaled polygons
spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     save_param = list(base_height = 2,
                                       save_name = '5.4_rna_color_polys'))
spatInSituPlotPoints(fov_join,
                     feat_type = 'neg_probe',
                     show_polygon = TRUE,
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     save_param = list(base_height = 2,
                                       save_name = '5.5_neg_color_polys'))

# 9. Dimension Reduction

# 9.1 Detect highly variable genes and generate PCA

# Detect highly variable genes using the pearson residuals method based on the ‘pearson’ expression matrix. T
# These results will be returned as a new ‘hvf’ column in the ‘rna’ feature metadata.

fov_join = calculateHVF(fov_join,
                        method = 'var_p_resid',
                        expression_values = 'pearson',
                        save_param = list(base_height = 5,
                                          save_name = '6.1_pearson_HVF'))

# print HVFs
gene_meta = fDataDT(fov_join)
gene_meta[hvf == 'yes', feat_ID]

# PCA generation will also be based on the ‘pearson’ matrix. 
# Scaling and centering of the PCA which is usually done by default will be skipped since the pearson matrix is already scaled.

fov_join = runPCA(fov_join,
                  scale_unit = FALSE,
                  center = FALSE,
                  expression_values = 'pearson')

# screeplot uses the generated PCA. No need to specify expr values
screePlot(fov_join, ncp = 30, save_param = list(save_name = '6.2_screeplot'))

plotPCA(fov_join,
        cell_color = 'nr_feats', # (from log norm statistics)
        color_as_factor = FALSE,
        point_size = 0.1,
        point_shape = 'no_border',
        save_param = list(save_name = '6.3_PCA'))

# 9.2 Run UMAP
# Generate UMAP from PCA
fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:30,
                    n_threads = 4)

plotUMAP(gobject = fov_join, save_param = list(save_name = '6.4_UMAP'))

# 9.3 Plot features on expression space
dimFeatPlot2D(gobject = fov_join,
              feat_type = 'rna',
              feats = c('MKI67', 'CD8A', 'CD4',
                        'COL1A1', 'MS4A1', 'MZB1'),
              expression_values = 'normalized',
              point_shape = 'no_border',
              point_size = 0.01,
              cow_n_col = 3,
              save_param = list(base_height = 5,
                                save_name = '6.5_UMAP_feats'))

# 10- Cluster
# 10.1- Visualize clustering
fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:30,
                                 k = 30)

fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 1,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         cell_color_code = pal10,
         show_NN_network = FALSE,
         point_size = 2,
         save_param = list(save_name = '7.1_UMAP_leiden'))

plotUMAP_3D(gobject = fov_join, cell_color = 'leiden_clus',
            cell_color_code = pal10,
            show_NN_network = FALSE,
            point_size = 2,
            save_param = list(save_name = '7.1.1_UMAP_leiden'))

# 10.2- Visualize clustering on expression and spatial space

# visualize UMAP and spatial results
spatDimPlot2D(gobject = fov_join,
              show_image = TRUE,
              image_name = image_names,
              cell_color = 'leiden_clus',
              cell_color_code = pal10,
              spat_point_size = 1,
              save_param = list(save_name = '7.2_spatdim_leiden'))

# 10.3 Map clustering spatially with RNA info on top or not 
spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.05,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     save_param = list(base_height = 5,
                                       save_name = '7.3_spatinsitu_leiden'))


spatInSituPlotPoints(fov_join,
                     feats = list('rna' = c('MMP2', 'VEGFA', 'IGF1R',
                                            'KRT18', 'EPCAM', 'CD74')),
                     point_size = 0.3,
                     feats_color_code = viv10,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.05,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     save_param = list(base_height = 5,
                                       save_name = '7.3_spatinsitu_leiden'))

# 11. Subcellular Visualization
#subset a Giotto object based on spatial locations
smallfov <- subsetGiottoLocs(fov_join,
                             x_max = 4000,
                             x_min = 2000,
                             y_max = 4000,
                             y_min = 2000)

#extract all genes observed in new object
smallfeats <- fDataDT(smallfov)[, feat_ID]

#plot all genes
spatInSituPlotPoints(smallfov,
                     feats = list(smallfeats),
                     point_size = 0.15,
                     polygon_line_size = 0.1,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     show_image = TRUE,
                     largeImage_name = 'composite',
                     show_legend = FALSE,
                     save_param = list(save_name = '8.1_smallfov_points'))

# plot only the polygon outlines
spatInSituPlotPoints(smallfov,
                     polygon_line_size = 0.1,
                     polygon_alpha = 0,
                     polygon_color = 'white',
                     show_polygon = TRUE,
                     show_image = TRUE,
                     largeImage_name = 'composite',
                     show_legend = FALSE,
                     save_param = list(save_name = '8.2_smallfov_poly'))

# plot polygons color labeled with leiden clusters
spatInSituPlotPoints(smallfov,
                     polygon_line_size = 0.1,
                     show_polygon = TRUE,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     show_image = TRUE,
                     largeImage_name = 'composite',
                     show_legend = FALSE,
                     save_param = list(save_name = '8.3_smallfov_leiden'))


# 12. Identify cluster differential expression genes

# 12.1 Violin plot
# Scran
markers = findMarkers_one_vs_all(gobject = fov_join,
                                 method = 'scran',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_clus',
                                 logFC = 0.25)
# First 10 results by cluster
markers[, head(.SD, 10), by = 'cluster']
write.csv(markers, "./Plots/2_Clustering_Annotation/Giotto_top_DGE_leiden_clusters_fov01_scran.csv")

markers <- read.csv("./Plots/2_Clustering_Annotation/Giotto_top_DGE_leiden_clusters_fov01_scran.csv")

# violinplot
topgini_genes = unique(markers[, head(.SD, 10), by = 'cluster']$feats)
violinPlot(fov_join,
           feats = topgini_genes,
           cluster_column = 'leiden_clus',
           strip_position = 'right',
           save_param = list(save_name = '10.1_gini_violin'))

# Heatmap
cluster_order = 1:11
plotMetaDataHeatmap(fov_join,
                    expression_values = 'normalized',
                    metadata_cols = c('leiden_clus'),
                    selected_feats = topgini_genes,
                    custom_cluster_order = cluster_order,  y_text_size = 5,
                    save_param = list(base_height = 5,
                                      save_name = '10.2_heatmap_scran'))

# 13- Annotate Giotto Object - see methods of integration with scRNAseq datasets

plotMetaDataHeatmap(fov_join,
                    expression_values = 'normalized',
                    metadata_cols = c('leiden_clus'),
                    selected_feats = c("CD8A", "CD8B"),
                    custom_cluster_order = cluster_order,  y_text_size = 5,
                    save_param = list(base_height = 5,
                                      save_name = 'heatmap_scran_CD8'))

plotMetaDataHeatmap(fov_join,
                    expression_values = 'normalized',
                    metadata_cols = c('leiden_clus'),
                    selected_feats = c("CD3E", "CD4"),
                    custom_cluster_order = cluster_order,  y_text_size = 5,
                    save_param = list(base_height = 5,
                                      save_name = 'heatmap_scran_CD3_CD4'))

# Match features from cosmx to the reference marker gene list
markers_ref <- read.csv("./E-MTAB-8410.marker_genes_inferred_cell_type_-_ontology_labels.csv")

ortho <- markers_ref
cluster <- markers
idx <- match(ortho$geneID, cluster$feats)
ortho$feats <- cluster$feats[ idx ]
ortho <- na.omit(ortho) 
write.csv(ortho, "marker_genes_reference_cosmx_matched.csv")

# add geneID from the matched reference marker gene list above in the cosmx top marker genes list
markers_ref <- read.csv("./marker_genes_reference_cosmx_matched.csv")
ortho <- markers_ref
cluster <- markers
idx <- match(ortho$geneID, cluster$feats)
cluster$geneID <- ortho$geneID[ idx ]
#ortho <- na.omit(ortho) 
write.csv(ortho, "marker_genes_reference_cosmx_matched.csv")

# UMAP
# low, mid, high
custom_scale = c('#440154', '#1F968B', '#FDE725')

dimFeatPlot2D(fov_join,
              expression_values = 'normalized',
              cell_color_gradient = custom_scale,
              gradient_midpoint = 5,
              feats = topgini_genes,
              point_shape = 'no_border',
              point_size = 0.001,
              cow_n_col = 4,
              save_param = list(base_height = 8,
                                save_name = '10.3_gini_genes'))

## add cell types ###
clusters_cell_types = c('NK_CD8 T cell', 'Monocte_Macrophage', 'Epithelial cell', 'Tuft cell',
                             'CD4 T cell', 'Neutrophils', 'Plasma cell', 'CD16+ Monocyte', 'CCL2+SPP1+ Monocyte', 'B cell', 'CXCL10+ Epithelial cell')

names(clusters_cell_types) = 1:11
fov_join = annotateGiotto(gobject = fov_join,
                          annotation_vector = clusters_cell_types,
                          cluster_column = 'leiden_clus',
                          name = 'cell_types')

# Visualise top expressed genes in annotated clusters 
plotMetaDataHeatmap(fov_join,
                    expression_values = 'normalized',
                    metadata_cols = c('cell_types'),
                    selected_feats = topgini_genes,
                    y_text_size = 5,
                    save_param = list(base_height = 5,
                                      save_name = 'heatmap_scran_leiden_annotated'))

# Visualise annotated clusters on UMAP
plotUMAP(fov_join,
         cell_color = 'cell_types',
         cell_color_code = pal10,
         point_size = 1.5,
         save_param = list(save_name = '11_anno_umap'))

# Visualise annotated clusters on UMAP and spatially
spatDimPlot2D(gobject = fov_join,
              show_image = FALSE,
              image_name = image,
              cell_color = 'cell_types',
              cell_color_code = pal10,
              spat_point_size = 1.5, plot_alignment = 'horizontal',
              save_param = list(save_name = '12_spatdim_type'))

# show only a subset of the clusters
spatDimPlot(fov_join,
            cell_color = 'cell_types',  select_cell_groups = c('Epithelial cell', 'Tuft cell'),
            plot_alignment = 'horizontal', spat_point_size = 2)

spatDimPlot(fov_join,
            cell_color = 'cell_types',  select_cell_groups = c('CXCL10+ Epithelial cell'),
            plot_alignment = 'horizontal', spat_point_size = 2)

# show only a subset of the clusters, excluding not selected cells
spatDimPlot(fov_join,
            cell_color = 'cell_types',  select_cell_groups = 'Epithelial cell', show_other_cells = F,
            plot_alignment = 'horizontal', spat_point_size = 2)

# create spatial plots with all leiden clustering result
spatPlot(gobject = fov_join,
         show_image = FALSE,
         image_name = image,
         cell_color = 'cell_types',
         cell_color_code = pal10,
         point_size = 2,
         save_param = list(save_name = '12_spatdim_type'))

# create spatial plots grouped by the leiden clustering result
spatPlot(gobject = fov_join,
         group_by = 'cell_types',
         cell_color = 'cell_types', show_other_cells = FALSE,
         point_size = 1, cow_n_col = 4, axis_text = 6, axis_title = 6, legend_text = 4,
         save_param = list(save_name = '12_spatdim_type'))

# Visualise annotated clusters on top of the segmentation mask

spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'grey',
                     polygon_line_size = 0.05,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     save_param = list(base_height = 2,
                                       save_name = '13_insitu_type'))

# Visualise annotated clusters on top of the segmentation mask plus specific transcripts

c('ICOS', 'EPCAM', 'CD24', 'CD33', 'JCHAIN', 'CXCL8',
  'MIF', 'CD4', 'FCGR3A', 'CCL2', 'CSF1','CXCL10')

spatInSituPlotPoints(fov_join,
                     feats = list('rna' = c('ICOS', 'CD33', 'JCHAIN', 'CXCL8',
                                            'CD4', 'FCGR3A', 'CCL2', 'CSF1')),
                     point_size = 0.5,
                     feats_color_code = viv10,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.05,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     save_param = list(base_height = 5,
                                       save_name = 'spatinsitu_leiden_annotated_immune_cell_transcripts'))

spatInSituPlotPoints(fov_join,
                     feats = list('rna' = c('EPCAM', 'CD24', 'MIF', 'CXCL10')),
                     point_size = 0.5,
                     feats_color_code = viv10,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.05,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     save_param = list(base_height = 5,
                                       save_name = 'spatinsitu_leiden_annotated_epithelial_transcripts'))


######################################################################################################################################################################################################

# Part II - Spatial Analysis

# Gene-level analysis

# 14- Build Spatial Graph - Spatial network

# visualize information about the default Delaunay network
# create a spatial Delaunay network (default)
# create a spatial kNN network

# create spatial network based on delaunay triangulation and physical distance of cell centroids
COVID_gobject = createSpatialNetwork(gobject = COVID_gobject,
                                minimum_k = 2, method = "Delaunay", delaunay_method = c("delaunayn_geometry"),
                                maximum_distance_delaunay = 20)

# create spatial network based on KNN and physical distance of cell centroids
COVID_gobject = createSpatialNetwork(gobject = COVID_gobject, minimum_k = 2,
                                    method = 'kNN', k = 4, maximum_distance_knn = 20)
showGiottoSpatNetworks(COVID_gobject)

HIVNEG_gobject = createSpatialNetwork(gobject = HIVNEG_gobject, minimum_k = 2,
                                     method = 'kNN', k = 4, maximum_distance_knn = 20)
showGiottoSpatNetworks(HIVNEG_gobject)

HIVPOS_gobject = createSpatialNetwork(gobject = HIVPOS_gobject, minimum_k = 2,
                                     method = 'kNN', k = 4, maximum_distance_knn = 20)
showGiottoSpatNetworks(HIVPOS_gobject)

# visualize the two different spatial networks
spatPlot(gobject = COVID_gobject, show_network = T, group_by = "ROI", group_by_subset = "0305_1_A",
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 1.5, cell_color = 'pheno_cluster', save_param = list(save_name = 'DelaunayNetwork'))

spatPlot(gobject = COVID_gobject, show_network = T, group_by = "ROI", group_by_subset = "0305_1_A",
         network_color = 'blue', spatial_network_name = 'kNN_network',
         point_size = 1.5, cell_color = 'pheno_cluster', save_param = list(save_name = 'kNN_network'))

## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = COVID_gobject, spatial_network_name = "kNN_network",
                                           cluster_column = 'pheno_cluster', adjust_method = "BH",
                                           number_of_simulations = 1000, set_seed = TRUE,
                                           seed_number = 1234)
## barplot
cellProximityBarplot(gobject = COVID_gobject, CPscore = cell_proximities, min_orig_ints = 25, min_sim_ints = 25,
                     save_param = c(save_name = 'COVID_barplot_cell_cell_enrichment')) + theme(axis.title.y = element_text(size = 2))
## heatmap
cellProximityHeatmap(gobject = COVID_gobject, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'), show_plot = T, return_plot = T,
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))
## network
cellProximityNetwork(gobject = COVID_gobject, CPscore = cell_proximities, remove_self_edges = F, only_show_enrichment_edges = T,
                     save_param = c(save_name = '12_c_network_cell_cell_enrichment')) + theme(text = element_text(size = 2))


## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = HIVNEG_gobject, spatial_network_name = "kNN_network",
                                           cluster_column = 'pheno_cluster', adjust_method = "BH",
                                           number_of_simulations = 1000, set_seed = TRUE,
                                           seed_number = 1234)
## barplot
cellProximityBarplot(gobject = HIVNEG_gobject, CPscore = cell_proximities, min_orig_ints = 25, min_sim_ints = 25,
                     save_param = c(save_name = 'COVID_barplot_cell_cell_enrichment')) + theme(axis.title.y = element_text(size = 2))
## heatmap
cellProximityHeatmap(gobject = HIVNEG_gobject, CPscore = cell_proximities, order_cell_types = T, scale = T, return_plot = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))
## network
cellProximityNetwork(gobject = HIVNEG_gobject, CPscore = cell_proximities, remove_self_edges = F, only_show_enrichment_edges = T,
                     save_param = c(save_name = '12_c_network_cell_cell_enrichment')) + theme(text = element_text(size = 2))


## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = HIVPOS_gobject, spatial_network_name = "kNN_network",
                                           cluster_column = 'pheno_cluster', adjust_method = "BH",
                                           number_of_simulations = 1000, set_seed = TRUE,
                                           seed_number = 1234)
## barplot
cellProximityBarplot(gobject = HIVPOS_gobject, CPscore = cell_proximities, min_orig_ints = 25, min_sim_ints = 25,
                     save_param = c(save_name = 'COVID_barplot_cell_cell_enrichment')) + theme(axis.title.y = element_text(size = 2))
## heatmap
cellProximityHeatmap(gobject = HIVPOS_gobject, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in')) + theme(text = element_text(size = 0.5))
## network
cellProximityNetwork(gobject = HIVPOS_gobject, CPscore = cell_proximities, remove_self_edges = F, only_show_enrichment_edges = T,
                     save_param = c(save_name = '12_c_network_cell_cell_enrichment')) + theme(text = element_text(size = 2))




# 15. Spatial Expression Patterns

# Find spatially organized gene expression by examining the binarized expression of cells and their spatial neighbors.
# perform Binary Spatial Extraction of genes - NOTE: Depending on your system this could take time
km_spatialgenes = binSpect(fov_join)
write.csv(km_spatialgenes, "km_spatialgenes.csv")

# visualize spatial expression of selected genes obtained from binSpect
spatFeatPlot2D(fov_join,
               expression_values = 'normalized',
               feats = km_spatialgenes$feats[1:10],
               point_shape = 'no_border',
               point_border_stroke = 0.01,
               point_size = 0.01,
               cow_n_col = 2,
               save_param = list(save_name = '9_binspect_genes'))

rank_spatialgenes = binSpect(fov_join, bin_method = 'rank')
write.csv(rank_spatialgenes, "rank_spatialgenes.csv")

spatFeatPlot2D(fov_join, expression_values = 'scaled',
               feats = rank_spatialgenes[1:6]$feats,
               point_shape = 'border', point_border_stroke = 0.1,
               show_network = T, network_color = 'lightgrey', point_size = 1.5,
               cow_n_col = 2, save_param = list(save_name = '15_Top_spatial_expression_genes'))


# 16- Spatial co-expression patterns

# Identify robust spatial co-expression patterns using the spatial network or grid and a subset of individual spatial genes.
# calculate spatial correlation scores
# cluster correlation scores

# 16.1- calculate spatial correlation scores for the top 500 spatially variable genes
ext_spatial_genes = km_spatialgenes[1:500]$feats

spat_cor_netw_DT = detectSpatialCorFeats(fov_join,
                                         method = 'network',
                                         spatial_network_name = 'kNN_network',
                                         subset_feats = ext_spatial_genes) #it can use either spatial graphs
# 16.2- cluster correlation scores
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT,
                                          name = 'spat_netw_clus', k = 8)

# visualise correlated genes
heatmSpatialCorFeats(fov_join, spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus')

# rank spatial correlated clusters and show genes for selected clusters
netw_ranks = rankSpatialCorGroups(fov_join,
                                  spatCorObject = spat_cor_netw_DT,
                                  use_clus_name = 'spat_netw_clus')

# See top 10 genes coexpressed per each gene in specific clusters
top_netw_spat_cluster8 = showSpatialCorFeats(spat_cor_netw_DT,
                                            use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 8,
                                            show_top_feats = 10)

# list of correlated genes and enrichment score for clusters
cluster_genes_DT = showSpatialCorFeats(spat_cor_netw_DT,
                                       use_clus_name = 'spat_netw_clus',
                                       show_top_feats = 10)

# create metagene enrichment score for clusters
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$feat_ID

fov_join = createMetafeats(fov_join,
                               feat_clusters = cluster_genes,
                               name = 'cluster_metagene')

# Plot clusters of co-expressed genes in the spatial image
spatCellPlot(fov_join,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, cow_n_col = 3)

# Identify and visualise most similar spatially correlated genes for one gene
EPCAM_top10_genes = showSpatialCorFeats(spat_cor_netw_DT, feats = 'ABL1', show_top_feats = 10)

spatFeatPlot2D(fov_join, expression_values = 'scaled',
               feats = c('ABL1', 'EPCAM'), point_size = 1)


######################################################################################################################################################################################################


# Cell-level analysis

# 17- Spatial HMRF Domains

#remotes::install_bitbucket(repo = 'qzhudfci/smfishhmrf-r', ref='master')
library(smfishHmrf)

hmrf_folder = paste0(temp_dir,'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
# perform hmrf
km_spatialgenes <- read.csv("./km_spatialgenes.csv")
my_spatial_genes = km_spatialgenes$feats
my_spatial_genes <- my_spatial_genes[1:500]

HMRF_spatial_genes = doHMRF(gobject = fov_join,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            spatial_network_name = 'Delaunay_network',
                            k = 4,
                            betas = c(28,2,2),
                            output_folder = "./Plots/5_HMRF")

# check and select hmrf
for(i in seq(28, 30, by = 2)) {
  viewHMRFresults2D(gobject = fov_join,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 4, betas_to_view = i,
                    point_size = 2)
}

fov_join = addHMRF(gobject = fov_join,
                       HMRFoutput = HMRF_spatial_genes,
                       k = 4, betas_to_add = c(30),
                       hmrf_name = 'HMRF')
saveRDS(HMRF_spatial_genes, "HMRF_spatial_genes.RDS")

# visualize selected hmrf result
giotto_colors = Giotto:::getDistinctColors(9)
names(giotto_colors) = 1:4
spatPlot(gobject = fov_join, cell_color = 'HMRF_k4_b.30',
         point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)
   

# 18- Cell neighborhood: cell-type/cell-type interactions
set.seed(seed = 2841)
cell_proximities = cellProximityEnrichment(gobject = fov_join,
                                           cluster_column = 'cell_types',
                                           spatial_network_name = 'kNN_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)
# barplot
cellProximityBarplot(gobject = fov_join,
                     CPscore = cell_proximities,
                     min_orig_ints = 5, min_sim_ints = 5, p_val = 0.5)

## heatmap
cellProximityHeatmap(gobject = fov_join, CPscore = cell_proximities,
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5),
                     color_names = c('blue', 'white', 'red'))

# network
cellProximityNetwork(gobject = fov_join, CPscore = cell_proximities,
                     remove_self_edges = T, only_show_enrichment_edges = T)


# network with self-edges
cellProximityNetwork(gobject = fov_join, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5))


# Visualisation of specific cell types - Option 1
spec_interaction = "CXCL10+ Epithelial cell--CD16+ Monocyte"
cellProximitySpatPlot2D(gobject = fov_join,
                        interaction_name = spec_interaction,
                        show_network = T,
                        cluster_column = 'cell_types',
                        cell_color = 'cell_types',
                        cell_color_code = c('CXCL10+ Epithelial cell' = 'lightblue', 'CD16+ Monocyte' = 'red'),
                        point_size_select = 4, point_size_other = 2)


# Option 2: create additional metadata
fov_join = addCellIntMetadata(fov_join,
                                  spat_unit = "cell",
                                  spatial_network = 'Delaunay_network',
                                  cluster_column = 'cell_types',
                                  cell_interaction = spec_interaction,
                                  name = 'Epithelial_Monocyte_Interactions')
spatPlot(fov_join, cell_color = 'Epithelial_Monocyte_Interactions', legend_symbol_size = 3, show_network = T, network_color = 'black',
         select_cell_groups =  c('other_CXCL10+ Epithelial cell', 'other_CD16+ Monocyte', 'select_CXCL10+ Epithelial cell', 'select_CD16+ Monocyte'))

############################################################################################################################################################################################################################################################

# 19- Interaction Changed Genes - after cell annotation of leiden clusters
future::plan('multisession', workers = 4) # NOTE: Depending on your system this could take time


# identify features (genes) expressed in cell A that are associated with proximity to other cell types
# to run this,m the data has to be normalised in giotto first. I used the normalisation code from the codex tutorial

COVID_gobject <- normalizeGiotto(gobject = COVID_gobject,
                              scalefactor = 6000,
                              verbose = T,
                              log_norm = FALSE,
                              library_size_norm = FALSE,
                              scale_feats = FALSE,
                              scale_cells = TRUE)

HIVNEG_gobject <- normalizeGiotto(gobject = HIVNEG_gobject,
                                 scalefactor = 6000,
                                 verbose = T,
                                 log_norm = FALSE,
                                 library_size_norm = FALSE,
                                 scale_feats = FALSE,
                                 scale_cells = TRUE)

HIVPOS_gobject <- normalizeGiotto(gobject = HIVPOS_gobject,
                                 scalefactor = 6000,
                                 verbose = T,
                                 log_norm = FALSE,
                                 library_size_norm = FALSE,
                                 scale_feats = FALSE,
                                 scale_cells = TRUE)


icf = findInteractionChangedFeats(gobject = COVID_gobject, spatial_network_name = 'kNN_network',
                                  cluster_column = 'pheno_cluster', nr_permutations = 1000,
                                  do_parallel = T,
                                  set_seed = TRUE,
                                  seed_number = 1234)

icf2 = findInteractionChangedFeats(gobject = HIVNEG_gobject, spatial_network_name = 'kNN_network',
                                  cluster_column = 'pheno_cluster', nr_permutations = 1000,
                                  do_parallel = T,
                                  set_seed = TRUE,
                                  seed_number = 1234)

icf3 = findInteractionChangedFeats(gobject = HIVPOS_gobject, spatial_network_name = 'kNN_network',
                                  cluster_column = 'pheno_cluster', nr_permutations = 1000,
                                  do_parallel = T,
                                  set_seed = TRUE,
                                  seed_number = 1234)

#diff_test = 'permutation', default
#adjust_method = 'bonferroni', default
saveRDS(icf, " Interaction_changed_features_COVID.rds")
saveRDS(icf2, " Interaction_changed_features_HIVNEG.rds")
saveRDS(icf3, " Interaction_changed_features_HIVPOS.rds")

icf$ICFscores[type_int == 'hetero']
icf2$ICFscores[type_int == 'hetero']
icf3$ICFscores[type_int == 'hetero']

## visualize all genes
plotCellProximityFeats(COVID_gobject, icfObject = icf2, method = 'dotplot')
plotCellProximityFeats(COVID_gobject, icfObject = icf, method = 'volcano', cell_color_code = viv10)
plotCellProximityFeats(COVID_gobject, icfObject = icf, method = 'cell-cell', cell_color_code = viv10, direction = c("both"))

## visualize subset of interaction changed genes (ICGs)
# here, select the cells to test which proteins expression changes when interacting with the source cell.
# run for CD4+, CD8+ and Tregs
ICF_genes = c('MHCII', 'MHCI', 'CD107a', 'GranzymeB','Ki67', 'Arginase1', 'CD38', 'CD45RO')
ICF_genes_types = c('SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil')

names(ICF_genes) = ICF_genes_types

plotICF(gobject = HIVPOS_gobject,
        icfObject = icf3,
        source_type = "CD8 T cell", #change to the source cell type, here would be SARSCoV2+ neutrophils
        source_markers = c("CD3", "CD8"), cell_color_code = "dark blue",
        ICF_feats = ICF_genes)

ICF_genes = c('MHCII', 'MHCI', 'CD107a', 'GranzymeB','Ki67', 'Arginase1', 'CD38', 'CD31')
"CD4 Treg cell", 'Proliferative CD4 T cell','EM CD4 T cell', "CD8 T cell",'Proliferative CD8 T cell', "EM CD8 T cell"

# Identify top 20 interaction changed features
icf$ICFscores[type_int == 'hetero']$feats[1:20]

icf_plotfeats = icf$ICFscores[type_int == 'hetero']$feats[1:10]

# Visualize ICF expression
spatInSituPlotPoints(COVID_gobject,
                     feats = list(icf_plotfeats),
                     point_size = 0.1,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal10,
                     save_param = list(base_height = 6,
                                       save_name = '14_ICF'))


############################################################################################################################################################################################################################################################

# 20- Ligand-receptor cell-cell communication

LR_data = data.table::fread("./human_lr_pair.txt")
LR_data[, ligand_det := ifelse(ligand_gene_symbol %in% fov_join@feat_ID[['rna']], T, F)]
LR_data[, receptor_det := ifelse(receptor_gene_symbol %in% fov_join@feat_ID[['rna']], T, F)]

LR_data_det = LR_data[ligand_det == T & receptor_det == T]
write.csv(LR_data_det, "Ligand-receptor_pairs_CosMX_fov01.csv")

select_ligands = LR_data_det$ligand_gene_symbol
select_receptors = LR_data_det$receptor_gene_symbol


## get statistical significance of gene pair expression changes based on expression ##
expr_only_scores = exprCellCellcom(gobject = fov_join,
                                   cluster_column = 'cell_types',
                                   random_iter = 50,
                                   feat_set_1 = select_ligands,
                                   feat_set_2 = select_receptors)

## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores = spatCellCellcom(fov_join,
                                     spat_unit = 'cell',
                                     feat_type = 'rna',
                                     spatial_network_name = 'Delaunay_network',
                                     cluster_column = 'cell_types',
                                     random_iter = 50,
                                     feat_set_1 = select_ligands,
                                     feat_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = 4,
                                     verbose = 'none')

## * plot communication scores ####
## select top LR ##
selected_spat = spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
data.table::setorder(selected_spat, -PI)
top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:33]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]
plotCCcomHeatmap(gobject = seqfish_mini,
                 comScores = spatial_all_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 show = 'LR_expr')
