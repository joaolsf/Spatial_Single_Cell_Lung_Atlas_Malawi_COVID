


# Ensure Giotto Suite is installed.
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}

library(Giotto)
# Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}

# set working directory
results_folder = "./Plots/"

# Package Accessibility
default_instrs <- createGiottoInstructions()

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

# 5- Creating Giotto object from Anndata

# Reference: https://giottosuite.readthedocs.io/en/latest/subsections/datasets/interoperability_04122023.html
# To convert an AnnData Object back into a Giotto object, it must first be saved as a .h5ad file. 
# The name of said file may then be provided to anndataToGiotto() for conversion.
# If a nearest neighbor network or spatial netowkr was created using the key_added argument, they may be provided to arguments n_key_added and/or spatial_n_key_added, respectively.

COVID_gobject <- anndataToGiotto(anndata_path = "./adata_COVID.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE)

saveGiotto(gobject = fov_integrated_gobject,
           dir = './COSMX_Giotto/from_Anndata', overwrite = TRUE)

fov_integrated <- loadGiotto('./COSMX_Giotto/from_Anndata')

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 5- Creating Giotto object from csv files

# Reference: https://giottosuite.readthedocs.io/en/latest/subsections/datasets/mouse_CODEX_spleen.html
# Export the tables from anndata.
# obs.csv = annotation file
# create coord file
# X.csv = expression file - I think it has to be transposed with cells as columns and proteins as rows

# 5.1- create giotto object from provided paths ####
expr_path = paste0(results_folder, "covid_expression.csv")
loc_path = paste0(results_folder, "covid_coord.csv")
meta_path = paste0(results_folder, "covid_annotation.csv")

expr_path = paste0(results_folder, "hivneg_expression.csv")
loc_path = paste0(results_folder, "hivneg_coord.csv")
meta_path = paste0(results_folder, "hivneg_annotation.csv")

expr_path = paste0(results_folder, "hivpos_expression.csv")
loc_path = paste0(results_folder, "hivpos_coord.csv")
meta_path = paste0(results_folder, "hivpos_annotation.csv")

# 5.2 read in data information

# expression info
codex_expression = readExprMatrix(expr_path, transpose = T)
# cell coordinate info
codex_locations = data.table::fread(loc_path)
# metadata
codex_metadata = data.table::fread(meta_path)

# create Giotto object
COVID_gobject <- createGiottoObject(expression = codex_expression,
                                 spatial_locs = codex_locations,
                                 instructions = instrs)

codex_metadata$cell_ID<- as.character(codex_metadata$cellID)
COVID_gobject <-addCellMetadata(COVID_gobject, new_metadata = codex_metadata,
                            by_column = T,
                            column_cell_ID = "cell_ID")


saveGiotto(gobject = COVID_gobject,
           dir = './COVID_Giotto_object', overwrite = TRUE)

HIVNEG_gobject <- createGiottoObject(expression = codex_expression,
                                    spatial_locs = codex_locations,
                                    instructions = instrs)

codex_metadata$cell_ID<- as.character(codex_metadata$cellID)
HIVNEG_gobject <-addCellMetadata(HIVNEG_gobject, new_metadata = codex_metadata,
                                by_column = T,
                                column_cell_ID = "cell_ID")


saveGiotto(gobject = HIVNEG_gobject,
           dir = './HIVNEG_Giotto_object', overwrite = TRUE)

HIVPOS_gobject <- createGiottoObject(expression = codex_expression,
                                    spatial_locs = codex_locations,
                                    instructions = instrs)

codex_metadata$cell_ID<- as.character(codex_metadata$cellID)
HIVPOS_gobject <-addCellMetadata(HIVPOS_gobject, new_metadata = codex_metadata,
                                by_column = T,
                                column_cell_ID = "cell_ID")


saveGiotto(gobject = HIVPOS_gobject,
           dir = './HIVPOS_Giotto_object', overwrite = TRUE)

