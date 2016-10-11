# LandformClassification
ArcGIS Python Toolbox containing tools for classifying landforms in low-relief agricultural fields as described by Evans et al., 2016.  This article can be downloaded at http://scholarworks.sfasu.edu/j_of_geospatial_applications_in_natural_resources/vol1/iss1/1/.

## Toolbox tools
### Correlate Scaled Terrain Indices/Surface Derivatives with Soil Parameters
The purpose of this tool is to find the optimal scale of terrain indices/surface derivatives for landform classification.  Two data layers required are 1) a digital elevation model of the area of interest and 2) a point shapefile/feature class of soil sampling locations throughout the area of interest that contains data for a certain soil parameter that is influenced by landform elements.  This tool calculates the Pearson or Spearman correlation coefficient between the soil parameter and terrain indices/surface derivatives derived at various scales.  The scale at which a terrain index/surface derivative is most correlated with the soil parameter is deemed the most suitable scale for deriving the individual terrain index/surface derivative.  
### Landform Classification (Evans et al., 2016)
This tools classifies landforms according to the modified Topographic Position Index method described in Evans et al., 2016.  Only a digital elevation model of the area of interest is needed as input for this tool.  Also, the scales or neighborhood radii (in cells) in which the topographic position index and convergence index is most correlated with a soil parameter as determined in the “Correlate Scaled Terrain Indices/Surface Derivatives with Soil Parameters” tool are used as input.  The output is a raster whose cells have been classified as convergent/divergent footslopes, backslopes, and shoulders as well as level.  
### Landform Classification (Pennock et al., 1987)
This tools classifies landforms according to the method described in Pennock et al., 1987.  Only a digital elevation model of the area of interest is needed as input for this tool.  Also, the scales or neighborhood radii (in cells) in which profile curvature and plan is most correlated with a soil parameter as determined in the “Correlate Scaled Terrain Indices/Surface Derivatives with Soil Parameters” tool are used as input.  The output is a raster whose cells have been classified as convergent/divergent footslopes, backslopes, and shoulders as well as level.  
### Scale Terrain Indices/Surface Derivatives
This tool generates terrain index/surface derivative rasters at various scales depending on the input neighborhood radius (in cells).  The purpose of this tool is to provide users with a way to generate the scaled terrain indices/surface derivatives that the above tools are using for analysis.  

## Citations:
Evans, D. A., Williard, K. W., & Schoonover, J. E. (2016). Comparison of Terrain Indices and Landform Classification Procedures in Low-Relief Agricultural Fields. Journal of Geospatial Applications in Natural Resources, 1(1), 1.

Pennock, D. J., Zebarth, B. J., and De Jong, E. (1987). Landform classification and soil distribution in hummocky terrain, Saskatchewan, Canada. Geoderma, 40, 297-315.

## Dependencies
Dependencies for this ArcGIS Python Toolbox include arcpy (comes with ArcGIS), pandas, numpy, matplotlib, and scipy python packages as well as SAGA (System for Automated Geoscientific Analysis) which can be downloaded here: https://sourceforge.net/projects/saga-gis/. 

Go to https://sourceforge.net/projects/scipy/files/scipy/ to get the latest python 2.7 version of scipy

Users can download the other python packages with pip installer. To install pip, go to https://pip.pypa.io/en/stable/installing/#id8, save the get-pip.py file, and run it. Next, include the path to pip in the "Path" system environment variable. In the command line, type "pip install package-name" to install specific package.

Install the version of pandas that is compatible with the numpy version with pip. ArcGIS 10.3 comes with numpy 1.7.1, latest compatible version of pandas is 0.13.1.  To install specific package version with pip from command line: "pip install pandas==0.13.1"

Users can also upgrade numpy with pip from command line: "pip install numpy --upgrade"

To access SAGA command line tools, include path to saga_cmd in system environmental Path variable.

