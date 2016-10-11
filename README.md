# LandformClassification
ArcGIS Python Toolbox containing tools for classifying landforms in low-relief agricultural fields as described by Evans et al., 2016.  This article can be downloaded at http://scholarworks.sfasu.edu/j_of_geospatial_applications_in_natural_resources/vol1/iss1/1/.

## Full citation for article:
Evans, D. A., Williard, K. W., & Schoonover, J. E. (2016). Comparison of Terrain Indices and Landform Classification Procedures in Low-Relief Agricultural Fields. Journal of Geospatial Applications in Natural Resources, 1(1), 1.

## Dependencies
Dependencies for this ArcGIS Python Toolbox include arcpy (comes with ArcGIS), pandas, numpy, matplotlib, and scipy python packages as well as SAGA (System for Automated Geoscientific Analysis) which can be downloaded here: https://sourceforge.net/projects/saga-gis/. 

Go to https://sourceforge.net/projects/scipy/files/scipy/ to get the latest python 2.7 version of scipy

Users can download the other python packages with pip installer. To install pip, go to https://pip.pypa.io/en/stable/installing/#id8, save the get-pip.py file, and run it. Next, include the path to pip in the "Path" system environment variable. In the command line, type "pip install package-name" to install specific package

Install the version of pandas that is compatible with the numpy version with pip. ArcGIS 10.3 comes with numpy 1.7.1, latest compatible version of pandas is 0.13.1To install specific package version with pip from command line: "pip install pandas==0.13.1"

Users can also upgrade numpy with pip from command line: "pip install numpy --upgrade"

To access SAGA command line tools, include path to saga_cmd in system environmental Path variable

