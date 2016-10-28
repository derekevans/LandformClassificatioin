#######################################################################################################################################
##This toolbox contains tools that can be use to classify landforms in low relief agriculture fields as described in Evans et al., 2016.
##
##Created by Derek Evans
##
##Full citation for article:
##Evans, D. A., Williard, K. W., & Schoonover, J. E. (2016). Comparison of Terrain Indices and Landform Classification
##Procedures in Low-Relief Agricultural Fields. Journal of Geospatial Applications in Natural Resources, 1(1), 1.
##
##Download full text article at http://scholarworks.sfasu.edu/j_of_geospatial_applications_in_natural_resources/vol1/iss1/1/
##
##Dependencies for this ArcGIS Python Toolbox include arcpy (comes with ArcGIS), pandas, numpy, matplotlib, and scipy python packages
##as well as SAGA (System for Automated Geoscientific Analysis) which can be downloaded here: https://sourceforge.net/projects/saga-gis/
##
##Go to https://sourceforge.net/projects/scipy/files/scipy/ to get the latest python 2.7 version of scipy
##
##Can download the other python packages with pip installer
##To install pip, go to https://pip.pypa.io/en/stable/installing/#id8, save the get-pip.py file, and run it
##Next, include the path to pip in the "Path" system environment variable
##In the command line, type "pip install package-name" to install specific package
##
##Install the version of pandas that is compatible with the numpy version with pip
##ArcGIS 10.3 comes with numpy 1.7.1, latest compatible version of pandas is 0.13.1
##
##To install specific package version with pip from command line: "pip install pandas==0.13.1"
##Can also upgrade numpy with pip from command line: "pip install numpy --upgrade"
##
##To access SAGA command line tools, include path to saga_cmd in system environmental Path variable
#######################################################################################################################################

import arcpy
from arcpy.sa import *
from arcpy import env
import pandas as pd
from os import path
import numpy
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr,spearmanr
import subprocess as sp

#this class is used to scale terrain indices/surface derivatives to the specified neighborhood radius in cells.
class Scale(object):
    def __init__(self):
        pass
    #scale profile curvature
    def profile_curvature(self,dem,cell_radius,out_raster):
        #generate profile curvature raster
        total_curv = Curvature(dem,1,path.join(env.scratchGDB, "prof_curv"),"#")
        prof_curv = Raster(path.join(env.scratchGDB, "prof_curv"))
        
        #scale by calculating focal mean or profile curvature using neighborhood radius
        prof_curv_fm = FocalStatistics(prof_curv,NbrCircle(cell_radius,"CELL"),"MEAN","DATA")
        
        #clip output to input dem
        prof_curv_fm_ext = ExtractByMask(prof_curv_fm, dem)
        prof_curv_fm_ext.save(out_raster)
        
        #clean up
        arcpy.Delete_management(prof_curv)
        arcpy.Delete_management(prof_curv_fm)
        arcpy.Delete_management(total_curv)
        
        return Raster(out_raster)

    #scale plan curvature
    def plan_curvature(self,dem,cell_radius,out_raster):
        #generate plan curvature raster
        total_curv = Curvature(dem,1,"#",path.join(env.scratchGDB, "plan_curv"))
        plan_curv = Raster(path.join(env.scratchGDB, "plan_curv"))

        #scale by calculating focal mean of plan curvature raster using neighborhood radius
        plan_curv_fm = FocalStatistics(plan_curv,NbrCircle(cell_radius,"CELL"),"MEAN","DATA")

        #clip output to input dem
        plan_curv_fm_ext = ExtractByMask(plan_curv_fm, dem)
        plan_curv_fm_ext.save(out_raster)

        #clean up
        arcpy.Delete_management(plan_curv)
        arcpy.Delete_management(plan_curv_fm)
        arcpy.Delete_management(total_curv)
        
        return Raster(out_raster)

    #scale slope
    def slope(self,dem,cell_radius,out_raster):

        #generate slope raster
        slope = Slope(dem,"PERCENT_RISE",1)

        #scale by calculating focal mean of slope raster using neighborhood radius
        slope_fm = FocalStatistics(slope,NbrCircle(cell_radius,"CELL"),"MEAN","DATA")

        #clip raster to input dem
        slope_fm_ext = ExtractByMask(slope_fm, dem)
        slope_fm_ext.save(out_raster)

        #clean up
        arcpy.Delete_management(slope)
        arcpy.Delete_management(slope_fm)
        return Raster(out_raster)

    #generate topographic position index
    def topographic_position_index(self,dem,cell_radius,out_raster):

        #generate focal mean of dem using input neighborhood radius
        dem_fm = FocalStatistics(dem,NbrCircle(cell_radius,"CELL"),"MEAN","DATA")

        #tpi = dem - dem focal mean
        tpi = Minus(dem,dem_fm)
        tpi.save(out_raster)

        #clean up
        arcpy.Delete_management(dem_fm)
        
        return Raster(out_raster)

    #generate convergence index 
    def convergence_index(self,dem,cell_radius,out_raster):

        #get the spatial reference of the input dem (used later to set the output ascii raster from SAGA)
        dem_sr = arcpy.Describe(dem).spatialReference

        #paths to scratch folder for dem and ci ascii and sgrd files
        dem_ascii_path = str(env.scratchFolder) + "\\dem.asc"
        dem_sgrd_path = str(env.scratchFolder) + "\\dem.sgrd"
        ci_sgrd_path = str(env.scratchFolder) + "\\ci.sgrd"
        ci_ascii_path = str(env.scratchFolder) + "\\ci.asc"

        #convert input dem to ascii
        dem_ascii = arcpy.RasterToASCII_conversion(dem,dem_ascii_path)

        #convert ascii to sgrd with SAGA
        #io_grid 1 - Import ESRI Arc/Info Grid
        dem_sgrd = sp.Popen(['saga_cmd','io_grid','1','-GRID='+dem_sgrd_path,'-FILE='+ dem_ascii_path],stdout=sp.PIPE, shell=True)
        dem_sgrd.wait()

        #derive the convergence index with SAGA
        #ta_morphometry 2 - Convergence Index (Search Radius)
        ci_sgrd = sp.Popen(['saga_cmd', 'ta_morphometry', '2', '-ELEVATION='+dem_sgrd_path,'-CONVERGENCE='+ci_sgrd_path,
                            '-RADIUS='+str(cell_radius)],stdout=sp.PIPE, shell=True)
        ci_sgrd.wait()

        #convert the convergence index sgrd to ascii with SAGA
        #io_grid 0 - Export ESRI Arc/Info Grid
        ci_ascii = sp.Popen(['saga_cmd', 'io_grid', '0', '-GRID='+ci_sgrd_path, '-FILE='+ci_ascii_path],
                  stdout=sp.PIPE, shell=True)
        ci_ascii.wait()

        #convert the converence index ascii to ESRI grid
        ci = arcpy.ASCIIToRaster_conversion(ci_ascii_path,out_raster,"FLOAT")

        #set the spatial reference of the convergence index raster to that of the input dem
        ci = arcpy.DefineProjection_management(ci,dem_sr)
        return Raster(out_raster)
            
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Landform Classification Tools (Evans et al., 2016)"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [ScaleTerrainIndices,CorrelateSoilParamsAndIndices,LandformClassificationEvans,LandformClassificationPennock]

#This is the tool that interfaces directly with the Scale class to output the scaled terrain indices/surface derivative
#Output from this tool are similar to Figure 5 in Evans et al., 2016
class ScaleTerrainIndices(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Scale Terrain Indices/Surface Derivatives"
        self.description = "This tool generates the specified terrain index/surface derivative at the specified neighborhood radius."
        self.canRunInBackground = False
        
    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
                displayName="Digital Elevation Model",
                name="dem",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Input")
        
        param1 = arcpy.Parameter(
                displayName="Terrain Index",
                name="terrain_index",
                datatype="String",
                parameterType="Required",
                direction="Input")
        
        param1.filter.type = "ValueList"
        param1.filter.list = ["TOPOGRAPHIC_POSITION_INDEX","CONVERGENCE_INDEX","PROFILE_CURVATURE","PLAN_CURVATURE","SLOPE"]

        param2 = arcpy.Parameter(
                displayName="Neighborhood Radius (cells)",
                name="cell_radius",
                datatype="Long",
                parameterType="Required",
                direction="Input")

        param3 = arcpy.Parameter(
                displayName="Output Raster",
                name="out_raster",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Output")
        
        param3.parameterDependencies = [param0.name]
        params = [param0, param1, param2, param3]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        #Make sure the neighborhood radius is greater than zero
        if parameters[2].altered:
            nr = parameters[2].value
            if nr <= 0:
                parameters[2].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[2].clearMessage()
        return

    def execute(self, parameters, messages):
        
        dem = parameters[0].valueAsText
        terrain_index = parameters[1].valueAsText
        cell_radius = parameters[2].value
        out_raster = parameters[3].valueAsText

        env.overwriteOutput = True

        scale = Scale()

        #generate the user defined terrain index at the neighborhood radius specified
        if terrain_index == "PROFILE_CURVATURE":
            return scale.profile_curvature(dem,cell_radius,out_raster)
        elif terrain_index == "PLAN_CURVATURE":
            return scale.plan_curvature(dem,cell_radius,out_raster)
        elif terrain_index == "SLOPE":
            return scale.slope(dem,cell_radius,out_raster)
        elif terrain_index == "TOPOGRAPHIC_POSITION_INDEX":
            return scale.topographic_position_index(dem,cell_radius,out_raster)
        elif terrain_index == "CONVERGENCE_INDEX":
            return scale.convergence_index(dem,cell_radius,out_raster)

#This tool correlates a soil parameter from an input point feature class/shapefile to terrain indices derived at various neighborhood radius
#Prints to the script tool output the neighborhood radius where the scaled terrain index and soil parameter have maximum correlation
#Optional output is a csv of correlation results and a graph of soil parameter correlation vs terrain index, like Fig. 4 in Evans et al., 2016
class CorrelateSoilParamsAndIndices(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Correlate Scaled Terrain Indices with Soil Parameter"
        self.description = ""
        self.canRunInBackground = False
        
    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
                displayName="Digital Elevation Model",
                name="dem",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Input")
        
        param1 = arcpy.Parameter(
                displayName="Soil Sampling Point Shapefile or Feature Class",
                name="point_shp",
                datatype="Feature Layer",
                parameterType="Required",
                direction="Input")
        
        param2 = arcpy.Parameter(
                displayName="Soil Parameter",
                name="soil_parameter",
                datatype="Field",
                parameterType="Required",
                direction="Input")

        param2.filter.list = ['Short', 'Long', 'Float','Single','Double']
        param2.parameterDependencies = [param1.name]

        param3 = arcpy.Parameter(
                displayName="Terrain Indices for Soil Parameter Correlation",
                name="terrain_indices",
                datatype="String",
                parameterType="Required",
                direction="Input",
                multiValue=True)
        
        param3.filter.type = "ValueList"
        param3.filter.list = ["TOPOGRAPHIC_POSITION_INDEX","CONVERGENCE_INDEX","PROFILE_CURVATURE","PLAN_CURVATURE"]

        param4 = arcpy.Parameter(
                displayName="Minimum Neighborhood Radius (cells)",
                name="min_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")

        param5 = arcpy.Parameter(
                displayName="Maximum Neighborhood Radius (cells)",
                name="max_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")

        param6 = arcpy.Parameter(
                displayName="Neighborhood Radius Interval (cells)",
                name="interval_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")
        
        param7 = arcpy.Parameter(
                displayName="Correlation Type",
                name="correlation_type",
                datatype="String",
                parameterType="Required",
                direction="Input")
        
        param7.filter.type = "ValueList"
        param7.filter.list = ["PEARSON","SPEARMAN"]

        param8 = arcpy.Parameter(
                displayName="Soil Parameter and Scaled Terrain Indices Correlation Results",
                name="output_csv",
                datatype="File",
                parameterType="Optional",
                direction="Output")
        
        param9 = arcpy.Parameter(
                displayName="Output Terrain Indices and Correlation Graph",
                name="output_png",
                datatype="File",
                parameterType="Optional",
                direction="Output")

        params = [param0, param1, param2, param3, param4, param5, param6, param7, param8, param9]
        
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
              
        #make sure the extension on the output is .csv
        if parameters[8].altered:
            output_csv = parameters[8].valueAsText
            path_str,ext = path.splitext(output_csv)
            if ext == "" or ext != ".csv":
                parameters[8].value = path_str + ".csv"

        #make sure the extension on the output is .png
        if parameters[9].altered:
            output_png = parameters[9].valueAsText
            path_str,ext = path.splitext(output_png)
            if ext == "" or ext != ".png":
                parameters[9].value = path_str + ".png"

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        #Make sure the minimum neighborhood radius is greater than zero
        if parameters[4].altered:
            min_nr = parameters[4].value
            if min_nr <= 0:
                parameters[4].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[4].clearMessage()

        #Make sure the maximum neighborhood radius is greater than zero and is greater than the minimum neighborhood radius
        if parameters[5].altered:
            max_nr = parameters[5].value
            if max_nr <= 0:
                parameters[5].setErrorMessage("Input must be greater than 0!")
            elif max_nr <= min_nr:
                parameters[5].setErrorMessage("Maximum neighborhood radius must be greater than minimum neighborhood radius!")
            else:
                parameters[5].clearMessage()
                
        #Make sure the neighborhood radius interval is greater than zero
        if parameters[6].altered:
            nr_interval = parameters[6].value
            if nr_interval <= 0:
                parameters[6].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[6].clearMessage()
                
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        point_shp = parameters[1].valueAsText
        soil_parameter = parameters[2].valueAsText
        terrain_indices = parameters[3].valueAsText.split(";")
        min_nr = parameters[4].value
        max_nr = parameters[5].value
        interval_nr = parameters[6].value
        correlation_type = parameters[7].valueAsText
        output_csv = parameters[8].valueAsText
        output_png = parameters[9].valueAsText
        
        scale = Scale()
        env.overwriteOutput = True

        #reduce the fields in the shapefile or feature class to only include the soil parameter
        def reduce_fields(shp,keep_field):
            fields = arcpy.ListFields(shp)
            delete_fields=[]
            for field in fields:
                if not (keep_field==format(field.name)or field.required):
                    delete_fields.append(format(field.name))
            return arcpy.DeleteField_management(shp,delete_fields)

        #convert feature class or shapefile attribute data to a dataframe
        def feature_to_df(fc):
            fields = arcpy.ListFields(fc)                         
            df_fields=[]
            for field in fields:
                if not field.required:
                    df_fields.append(format(field.name))
                
            df_dict = dict()
            for field in df_fields:
                fc_cursor = arcpy.SearchCursor(fc)
                df_dict[str(field)] = []
                for row in fc_cursor:
                    df_dict[field].append(row.getValue(field))

            df = pd.DataFrame(df_dict)
            return df

        #format p value to make it more readable
        def format_p_value(p):
            if p < 0.0001:
                return "<0.0001"
            else:
                return str("{:.4f}".format(p))
        
        #make a copy of the input point feature class or shapefile and reduce fields
        arcpy.AddMessage("Copying shapefile...")
        point_copy = arcpy.CopyFeatures_management(point_shp,path.join(env.scratchGDB, "point_copy"))
        point_copy = reduce_fields(point_copy,soil_parameter)
        
        #iterate over the indices the user chose, scale them based on the min radius, max radius, and interval,
        #extract values from raster to point copy, and delete the raster to free up memory
        for index in terrain_indices:
            for i in range(min_nr,max_nr+interval_nr,interval_nr):
                arcpy.AddMessage("Calculating " + index + " " + str(i) + " cell radius...")
                if index == "PROFILE_CURVATURE":
                    prof = scale.profile_curvature(dem,i,path.join(env.scratchGDB, "prof_"+str(i)))
                    ExtractMultiValuesToPoints(point_copy,prof)
                    arcpy.Delete_management(prof)
                elif index == "PLAN_CURVATURE":
                    plan = scale.plan_curvature(dem,i,path.join(env.scratchGDB, "plan_"+str(i)))
                    ExtractMultiValuesToPoints(point_copy,plan)
                    arcpy.Delete_management(plan)
                elif index == "TOPOGRAPHIC_POSITION_INDEX":
                    tpi = scale.topographic_position_index(dem,i,path.join(env.scratchGDB, "tpi_"+str(i)))
                    ExtractMultiValuesToPoints(point_copy,tpi)
                    arcpy.Delete_management(tpi)
                elif index == "CONVERGENCE_INDEX":
                    ci = scale.convergence_index(dem,i,path.join(env.scratchGDB, "ci_"+str(i)))
                    ExtractMultiValuesToPoints(point_copy,ci)
                    arcpy.Delete_management(ci)
                    
        #convert the point copy shapefile attributes to a dataframe
        df = feature_to_df(point_copy)

        #iterate over the scaled terrain indices values, correlate with soil parameter, append results to new dataframe
        cols = list(df.columns)
        cols.remove(soil_parameter)
        df_corr = pd.DataFrame()
        for col in cols:
            if correlation_type == "PEARSON":
                r,p = pearsonr(df[soil_parameter],df[col])
            else:
                r,p = spearmanr(df[soil_parameter],df[col])
            index = col.split("_")[0]
            nr = col.split("_")[1]
            df_corr = df_corr.append({"Index": index, "Neighborhood Radius (cells)": int(nr), "r": r, "p-value": p},ignore_index = True)

        #map the df_corr indices to new names for output
        col_name_map = {"ci":"Convergence Index","prof":"Profile Curvature","plan":"Plan Curvature", "tpi":"Topographic Position Index", "slope": "Slope"}
        df_corr["Index"] = df_corr["Index"].map(col_name_map)

        #take absolute value of r to use for finding max r
        df_corr["abs_r"] = abs(df_corr["r"])

        #iterate over indices in dataframe, find row with max r, and print info to output
        indices = pd.unique(df_corr["Index"])
        for index in indices:
            max_row = df_corr.iloc[df_corr[df_corr["Index"]==index]["abs_r"].idxmax()]
            arcpy.AddMessage("--------------------------------------------------------------------")
            arcpy.AddMessage(soil_parameter + " and " + index + " have maximum correlation at " + str(max_row["Neighborhood Radius (cells)"]) + " cell neighborhood radius.")
            if correlation_type == "PEARSON":
                arcpy.AddMessage("Pearson Correlation Coefficient (r): " + "{0:.3f}".format(max_row["r"]))
                arcpy.AddMessage("p-value: " + format_p_value(max_row["p-value"]))
                arcpy.AddMessage("--------------------------------------------------------------------")
            else:
                arcpy.AddMessage("Spearman Correlation Coefficient (r): " + "{0:.3f}".format(max_row["r"]))
                arcpy.AddMessage("p-value: " + format_p_value(max_row["p-value"]))
                arcpy.AddMessage("--------------------------------------------------------------------")

        #if user chooses so, output correlation results to csv
        if output_csv:
            df_corr = df_corr.drop(labels=["abs_r"],axis=1)
            df_corr.to_csv(output_csv,",",index=False)

        #if user chooses so, output plot correlation results and output to png
        if output_png:
            #get min and max r before pivoting
            r_min = df_corr["r"].min()
            r_max = df_corr["r"].max()
             
            df_corr = df_corr.pivot(index='Neighborhood Radius (cells)', columns='Index').r
            
            #set pyplot properties
            plt.xticks(size=12)
            plt.yticks(size=12)
            plt.xlabel("Neighborhood Radius $cells$",fontsize=16, labelpad=15)
            plt.ylabel("Correlation Coefficient $r$",fontsize=16, labelpad=15)
            plt.xlim([min_nr-1,max_nr+1])
            plt.ylim([r_min-0.05,r_max+0.05])

            #styles for lines and markers
            styles = {"Convergence Index":['-','v','k'],"Topographic Position Index": ['--','o','k'], "Profile Curvature": [':','o','w'], "Plan Curvature": ['-.','v','w']}

            #plot the data for each terrain index
            for col in df_corr.columns:
                plt.plot(df_corr.index,df_corr[col],linestyle=styles[col][0],marker=styles[col][1],
                         color='0.75',linewidth=2.0,markeredgecolor='k',label=col,
                         markeredgewidth=1.5,markersize=8,markerfacecolor=styles[col][2])
                
            #set legend position and styles
            plt.legend(loc='upper center', bbox_to_anchor=(0., 1.09, 1., .102),ncol=2,borderaxespad=0.0, fontsize=14)

            #save plot
            plt.savefig(output_png,bbox_inches='tight')

            #clear plot from memory
            plt.close("all")
            return
        
#Classify landforms using the method in Evans et al., 2016.  Output is similar to Evans et al., 2016 Fig. 6.       
class LandformClassificationEvans(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Landform Classification (Evans et al., 2016)"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
                displayName="Digital Elevation Model",
                name="dem",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Input")

        param1 = arcpy.Parameter(
                displayName="Topographic Position Index Neighborhood Radius (cells)",
                name="tpi_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")

        param2 = arcpy.Parameter(
                displayName="Convergence Index Neighborhood Radius (cells)",
                name="ci_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")
        
        param3 = arcpy.Parameter(
                displayName="Output Raster",
                name="out_raster",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Output")
        
        params = [param0,param1,param2,param3]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        #Make sure the tpi neighborhood radius is greater than zero
        if parameters[1].altered:
            tpi_nr = parameters[1].value
            if tpi_nr <= 0:
                parameters[1].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[1].clearMessage()

        #Make sure the ci neighborhood radius is greater than zero
        if parameters[2].altered:
            ci_nr = parameters[2].value
            if ci_nr <= 0:
                parameters[2].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[2].clearMessage()

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        dem = parameters[0].valueAsText
        tpi_nr = parameters[1].value
        ci_nr = parameters[2].value
        out_raster = parameters[3].valueAsText

        scale = Scale()
        #generate tpi
        tpi = scale.topographic_position_index(dem,tpi_nr,path.join(env.scratchGDB, "tpi" + "_" + str(tpi_nr)))

        #get the tpi mean and standard deviation
        tpi_mean = float(arcpy.GetRasterProperties_management(tpi, "MEAN").getOutput(0))
        tpi_std = float(arcpy.GetRasterProperties_management(tpi, "STD").getOutput(0))

        #get the lower and upper bounds for landform classification
        tpi_lower = tpi_mean-(tpi_std*0.75)
        tpi_upper = tpi_mean+(tpi_std*0.75)

        #classify tpi < tpi_lower as footslope
        foot = Con(tpi < tpi_lower,10,0)

        #classify tpi >= tpi_upper as shoulder
        shoulder = Con(tpi >= tpi_upper,30,0)

        #classify tpi_lower <= tpi < tpi_upper as linear areas
        linear = Con((tpi >= tpi_lower) & (tpi < tpi_upper),1,0)

        #generate slope at neighborhood radius as tpi and reclassify areas 0-1% and 1-100%
        slope_map = RemapRange([[0, 1, 1], [1, 100, 2]])
        slope = scale.slope(dem,tpi_nr,path.join(env.scratchGDB, "slope" + "_" + str(tpi_nr)))
        slope_reclass = Reclassify(slope, "VALUE", slope_map)

        #if tpi cell is linear and on 0-1% slope, then classify as level
        level = Con((linear == 1) & (slope_reclass == 1),40,0)

        #if tpi cell is linear and on 1-100% slope, classify as backslope
        back = Con((linear == 1) & (slope_reclass == 2),20,0)

        #generate the convergence indes raster
        ci = scale.convergence_index(dem,ci_nr,path.join(env.scratchGDB, "ci" + "_" + str(ci_nr)))

        #classify cells -100 - 0 as convergent and cells 0 - 100 as divergent
        ci_map = RemapRange([[-100, 0, 1], [0, 100, 2]])
        ci_reclass = Reclassify(ci, "VALUE", ci_map)

        #add all the rasters together to combine flow convergent/divergent and acceleration/deceleration indices. Values map to:
        #11 - Convergent Footslope
        #12 - Divergent Footslope (uncommon in nature landscapes)
        #21 - Convergent Backslope
        #22 - Divergent Backslope
        #31 - Convergent Shoulder (uncommon in nature landscapes)
        #32 - Divergent Shoulder
        #41 - Level
        #42 - Level
        landforms = foot + back + shoulder + level + ci_reclass
        landforms_map = RemapValue([[11,11], [12,12],[21,21],
                             [22,22],[31,31],[32,32],
                             [41,40],[42,40]])
        landforms_reclass = Reclassify(landforms,"VALUE",landforms_map)

        #Add a new field to the raster that contains the landform classes
        arcpy.AddField_management(landforms_reclass,"Landform","TEXT")
        landform_names = """def set_landform_name(x):
          if x == 11:
            return 'Convergent Footslope'
          elif x == 12:
            return 'Divergent Footslope'
          elif x == 21:
            return 'Convergent Backslope'
          elif x == 22:
            return 'Divergent Backslope'
          elif x == 31:
            return 'Convergent Shoulder'
          elif x == 32:
            return 'Divergent Shoulder'
          elif x == 40:
            return 'Level'"""
            
        arcpy.CalculateField_management(landforms_reclass, "Landform", "set_landform_name(!VALUE!)","PYTHON",landform_names)

        landforms_reclass.save(out_raster)
        return Raster(out_raster)

#Classify landforms using the method in Pennock et al., 1987.  Output is similar to Evans et al., 2016 Fig. 6. 
class LandformClassificationPennock(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Landform Classification (Pennock et al., 1987)"
        self.description = "This tool classifies landforms according to Pennock et al., 1987 but using plan and profile curvature at the "
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
                displayName="Digital Elevation Model",
                name="dem",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Input")

        param1 = arcpy.Parameter(
                displayName="Profile Curvature Neighborhood Radius (cells)",
                name="prof_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")

        param2 = arcpy.Parameter(
                displayName="Plan Curvature Neighborhood Radius (cells)",
                name="plan_nr",
                datatype="Long",
                parameterType="Required",
                direction="Input")
        
        param3 = arcpy.Parameter(
                displayName="Output Raster",
                name="out_raster",
                datatype="Raster Layer",
                parameterType="Required",
                direction="Output")
        
        params = [param0,param1,param2,param3]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        #Make sure the profile curvature neighborhood radius is greater than zero
        if parameters[1].altered:
            prof_nr = parameters[1].value
            if prof_nr <= 0:
                parameters[1].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[1].clearMessage()

        #Make sure the plan curvature neighborhood radius is greater than zero
        if parameters[2].altered:
            plan_nr = parameters[2].value
            if plan_nr <= 0:
                parameters[2].setErrorMessage("Input must be greater than 0!")
            else:
                parameters[2].clearMessage()
                
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        dem = parameters[0].valueAsText
        prof_nr = parameters[1].value
        plan_nr = parameters[2].value
        out_raster = parameters[3].valueAsText

        scale = Scale()

        #generate the profile curvature raster
        prof = scale.profile_curvature(dem,prof_nr,path.join(env.scratchGDB, "prof" + "_" + str(prof_nr)))
        
        #invert values so convex is positive and concave is negative
        prof = prof * -1
        
        #classify footslopes as profile curvature < -0.1
        foot = Con(prof < -0.1,10,0)

        #classify shoulder as profile curvature >= 0.1
        shoulder = Con(prof >= 0.1,30,0)

        #classify linear areas as -0.1 <= profile curvature < 0.1
        linear = Con((prof >= -0.1) & (prof < 0.1),1,0)

        #generate slope at neighborhood radius as profile curvature and reclassify areas 0-1% and 1-100%
        slope_map = RemapRange([[0, 1, 1], [1, 100, 2]])
        slope = scale.slope(dem,prof_nr,path.join(env.scratchGDB, "slope" + "_" + str(prof_nr)))
        slope_reclass = Reclassify(slope, "VALUE", slope_map)

        #if profile curvature cell is linear and on 0-1% slope, then classify as level
        level = Con((linear == 1) & (slope_reclass == 1),40,0)

        #if profile curvature cell is linear and on 1-100% slope, then classify as level
        back = Con((linear == 1) & (slope_reclass == 2),20,0)

        #generate the plan curvature raster
        plan = scale.plan_curvature(dem,plan_nr,path.join(env.scratchGDB, "plan" + "_" + str(plan_nr)))

        #classify cells < 0 as convergent and cells > 0 as divergent
        plan_map = RemapRange([[-99999999, 0, 1], [0, 99999999, 2]])
        plan_reclass = Reclassify(plan, "VALUE", plan_map)

        #add all the rasters together to combine flow convergent/divergent and acceleration/deceleration indices. Values map to:
        #11 - Convergent Footslope
        #12 - Divergent Footslope (uncommon in natural landscapes)
        #21 - Convergent Backslope
        #22 - Divergent Backslope
        #31 - Convergent Shoulder (uncommon in natural landscapes)
        #32 - Divergent Shoulder
        #41 - Level
        #42 - Level
        landforms = foot + back + shoulder + level + plan_reclass
        landforms_map = RemapValue([[11,11], [12,12],[21,21],
                             [22,22],[31,31],[32,32],
                             [41,40],[42,40]])
        landforms_reclass = Reclassify(landforms,"VALUE",landforms_map)

        #Add a new field to the raster that contains the landform classes
        arcpy.AddField_management(landforms_reclass,"Landform","TEXT")
        landform_names = """def set_landform_name(x):
          if x == 11:
            return 'Convergent Footslope'
          elif x == 12:
            return 'Divergent Footslope'
          elif x == 21:
            return 'Convergent Backslope'
          elif x == 22:
            return 'Divergent Backslope'
          elif x == 31:
            return 'Convergent Shoulder'
          elif x == 32:
            return 'Divergent Shoulder'
          elif x == 40:
            return 'Level'"""
            
        arcpy.CalculateField_management(landforms_reclass, "Landform", "set_landform_name(!VALUE!)","PYTHON",landform_names)

        landforms_reclass.save(out_raster)
        return Raster(out_raster)           
