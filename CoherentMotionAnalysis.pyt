import arcpy as ap
import numpy as np
import math
import ctypes
import os
import sys




def sort2D (Data):
    temp_loc = np.argsort(Data[:,0])
    Data = Data[temp_loc,:]
    tmpX = Data[0,0]
    count = 0
    arY = np.empty((0),dtype=np.float32)
    while Data[count,0] == tmpX:
       arY = np.append(arY,Data[count,1])
       count = count+1
   
    del temp_loc
    if len(arY) > 1:
        temp_loc = np.argsort(arY);
        temp = Data[0,:]
        Data[0,:] = Data[temp_loc[0],:]
        Data[temp_loc[0],:] = temp
        del temp_loc
        del temp        
    
    return Data


def hashFinder(Data, X, Y):
    ff = (Data[:,0]==X)&(Data[:,1]==Y)
    gg = np.where(ff==True)

    return gg


def breakStrnLoc(strn, pos):
    start = 0
    return (strn[start:start+pos[0]+1], strn[start+pos[0]+1:start+pos[1]+1], strn[start+pos[1]+1:len(strn)])



def dPartition(Data, scl, critic):
    tmparray = np.array([Data[0,2]])
    clusterData = np.array([[ Data[0,0],Data[0,1] ]])
    restData = np.empty((0),dtype=np.float_)
    
    for r in np.arange(Data.shape[0]-1):
        if np.linalg.norm(np.array([Data[0,0], Data[0,1]]) - np.array([Data[r+1,0], Data[r+1,1]])) <= scl :
            clusterData = np.concatenate((clusterData, np.array([[ Data[r+1,0],Data[r+1,1] ]]) ),axis=0)
            tmparray = np.append(tmparray,Data[r+1,2])
        else:
            restData = np.append(restData, Data[r+1,:])

    restData = np.resize(restData,[len(restData)/3,3])
    SLP = 4.0
    for r in np.arange(len(tmparray)):
        if tmparray[r] >= critic:
            SLP = 8.0
            

    return (restData, clusterData,Data.shape[0]-restData.shape[0], restData.shape[0], SLP)
    
    


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Coherent motion analysis"
        self.alias = "CoherentMotionAnalysis"

        # List of tool classes associated with this toolbox
        self.tools = [CoherentMotionAnalysis]



class CoherentMotionAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Coherent motion analysis"
        self.description = "This tool evaluates the slope of the time-series associated with each point"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #First parameter: the feature layer containing the selected data
        in_features = ap.Parameter(
            displayName = "Layer containing selected SqueeSAR data",
            name = "in_features",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
##        mxd = arcpy.mapping.MapDocument(in_features)
        
        
        yr = ap.Parameter(
            displayName="Year to compute the trend in a time series",
            name="Year",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        yr.value = 2014

        # Parameter: threshold for warning situation
        warn_thresh = ap.Parameter(
            displayName="Median slope to use as warning threshold [mm/day]",
            name="warn_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        warn_thresh.value = 0.02

        # Parameter: threshold for critical situation
        critical_thresh = ap.Parameter(
            displayName="Median slope to use as critical threshold [mm/day]",
            name="critical_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        critical_thresh.value = 0.05

        # Parameter: threshold for critical situation
        scale = ap.Parameter(
            displayName="Buffer to construct polygon(in Feet)",
            name="scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        scale.value = 1000

        #Parameter: output layer will contain the trend
        out_features = ap.Parameter(
            displayName = "Coherent motion map",
            name = "out_layer",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Output")

        #Parameter: output layer will contain the thresholded trend
        thresh_features = ap.Parameter(
            displayName = "Active region(s) having coherent motion",
            name = "thresh_layer",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Output")



        params = [in_features, yr, warn_thresh, critical_thresh, scale, out_features, thresh_features]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
##        ap.env.overwriteOutput = True
        
        
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        TV_FLD1 = u"ThSen"
       
        
        # Tell ArcGIS that he can overwrite outputs
        ap.env.overwriteOutput = True

        # Get inputs
        in_layer = parameters[0].valueAsText
        ap.AddMessage("Input layer: {0}".format(in_layer))

        # Get the year of time series
        yr = parameters[1].valueAsText
        ap.AddMessage("Time Series year: {0}".format(yr))


        # Get the warning threshold
        warn_thresh = parameters[2].valueAsText
        ap.AddMessage("Warning threshold: {0}".format(warn_thresh))
        warn_thresh = np.float_(warn_thresh)

        # Get the critical threshold
        critical_thresh = parameters[3].valueAsText
        ap.AddMessage("Critical threshold: {0}".format(critical_thresh))
        critical_thresh = np.float_(critical_thresh)

        # Get the scale for polygon
        scl = parameters[4].valueAsText
        ap.AddMessage("Scale: {0}".format(scl))
        scl = np.float_(scl)

        # Get outputs
        out_layer = parameters[5].valueAsText
        ap.AddMessage("Output layer: {0}".format(out_layer))        

        # Get outputs
        thresh_layer = parameters[6].valueAsText
        ap.AddMessage("Polygon layer: {0}".format(thresh_layer))


        
        
        # Copy the input layer to the output layer
        try:
            ap.management.MakeFeatureLayer(in_layer, out_layer)
        except:
            ap.AddError("Failure to copy {0} to {1}".format(in_layer, out_layer))
            return

        # Check if the TV_FLD1 field already exists in the output layer
        if ap.ListFields(out_layer, TV_FLD1):
            ap.AddWarning("Overwriting existing {0} field values".format(TV_FLD1))
        else:            
            try:
                ap.management.AddField(out_layer, TV_FLD1, "DOUBLE")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(TV_FLD1, out_layer))
                return


        # Symbology
        symbol_layer = os.path.dirname(os.path.realpath(__file__)) +u"\CoherentMotionSymbology.lyr"



        
        # Data Acquisition
        ap.AddMessage("Start Data Acquisition: {0}".format(TV_FLD1))    
        count = 0
        locat = np.empty((0),dtype=np.float_)
        for row in ap.da.SearchCursor(in_layer,["SHAPE@XY"]):
            x,y = row[0]
            locat = np.append(locat,row[0])
            count+=1
        locat = np.resize(locat,[len(locat)/2, 2])
        del row

        # Get the list of field with date format in the output layer
        fld_name = "D"+yr+"*"
        field_names = [f.name for f in ap.ListFields(in_layer, fld_name)]
        ap.AddMessage("fields : {0}".format(field_names[0]))

        # Dataset preparation
        cursor = ap.da.SearchCursor(in_layer, field_names)
        disp = np.empty((0),dtype=np.float_)
        for row in cursor:
            for i in range(0, len(field_names)):
                disp = np.append(disp,row[i])

        disp = np.resize(disp,[len(disp)/len(field_names), len(field_names)])
        ap.AddMessage("time series dim: {0}".format(disp.shape))

        Data = np.concatenate((locat,disp), axis=1)
        ap.AddMessage("Data dim: {0}".format(Data.shape))
        del row
        del cursor
        del disp

        finalData = np.empty((0),dtype=np.float_)

        #Day count using field names
        strnDate = [s.replace("D","") for s in field_names]
        pos = np.array([3,5])
        intDate = [breakStrnLoc(s,pos) for s in strnDate]
        intDate = np.int_(intDate)
        ap.AddMessage("Date : {0}".format(intDate))

        # Theil-Sen linear trend algorithm
        import datetime        
        for r in np.arange(Data.shape[0]):
            tmp = np.empty((0),dtype=np.float_)
            for c in np.arange(Data.shape[1]-2):
                tt1 = datetime.date(intDate[c,0],intDate[c,1],intDate[c,2])
                for s in np.arange(Data.shape[1]-c-3):
                    tt2 = datetime.date(intDate[c+s+1,0],intDate[c+s+1,1],intDate[c+s+1,2])
                    difff = tt2-tt1
                    ff = np.true_divide((Data[r,c+s+3]-Data[r,c+2]),difff.days)
                    tmp = np.append(tmp,ff)

            finalData = np.append(finalData,np.median(tmp))
            

        ap.AddMessage("Start Update: {0}".format(TV_FLD1))
        # Update the values in TV_FLD1
        rows = ap.UpdateCursor(out_layer)
        count = 0
        for r in rows:
            r.setValue(TV_FLD1, np.fabs(finalData[count]))
            rows.updateRow(r)
            count=count+1

        # Properly close the cursor
        del r
        del rows

        ap.AddMessage("Finish Update: {0}".format(TV_FLD1))      
                             

        
        # use symbology layer to apply to your output layer

        sym_layer_obj = ap.mapping.Layer(symbol_layer)

        
        
        # Fixing broken link of the symbology layer
        data_loc = ap.Describe(in_layer).catalogPath
        (new_path, file) = os.path.split(data_loc)
        (file_name, file_ext) = os.path.splitext(file)

        sym_layer_obj.replaceDataSource(str(new_path), "SHAPEFILE_WORKSPACE", str(file_name))
        # The above lines can be wrapped up in a single line using the command "ParsePath_mb"

        
        ap.management.ApplySymbologyFromLayer(out_layer, sym_layer_obj)
        out_layer_obj = ap.mapping.Layer(out_layer)
        if out_layer_obj.symbologyType == "GRADUATED_COLORS":
            out_layer_obj.symbology.valueField = TV_FLD1
            out_layer_obj.symbology.classBreakValues = [0, warn_thresh, critical_thresh, 10.0]
            out_layer_obj.symbology.classBreakLabels = ["NORMAL", "WARNING", "CRITICAL"]


        #polygon feature class construction
        threshData = np.empty((0),dtype=np.float_)
        count = 0

        for r in np.arange(Data.shape[0]):
            if np.fabs(finalData[r]) >= warn_thresh:
                ax = np.array([Data[r,0], Data[r,1], np.fabs(finalData[r]) ])
                count+=1
                threshData = np.append(threshData, ax)

        threshData = np.resize(threshData,[len(threshData)/3,3])
        ap.AddMessage("thresholded data: {0}".format(threshData.shape))

        
######
######
######
######
######
######
######
######
######
######
######

        features = []
        ap.env.overwriteOutput = True
        count = 0
        # SLP = 2.0(Green), 4.0(Yellow), 8.0(Red)
        arSLP = np.empty((0),dtype=np.float_)
        while len(threshData) > 0:
            [threshData, cluster, count_clstr, countData, SLP] = dPartition(threshData, scl, critical_thresh)                  
            if count_clstr >= 3:
                features.append(ap.Polygon(ap.Array([ap.Point(*coords) for coords in cluster])))
                count+=1
                arSLP = np.append(arSLP,SLP)
                ap.AddMessage("cluster: {0}".format(count_clstr))
            if countData <= 2:
                threshData = np.empty((0),dtype=np.float_ )

#        new_class = os.path.dirname(os.path.realpath(__file__)) +"\polygonFile" +u"\polygonShape.shp"
        # Create a directory under new_path called "~tmpPolyDir"
        dr = os.path.dirname(os.path.realpath(__file__))+ "\\tmpPolyDir"
        try:
            os.stat(dr)
        except:
            os.mkdir(dr)
            

        new_class = dr+u"\\polygonShape.shp"
        
        ap.CopyFeatures_management(features, new_class)
        ap.management.MakeFeatureLayer(new_class, thresh_layer)
        
       
        #thresh_layer adjustment with equivalent slope

        TV_FLD2 = u"EQSLP"
        ap.env.overwriteOutput = True
        if ap.ListFields(thresh_layer, TV_FLD2):
            ap.AddWarning("Overwriting existing {0} field values".format(TV_FLD2))
        else:            
            try:
                ap.management.AddField(thresh_layer, TV_FLD2, "FLOAT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(TV_FLD2, thresh_layer))
                return


        # Symbology
        symbol_thresh_layer = os.path.dirname(__file__)+ u"\CoherentMotionPolygonSymbology.lyr"


        

        ap.AddMessage("Start Polygon construction: {0}".format(TV_FLD2))
        #Update SLP values in thresh_layer
        rows = ap.UpdateCursor(thresh_layer)
        count = 0
        for r in rows:
            r.setValue(TV_FLD2, arSLP[count])
            rows.updateRow(r)
            count+=1


          
        thresh_sym_layer_obj = ap.mapping.Layer(symbol_thresh_layer)
        
        # Fixing broken link of the polygon symbology layer
        
##        new_path1 = os.path.dirname(os.path.realpath(__file__)) + "\\CoherentMotionAnalysisDemoData\\"
##        thresh_sym_layer_obj.replaceDataSource(new_path1,"SHAPEFILE_WORKSPACE", "polygonData")


        new_path1 = dr+"\\"
        thresh_sym_layer_obj.replaceDataSource(new_path1,"SHAPEFILE_WORKSPACE", "polygonShape")
   
        ap.management.ApplySymbologyFromLayer(thresh_layer, thresh_sym_layer_obj)
        thresh_layer_obj = ap.mapping.Layer(thresh_layer)
        if thresh_layer_obj.symbologyType == "GRADUATED_COLORS":
            thresh_layer_obj.symbology.valueField = TV_FLD2
            thresh_layer_obj.symbology.classBreakValues = [0, 2.1, 4.1, 10.0]
            thresh_layer_obj.symbology.classBreakLabels = ["NORMAL", "WARNING", "CRITICAL"]

##        ap.AddMessage("Finish symbology: {0}".format(symbol_thresh_layer))

##        mxd = arcpy.mapping.MapDocument("CURRENT")
##        upBound = 4.0
##        ap.SelectLayerByAttribute_management(thresh_layer,"NEW_SELECTION",u"EQSLP > %f"%upBound)
##        df = ap.mapping.ListDataFrames(mxd, "Layers") [0]
##        df.zoomToSelectedFeatures()
##        ap.RefreshActiveView()
        
  
            


        return



        



        