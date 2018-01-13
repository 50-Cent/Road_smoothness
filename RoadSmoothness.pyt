import arcpy as ap
import numpy as np
import math
import ctypes
import os
import sys


#### Function_____________________________

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


def dPartition(Data, hashData, scale):
    restData = np.empty((0),dtype=np.float32)
    hashLoc = np.empty((0),dtype=np.int32)
    clusterData = np.array([[ Data[0,0],Data[0,1] ]])
    tempData = Data[0,:]
    hashLoc = np.append(hashLoc, hashFinder(hashData, Data[0,0], Data[0,1]))
    count_tmp = 1;
    count_rst = 0;


        
    for count in np.arange(len(Data)-1):
        dst = np.linalg.norm(np.array([Data[0,0], Data[0,1]]) - np.array([Data[count+1,0], Data[count+1,1]]))
        if dst <= scale:
            tempData = np.append(tempData,Data[count+1,:])
            hashLoc = np.append(hashLoc, hashFinder(hashData, Data[count+1,0], Data[count+1,1]))
            clusterData = np.concatenate((clusterData, np.array([[ Data[count+1,0],Data[count+1,1] ]]) ),axis=0)
            count_tmp = count_tmp+1           
        elif dst > scale:
            restData = np.append(restData,Data[count+1,:])
            count_rst = count_rst+1

    restData = np.resize(restData,[len(restData)/Data.shape[1], Data.shape[1]])
    tempData = np.resize(tempData,[len(tempData)/Data.shape[1], Data.shape[1]])

    eqsmV = 0
    if count_tmp >= 3:
        eqsmV = np.max(smooPlot(tempData))


    return (restData ,tempData, clusterData, hashLoc, count_tmp, count_rst, eqsmV) 

def smooPlot(Data):
    ##print(Data.shape[0]) ## DO delete
    row = Data.shape[0]
    col = Data.shape[1]
    AA = np.empty((0),dtype=np.float32)
    for r in np.arange(row):
        for c in np.arange(row):
            t1 = np.array([Data[r,0], Data[r,1]])
            t2 = np.array([Data[c,0], Data[c,1]])
            dst = np.linalg.norm(t1-t2)
            AA = np.append(AA,1./(1+np.square(dst))) ## function = 1/(1+distance^2)
    AA = np.resize(AA,[len(AA)/row, row])
    DD = np.zeros((row,row))
    dsum = np.sum(AA, axis = 1)
    for r in np.arange(row):
        DD[r][r] = dsum[r]


    L = DD-AA
    dt = Data[:,2]
    tmp1 = np.dot(np.dot(dt, L), dt)
    tmp2 = np.dot(dt,dt)
    tmp = np.true_divide(tmp1,tmp2)
    tmp = np.fabs(np.true_divide(tmp,np.sqrt(row)))
    smVc = np.multiply(tmp,np.ones((row,1)))

    return smVc


def polyForm(Data, scl, critic):
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
    smV = 4.0
    for r in np.arange(len(tmparray)):
        if tmparray[r] >= critic:
            smV = 8.0
            

    return (restData, clusterData,Data.shape[0]-restData.shape[0], restData.shape[0], smV)



class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Road Smoothness"
        self.alias = "RoadSmoothness"

        # List of tool classes associated with this toolbox
        self.tools = [RoadSmoothness]


class RoadSmoothness(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Road Smoothness"
        self.description = "Road Smoothness"
        self.canRunInBackground = False
        
        return
    

    def getParameterInfo(self):
       #First parameter: the feature layer containing the selected data
        in_features = ap.Parameter(
            displayName = "Layer containing selected SqueeSAR data",
            name = "in_features",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        #Scale
        scale = ap.Parameter(
            displayName = "Scale in Meter",
            name = "Scale",
            datatype = "GPLong",
            parameterType = "Required",
            direction = "Input")
        scale.value = 800

        #Time:
        time_pt = ap.Parameter(
            displayName = "Smoothness at time",
            name = "time_pt",
            datatype = "Field",
            parameterType = "Required",
            direction = "Input")

        time_pt.parameterDependencies = [in_features.name]
        time_pt.value = "D20140601"

        # Parameter: threshold for warning situation
        warn_thresh = ap.Parameter(
            displayName="Smoothness value to use as warning threshold [a value/day/location]",
            name="warn_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        warn_thresh.value = 0.015

        # Parameter: threshold for critical situation
        critical_thresh = ap.Parameter(
            displayName="Smoothness value to use as critical threshold [a value/day/location]",
            name="critical_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        critical_thresh.value = 0.1

        #Parameter: output layer will contain the smoothness
        out_features = ap.Parameter(
            displayName = "Road/pavement quality map",
            name = "out_layer",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Output")

        #Parameter: output layer will contain the thresholded trend
        thresh_features = ap.Parameter(
            displayName = "Regions with poor health",
            name = "polygon_layer",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Output")


        params = [in_features, scale,time_pt, warn_thresh, critical_thresh, out_features, thresh_features]

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

        return

    def execute(self, parameters, messages):
        """Evaluate the top 0-100 percentile total variation"""

        TV_FLD1 = u"SMOO"
        
        # Tell ArcGIS that he can overwrite outputs
        ap.env.overwriteOutput = True

        # Get inputs
        in_layer = parameters[0].valueAsText
        ap.AddMessage("Input layer: {0}".format(in_layer))
        
        #Get scale        
        scale = parameters[1].valueAsText
        ap.AddMessage("Scale: {0}".format(scale))
        scale = np.int_(scale)

        #time_pt
        time_pt = parameters[2].valueAsText
        ap.AddMessage("time: {0}".format(time_pt))
        #time_pt = int(time_pt)        

        # Get the warning threshold
        warn_thresh = parameters[3].valueAsText
        ap.AddMessage("Warning threshold: {0}".format(warn_thresh))
        warn_thresh = np.float_(warn_thresh)

        # Get the critical threshold
        critical_thresh = parameters[4].valueAsText
        ap.AddMessage("Critical threshold: {0}".format(critical_thresh))
        critical_thresh = np.float_(critical_thresh)

        # Get outputs
        out_layer = parameters[5].valueAsText
        ap.AddMessage("Output layer: {0}".format(out_layer))

        # Get outputs for polygon
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
                ap.management.AddField(out_layer, TV_FLD1, "FLOAT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(TV_FLD1, out_layer))
                return



        # Symbology
        symbol_layer = os.path.dirname(os.path.realpath(__file__)) +u"\RoadSmoothness.lyr"


       

        # Data Acquisition
        ap.AddMessage("Start Data Acquisition: {0}".format(TV_FLD1))    
        count = 0
        locat = np.empty((0),dtype=np.float32)
        for row in ap.da.SearchCursor(in_layer,["SHAPE@XY"]):
            x,y = row[0]
            locat = np.append(locat,row[0])
            count+=1
        locat = np.resize(locat,[len(locat)/2, 2])
        del row

        # Get the list of field with date format in the output layer
        field_names = time_pt
        # Dataset preparation
        cursor = ap.da.SearchCursor(in_layer, field_names)
        disp = np.empty((0),dtype=np.float32)
        for row in cursor:
            disp = np.append(disp, row)

        ap.AddMessage("time series dim: {0}".format(disp.shape))

        Data = np.insert(locat,locat.shape[1],disp,axis=1)

        del row
        del cursor
        del disp

        nsData = Data
        Data = sort2D(Data)
        repData = Data   
        
        ap.AddMessage("Start Computing: {0}".format(TV_FLD1))

        # Main Program
        hashData = np.append(nsData[:,0], nsData[:,1])
        hashData = (np.resize(hashData, [2,nsData.shape[0]])).T
        finalData = np.ones((Data.shape[0]))
        count_rst = Data.shape[0]
        
        ap.env.overwriteOutput = True
        
        while len(Data) > 0:
            (Data, newData, clusterData, hashLoc, count_tmp, count_rst, lmn) = dPartition(Data, hashData, scale)
            if count_tmp > 1:
                finalData[hashLoc] = smooPlot(newData)
            elif count_tmp == 1:
                finalData[hashLoc] = 0.0

            if count_rst == 1:
                hashEnt = hashFinder(hashData, Data[0,0], Data[0,1])
                finalData[hashEnt] = 0.0
                Data = np.empty((0),dtype=np.float32)


        del hashLoc
        del count_tmp
        del count_rst
        del newData

        ap.AddMessage("Start Update: {0}".format(TV_FLD1))
        # Update the values in TV_FLD1
        rows = ap.UpdateCursor(out_layer)
        count = 0
        for r in rows:
            r.setValue(TV_FLD1, finalData[count])
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

        ap.AddMessage("Finish symbology: {0}".format(symbol_layer))


##########
##########
##########
##########
##########
########## Polygon Construction       
##########
##########
##########
##########
##########
##########


        polyData = np.empty((0),dtype=np.float_)
        count = 0
        for r in np.arange(nsData.shape[0]):
            if np.fabs(finalData[r]) >= warn_thresh:
                ap.AddMessage("smoothness : {0}".format(finalData[r]))
                ax = np.array([nsData[r,0], nsData[r,1], np.fabs(finalData[r])])
                count+=1
                polyData = np.append(polyData, ax)
                

        polyData = np.resize(polyData,[len(polyData)/3,3])
        ap.AddMessage("smoo value: {0}".format(finalData.shape))
        ap.AddMessage("thresholded data: {0}".format(polyData.shape))

        arSmoo = np.empty((0),dtype=np.float_) 
        features = []
        ap.env.overwriteOutput = True

        while len(polyData) > 0:
            [polyData, cluster, count_clstr, countData, smV] = polyForm(polyData, scale, critical_thresh)                  
            if count_clstr >= 3:
                features.append(ap.Polygon(ap.Array([ap.Point(*coords) for coords in cluster])))
                count+=1
                arSmoo = np.append(arSmoo,smV)
                ap.AddMessage("cluster: {0}".format(count_clstr))
                ap.AddMessage("smoothness : {0}".format(smV))
            if countData <= 2:
                polyData = np.empty((0),dtype=np.float_ )
        
        



##        new_class = os.path.dirname(os.path.realpath(__file__)) +"\polygonSmooFile" +u"\polygonSmooShape.shp"
        dr = os.path.dirname(os.path.realpath(__file__))+ "\\tmpPolyDir"
        try:
            os.stat(dr)
        except:
            os.mkdir(dr)
            

        new_class = dr+u"\\polygonShape.shp"
        
        ap.CopyFeatures_management(features, new_class)
        ap.management.MakeFeatureLayer(new_class, thresh_layer)





        # Symbology
        symbol_thresh_layer = os.path.dirname(__file__)+ u"\PolygonRoadSmoothness.lyr"


        TV_FLD2 = u"EQSMOO"
        ap.env.overwriteOutput = True
        if ap.ListFields(thresh_layer, TV_FLD2):
            ap.AddWarning("Overwriting existing {0} field values".format(TV_FLD2))
        else:            
            try:
                ap.management.AddField(thresh_layer, TV_FLD2, "FLOAT")
            except:
                ap.AddError("Failure to create {0} field in {1}".format(TV_FLD2, thresh_layer))
                return
        



        ap.AddMessage("Start Polygon construction: {0}".format(TV_FLD2))
        #Update SLP values in thresh_layer
        rows = ap.UpdateCursor(thresh_layer)
        count = 0
        for r in rows:
            r.setValue(TV_FLD2, arSmoo[count])
            rows.updateRow(r)
            count+=1
          
        thresh_sym_layer_obj = ap.mapping.Layer(symbol_thresh_layer)
        # Fixing broken link of the polygon symbology layer
##        new_path1 = os.path.dirname(os.path.realpath(__file__)) + "\\DemoRoadQualityData\\"
##        thresh_sym_layer_obj.replaceDataSource(new_path1,"SHAPEFILE_WORKSPACE", "polygonData")
        new_path1 = dr+"\\"
        thresh_sym_layer_obj.replaceDataSource(new_path1,"SHAPEFILE_WORKSPACE", "polygonShape")
        
        ap.management.ApplySymbologyFromLayer(thresh_layer, thresh_sym_layer_obj)
        thresh_layer_obj = ap.mapping.Layer(thresh_layer)
        if thresh_layer_obj.symbologyType == "GRADUATED_COLORS":
            thresh_layer_obj.symbology.valueField = TV_FLD2
            thresh_layer_obj.symbology.classBreakValues = [0, 3.9, 7.9, 10.0]
            thresh_layer_obj.symbology.classBreakLabels = ["NORMAL", "WARNING", "CRITICAL"]

        ap.AddMessage("Finish symbology: {0}".format(symbol_thresh_layer))






        return

   

