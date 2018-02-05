# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus into a surface in DAPI channel,
#   Get mean and sum DAPI intensity inside surface
#   Get volume of surface
#   Cluster nucleus into ploidy
#   Exports result into .csv tables
# Note: This script is calibrated for 3D images of plant nuclei conterstained in DAPI obtained using Leica TCS SP8
# Creator: Mariamawit S. Ashenafi, Célia Baroux, UZH
# Published on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="NucleiPloidyPollen3" icon="Python" tooltip="NucleiPloidyPollen3">
#         <Command>PythonXT::XTNucleiPloidyPollen3(%i)</Command>
#       </Item>
#      </Menu>
#    </CustomTools>

#python library
import numpy as np
import pandas as pd
import logging
import tkMessageBox
import os
from Tkinter import *
import tkFileDialog
import time
import datetime
import itertools
#==============================================================================
# The following python packages are only installed on my home directory on the virtual machines
# TODO: Tell user to add any additional directoy containing packages in Preferences/Custom Tool
# Note this adds more time to load  all libs required
# Keep this option for zmb users or edit it out before publishing 
#==============================================================================
import sys #TODO: edit out?
sys.path.insert(1,'H:\Python') #TODO: edit out?
#import scipy #TODO: edit out?
#from sklearn.cluster import KMeans #both numpy and scipy needs to be imported before importing sklearns 
#import matplotlib.pyplot as plt
#==============================================================================

#Imaris library
import ImarisLib

#==============================================================================
# Functions required to log and display progress of the plugin
# takes a string as an argument and
#==============================================================================
def logtime(aTitle):
  # type: (object) -> object
  global gLasttime
  curtime = datetime.datetime.now()
  if (gLasttime is not None):
    diff = (curtime-gLasttime).total_seconds()
  else:
    diff = '-'
  gLasttime = curtime
  print (curtime.ctime()+ ' [ '+ str(diff)+ ' ] '+ aTitle)
  
#==============================================================================
# Pop-up windows to ask user to set folder pathway and channels to segment
#==============================================================================
class Checkbar(Frame):
  def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
    Frame.__init__(self, parent)
    self.vars = []
    for pick in picks:
      var = IntVar()
      chk = Checkbutton(self, text=pick, variable=var)
      chk.pack(side=side, anchor=anchor, expand=YES)
      self.vars.append(var)
  def state(self):
    return map((lambda var: var.get()), self.vars)

def allstates(): 
    global User_selection
    global root
    if sum(list(lng.state()))>0:
        User_selection=list(lng.state())
        root.destroy()
    else:
        Message="Please select one of the options."
        Label(root, text=Message).grid(row=1)
 
def PopUpMessage(OPTIONS, Messge):
  global root
  global lng
  root = Tk()
  label_text	=	Messge
  option		=	OPTIONS
  Label(root, text=label_text).grid(row=0)
  lng = Checkbar(root, option)
  lng.grid(row=2)
  lng.config(relief=GROOVE, bd=2)
  Button(root, text='Quit', fg="red", command=quit).grid(row=4)
  Button(root, text='Submit', fg="darkgreen", command=allstates).grid(row=5)
  root.mainloop()

#This function asks user to set some parameters FISH and DAPI channels.
def Ask_user(numberIndex):
    OPTIONS                         =       range(numberIndex)
    Messge="Please select the DAPI channel: \n Please only choose one channel."
    PopUpMessage(OPTIONS, Messge)
    DAPIChannel	=	[i for i, x in enumerate(User_selection) if x == 1][0]
    return DAPIChannel

#Segments a channel into surface and returns the surface
def Segment_Surface(vImage, ch, vSFW, name,vLCF, AreaThreshold):
    vROI = None
    vATA = 1
    vATM = 0
    vSFS = '"Area" above '+AreaThreshold+' um^2'
    vDCI = ch
    vSurface2 = vImaris.GetImageProcessing().DetectSurfaces(vImage, vROI, vDCI,vSFW, vLCF, vATA, vATM, vSFS)
    vSurface2.SetName(name)
    return vSurface2

#Funtion to get sum intensity and volume inside surfaces
def GetNucleusIntensity(surf, TypeSurface):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    if TypeSurface==1:
        VolumeList        = 	[float(vValues[a]) for a, x in enumerate(vNames) if x == "Volume"]
    else: 
        VolumeList        = 	[float(vValues[a]) for a, x in enumerate(vNames) if x == "Area"]
    IntenityList        =   [float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Sum"]
    return VolumeList, IntenityList

#This function returns the image in the form of float value, and get surface table (inside table value=1, outside table value=0)
def Get_Mask_data():
    vImage			=	vImaris.GetDataSet()
    vImageSizeX 	=	vImage.GetSizeX()
    vImageSizeY 	=	vImage.GetSizeY()
    vImageSizeZ 	=	vImage.GetSizeZ()
    vExtentMinX		=	vImage.GetExtendMinX()
    vExtentMinY		=	vImage.GetExtendMinY()
    vExtentMinZ		=	vImage.GetExtendMinZ()
    vExtentMaxX		=	vImage.GetExtendMaxX()
    vExtentMaxY		=	vImage.GetExtendMaxY()
    vExtentMaxZ		=	vImage.GetExtendMaxZ()
    mask_min		=	[vExtentMinX, vExtentMinY, vExtentMinZ]
    mask_max		=	[vExtentMaxX, vExtentMaxY, vExtentMaxZ]
    mask_size		=	[vImageSizeX, vImageSizeY, vImageSizeZ]
    SixeX                   =   vExtentMaxX -   vExtentMinX
    SixeY                   =   vExtentMaxY -   vExtentMinY
    SixeZ                   =   vExtentMaxZ -   vExtentMinZ
    VoxelSizeX              =   vImageSizeX/ SixeX
    VoxelSizeY              =   vImageSizeY/ SixeY
    VoxelSizeZ              =   vImageSizeZ/ SixeZ
    return mask_min, mask_max, mask_size, VoxelSizeX, VoxelSizeY, VoxelSizeZ, vImageSizeX, vImageSizeY, vImageSizeZ

# Function to convert spots coordinates into voxel index
def ConvertCoordianteIntoVoxelId( SpotPosition, VoxelSizeX, VoxelSizeY, VoxelSizeZ):
    Result      =   []
    for spot in SpotPosition:
        VoxelID     =   [int(round(x*y, 0)) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
        Result.append(VoxelID)
    return Result
   
#For a surface containing several surface, this function selectes the most dense surface and creates a surface containing only this selected surface
def SelectSurface(surfPollen, surfNuclei):
    NumberOfNuclei=surfNuclei.GetNumberOfSurfaces()
    ListOfNucleiPosition=[]#this list will contain the pollen ID for each nuclei 
    for nuc in range(NumberOfNuclei):
        CoordOfCenterOfMass     =   surfNuclei.GetCenterOfMass(nuc)[0]
        ListOfNucleiPosition.append(CoordOfCenterOfMass)
    #Get IDs of voxels inside pollen surface
    NumberOfPollen=surfPollen.GetNumberOfSurfaces()
    mask_min, mask_max, mask_size, VoxelSizeX, VoxelSizeY, VoxelSizeZ, vImageSizeX, vImageSizeY, vImageSizeZ             =   Get_Mask_data()
    ListPollenVoxelId=[] #Get voxels inside each pollen grain surface
    for PollenId in range(NumberOfPollen):
        MasKPollen			=	surfPollen.GetSingleMask(PollenId, mask_min[0], mask_min[1], mask_min[2], mask_max[0], mask_max[1], mask_max[2],mask_size[0],mask_size[1], mask_size[2])
        VoxelsInsidePollen		=	MasKPollen.GetDataVolumeFloats(DAPIChannel,0)
        VoxelsInsideSurface    =   [[x,y,z] for x in range(vImageSizeX) for y in range(vImageSizeY) for z in range(vImageSizeZ) if VoxelsInsidePollen[x][y][z]]   
        ListPollenVoxelId.append(VoxelsInsideSurface)
    #Convert nuclei position coordinate into voxelID
    VoxelIDOfNuclei     =   ConvertCoordianteIntoVoxelId(ListOfNucleiPosition, VoxelSizeX, VoxelSizeY, VoxelSizeZ)       
    PollenIDList = [z for x in VoxelIDOfNuclei for z in range(len(ListPollenVoxelId)) if x in ListPollenVoxelId[z]]
    return PollenIDList
    

# ==============================================================================
# Functions required to ask user set immuno signal spot diameter for all channels
# ==============================================================================
def fetch(entries, root):
    global SmothingFactor
    text = entries.get()
    try:
        x = float(text)
    except ValueError:
        print ("You must enter a number in field")
        return
    root.destroy()
    SmothingFactor = x


def makeform(root):
    field = "nucleus segmentation"
    row = Frame(root)
    lab = Label(row, width=30, text=field, anchor='w')
    ent = Entry(row)
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent.pack(side=RIGHT, expand=YES, fill=X)
    return ent

#
#def AskUserSmoothingFactor():
#    root = Tk()
#    label_text = "Set smooth surface detail for "
#    row = Frame(root)
#    lab = Label(row, width=30, text=label_text, anchor='w')
#    row.pack(side=TOP, fill=X, padx=5, pady=5)
#    lab.pack(side=LEFT)
#    ents = makeform(root)
#    root.bind('<Return>', (lambda event, e=ents: fetch(e, root)))
#    b1 = Button(root, fg="darkgreen", text='Submit',
#                command=(lambda e=ents: fetch(e, root)))
#    b1.pack(side=LEFT, padx=5, pady=5)
#    b2 = Button(root, fg="darkred", text='Quit', command=root.quit)
#    b2.pack(side=LEFT, padx=5, pady=5)
#    root.mainloop()


#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTCountSpotPerShell_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)
  

#==============================================================================
# Function to get parameters for this plugin from file:XTSegmentNucleus_Parmaeters.csv
#==============================================================================
def GetPluginParameters():
    currentDirectory                         =       os.getcwd()
    AllFilesInDirectory                      =       os.listdir(currentDirectory)
    ParametersList                           =       None
    if "XTNucleiPloidyPollen_Parameters.csv" in AllFilesInDirectory:
        ParameterData                            =       pd.read_csv("XTNucleiPloidyPollen_Parameters.csv", sep=";", header='infer',decimal='.')
        if "Value" in ParameterData.columns:
            ParametersList                           =       list(ParameterData["Value"])
        else:
            tkMessageBox.showinfo(title="Error", message="Please make sure the 'XTNucleiPloidyPollen_Parameters.csv' file contains a column 'Value' containing the values necessary for this plugin.")
            quit()
    else:
        tkMessageBox.showinfo(title="Error", message="Please make sure there is a 'XTNucleiPloidyPollen_Parameters.csv' in the folder containing the 'XTSegmentNucleus.py'.")
        quit()
    return ParametersList

#==============================================================================
# #Function to export output in .csv files
#==============================================================================
def CreateTables(FileNameList, vData, ExpFileName, Result_pathway,ChID, FormatTable, FileIndex):
    if len(vData)>0:
        SaveTable         =       pd.DataFrame(vData)    #Save position table
        if FormatTable   :
            SaveTable     =       SaveTable.T
        if FileNameList != None:
            data                            =       pd.DataFrame(FileNameList)
            SaveTable                      =       pd.concat([SaveTable, data], axis=1)
        vPathToSaveTables = os.path.join(Result_pathway, ExpFileName+str(ChID))
        SaveTable.to_csv(path_or_buf=vPathToSaveTables+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    logtime('Data saved- image_'+str(FileIndex))
     
#==============================================================================
# This function:
# Segments DAPI chanenl into nucleus surface and gets the mean and max DAPI intensity inside nucleus and the volume of the nucleus
# All data is then registered into IntensityMax.csv,IntensityMean.csv and VolumeNuclei.csv
#==============================================================================
def GetImageFeatures(FileIndex, Result_pathway, vFileName, vFullFileName):
    global vImaris
    global DAPIChannel
    global ParametersList
#    global VolumePollenResult
#    global VolumeNucleiResult
#    global IntensityPollenResult  
#    global IntensityNucleiResult
#    global PollenIndexList
    vImage               =   	vImaris.GetDataSet()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        if FileIndex    == 1:
            ParametersList = GetPluginParameters()
            DAPIChannel= Ask_user(numberIndex)
#            AskUserSmoothingFactor()
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        vImage                    =    vImaris.GetDataSet()
        logtime('START - image_'+str(FileIndex))
#==============================================================================
#         STEP1: Segment nuclei and pollen grains
#==============================================================================
        vScene                      =    vImaris.GetSurpassScene()
        vFactory		              =	    vImaris.GetFactory()
        GroupOfObjects	           =	    vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
        SurfPollen=Segment_Surface(vImage, DAPIChannel, float(ParametersList[0]), "Pollen Grain",0, str(ParametersList[1])) # Segment pollen
        vSigma                      =    ParametersList[5]
        vImaris.GetImageProcessing().GaussFilterChannel(vImage,DAPIChannel,vSigma) #Apply Gaussian smoothing filter
        SurfNucleus=Segment_Surface(vImage, DAPIChannel, float(ParametersList[2]), "Nuclei",float(ParametersList[3]), str(ParametersList[4])) # Segment nuclei
        GroupOfObjects.AddChild(SurfNucleus, -1)
        GroupOfObjects.AddChild(SurfPollen, -1)
        vScene.AddChild(GroupOfObjects, -1)
        logtime('Nucleus surface segmented - Image '+str(FileIndex))
#==============================================================================
#         STEP2: Get volume and intensity sum of every surface created
#==============================================================================
        VolumeListPollen, IntensityListPollen   =   GetNucleusIntensity(SurfPollen, 0)
        VolumeListNuclei, IntensityListNuclei   =   GetNucleusIntensity(SurfNucleus, 1)
        ListOuput=[VolumeListPollen, IntensityListPollen, VolumeListNuclei, IntensityListNuclei]
        OuputTtile=["PollenArea", "PollenIntensity", "NucleiVolume", "NucleiIntensity"]
        for outName, outValue in zip(OuputTtile, ListOuput):
            CreateTables(None, outValue, outName, Result_pathway, "", 0, FileIndex)        
#        VolumePollenResult.append(VolumeListPollen)
#        VolumeNucleiResult.append(VolumeListNuclei)
#        IntensityPollenResult.append(IntensityListPollen)  
#        IntensityNucleiResult.append(IntensityListNuclei)
        logtime('Get intensity and volume for each surface - Image '+str(FileIndex))
#==============================================================================
#         STEP3: Get ID of pollen for each nuclei (group nuclei by pollen grain surface)
#==============================================================================
#        PollenIDList=SelectSurface(SurfPollen, SurfNucleus) 
#        CreateTables(None, PollenIDList, "PollenIDList", Result_pathway, "", 0, FileIndex)
##        PollenIndexList.append(PollenIDList)
#        logtime('Get pollen id for each nucleus surface - Image '+str(FileIndex))
#==============================================================================
#         STEP4: Save file and snapshot, and remove orginal image
#==============================================================================
        vPath = os.path.join(Result_pathway, vFileName+".ims") #Save image file containing nucleus surface
        vImaris.FileSave(vPath, "")
        vScene.RemoveChild(GroupOfObjects) #remove containter from scene
        os.remove(vFullFileName)
        logtime('Image saved - Image '+str(FileIndex))
        logtime('END - image_'+str(FileIndex))

def BarplotFigure(ClusterIndex, SumIntensityList,Result_pathway):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ## the data
    UniqueValuesClusters        = list(set(ClusterIndex))
    NumberOfImages              = len(SumIntensityList)
    NumberOfNuclei              = [float(ClusterIndex.count(n)) for n in UniqueValuesClusters]
    PercentageOfNucleiPerCluster    = [(float(ClusterIndex.count(n)) / float(NumberOfImages)) * 100.0 for n in
                                UniqueValuesClusters]
    AverageIntensity            = []
    IntensityStd                = []
    for clusterID in UniqueValuesClusters:
        SelectedValuesIntensity = [SumIntensityList[x] for x in range(NumberOfImages) if ClusterIndex[x] == clusterID]
        AverageIntensity.append(np.mean(SelectedValuesIntensity))
        IntensityStd.append(np.std(SelectedValuesIntensity))
    #Order clusters according to DAPI intensity average value 
    SortedClusterIndexValue            =   [x for (y,x) in sorted(zip(AverageIntensity,UniqueValuesClusters))]
    n=0
    for i in SortedClusterIndexValue:
        UniqueValuesClusters[i]                    =   "P"+str(n)
        n+=1
    ind = np.array(range(len(UniqueValuesClusters)))  # the x locations for the groups
    width = 0.35  # the width of the bars
    ## the bars
    rects1 = ax2.bar(ind, AverageIntensity, width,
                     color='olive',
                     yerr=IntensityStd,
                     error_kw=dict(elinewidth=2, ecolor='olive'))
    rects2 = ax.bar(ind + width, PercentageOfNucleiPerCluster, width,
                    color='mediumslateblue')
    # axes and labels
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0.0, 100.0)
    ax.set_ylabel('Percentage of nuclei', color="mediumslateblue")
    ax.tick_params(axis='y', colors='mediumslateblue')
    ax.tick_params(axis='x', colors='white')
    ax2.set_ylabel('Average DAPI intensity sum', color="olive")
    ax2.tick_params(axis='y', colors='olive')
    ax.set_title('DAPI intensity sum and percentage of nuclei by cluster', color="white")
#    xTickMarks = ['Cluster' + str(i) for i in range(1, len(UniqueValuesClusters) + 1)]
    xTickMarks = UniqueValuesClusters
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ## add a legend
    ax.legend((rects1[0], rects2[0]), ('Intensity sum average', 'Number of nuclei'), fontsize=8)
    PlotFile = os.path.join(Result_pathway, "PloidyPlot.png") 
    fig.savefig(PlotFile, transparent=True)
    plt.close(fig)
    dataR       =       pd.DataFrame({"Cluster": UniqueValuesClusters, "SumIntensity Average":  AverageIntensity, 'SumIntensity STD': IntensityStd, 'Number of Nuclei': NumberOfNuclei})
    vPathToSaveTables = os.path.join(Result_pathway,"ClusterIntensityAverage.csv")
    dataR.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    return UniqueValuesClusters

#==============================================================================
# Function required to concatenate files containing the same information into one table
#==============================================================================
def ConcatenatedTables(FileNameList, Result_pathway,AllFilesWithParam,FileType):
    DataResult                          =       pd.DataFrame()
    for File in AllFilesWithParam:
        SelectedFile                    =       os.path.join(Result_pathway, File)
        data                            =       pd.read_csv(SelectedFile, header='infer',decimal='.')
        DataResult                      =       pd.concat([DataResult, data], axis=1)
        os.remove(SelectedFile)
    if len(DataResult)>0:
        NewColnames                     =       [x+str(y) for x in FileNameList for y in data.columns] 
        DataResult.columns              =       NewColnames
        vPathToSaveTables               =       os.path.join(Result_pathway, FileType)
        DataResult.to_csv(path_or_buf=vPathToSaveTables + ".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
         
def PullSimilarDataIntoOneFile(FileNameList, Result_pathway, ListTable):  
    AllFiles                =     os.listdir(Result_pathway) #get all files in the Result_pathway directory
    for typeTable in ListTable:
        FileType                            =       typeTable
        AllFilesWithParam                   =       [i for i in AllFiles if FileType in i]
        ConcatenatedTables(FileNameList, Result_pathway,AllFilesWithParam, FileType)
            
#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTNucleiPloidyPollen3(aImarisId):
    logging.basicConfig(level=logging.DEBUG, filename= "log[XTNucleiPloidyPollen3].log")
    try:
        #Declare global variables
        global gLasttime
        global vImaris
        global DAPIChannel
#        global VolumePollenResult
#        global VolumeNucleiResult
#        global IntensityPollenResult  
#        global IntensityNucleiResult
#        global PollenIndexList
#        VolumePollenResult          = []
#        VolumeNucleiResult          = []
#        IntensityPollenResult       = []
#        IntensityNucleiResult       = []
#        PollenIndexList             = []
        gLasttime                   = None
        OuputTtile=["PollenArea", "PollenIntensity", "NucleiVolume", "NucleiIntensity", "PollenIDList"]
        logtime('Extension XTNucleiPloidyPollen START')
        FileNameList            =   []
        #            Step1: Connect to Imaris
        #==============================================================================
        vImarisLib			=	ImarisLib.ImarisLib()
        # Get an imaris object with id aImarisId
        vImaris             =   vImarisLib.GetApplication(aImarisId)
        # """BEGIN LOOP OVER ALL IMAGES"""
        # Open File and get filename
        if vImaris is not None :
            logtime('Connected to Imaris')
            FileIndex                       =       1 #This  variable is used to count the number of files analysed
#==============================================================================
#             STEP1: Ask user to set input parameters
#==============================================================================
            root1				        =	        Tk()
            Image_folder	           =	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
            root1.destroy()
            root1				        =	        Tk()
            Result_pathway	           =	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
            root1.destroy()
#            CreateDirectoryToSaveFiles(Result_pathway)
            AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
            logtime('Get all files')
            AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
            TotalNumberFile         =           len(AllFilesToBeProcessed)
            logtime('Start bach processing')
            if TotalNumberFile  >   0:
#                for vFileName in AllFilesToBeProcessed:
##==============================================================================
##                     STEP2: Open one image, segement nuclei and pollen grains and get intensity sum and volume
##                     for each each surface created
##==============================================================================
#                    vFullFileName   =           os.path.join(Image_folder, vFileName)
#                    vImaris.FileOpen(vFullFileName, "")
#                    vImage               =   	vImaris.GetDataSet()
#                    if vImage is not None :
##                        GetImageFeatures(FileIndex, Result_pathway, vFileName, vFullFileName)
##                        FileIndex	               +=  		1
#                        FileNameList.append(vFileName)
                CreateTables(None, AllFilesToBeProcessed, "FileNameList", Result_pathway, "", 0, FileIndex)        
#==============================================================================
#                 STEP3: Plot number of nuclei per ploidy, I need to check quality of images first
#==============================================================================
#                VolumeList   =    [item for sublist in VolumeNucleiResult for item in sublist]
#                IntensityList   =    [item for sublist in IntensityNucleiResult for item in sublist]
#                IntensityVolumeArray=np.array([[x,y] for x, y in zip(IntensityList,VolumeList)])
#                kmeans                      =    KMeans(n_clusters=3, random_state=0).fit(IntensityVolumeArray) #I set the number of clusters to be created to 3 to have an intensity value that doubles in each cluster
#                ClusterIndex                =    list(kmeans.labels_)
#                BarplotFigure(ClusterIndex, IntensityList, Result_pathway)
#==============================================================================
#                 STEP4: Save intensity and volume data into .csv files
#==============================================================================
#                Result=[PollenIndexList, VolumePollenResult, VolumeNucleiResult, IntensityPollenResult, IntensityNucleiResult]
#                NameDataFrame=["PollenID", "AreaPollen", "VolumeNuclei", "IntensityPollen", "IntensityNuclei"]
#                for i, j in zip(Result, NameDataFrame):
#                    dataR       =       pd.DataFrame(i).T
#                    dataR.columns=FileNameList
#                    vPathToSaveTables = os.path.join(Result_pathway, j+".csv")
#                    dataR.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')    
#                if FileIndex>1:
#                    PullSimilarDataIntoOneFile(FileNameList, Result_pathway, OuputTtile)          
            else:
                tkMessageBox.showinfo(title="Alert", message="There is no .ims or .ics file detected in the selected folder.")
                quit()
            logtime('Extension XTNucleiPloidyPollen END')
            print ("All tasks have been completed successfully. \n Resulting files are saved in the folder XTNucleiPloidyPollen_Result")
            raw_input("Press Enter to terminate.")
        else:
            tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
            logtime('Extension XTNucleiPloidyPollen END')
    except:
        logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
