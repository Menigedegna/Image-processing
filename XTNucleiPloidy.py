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
# Creator: Mariamawit S. Ashenafi, UZH
# Published on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTNucleiPloidy" icon="Python" tooltip="NucleiPloidy">
#         <Command>PythonXT::XTNucleiPloidy(%i)</Command>
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
import scipy #TODO: edit out?
from sklearn.cluster import KMeans #both numpy and scipy needs to be imported before importing sklearns 
import matplotlib.pyplot as plt
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
  print curtime.ctime(), '[', str(diff), ']', aTitle
  
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
def Segment_Surface(vImage, ch, vSFW, name,vLCF):
    vROI = None
    vATA = 1
    vATM = 0
    vSFS = ''
    vDCI = ch
    vSurface2 = vImaris.GetImageProcessing().DetectSurfaces(vImage, vROI, vDCI,vSFW, vLCF, vATA, vATM, vSFS)
    vSurface2.SetName(name)
    return vSurface2


#Funtion to get mean intensity inside nucleus surface for all channels selected by the user
def GetNucleusIntensity(surf):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    IntenityList        = 	[float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Mean"]
    IntensityMean       =   IntenityList[DAPIChannel]
    IntenityList        =   [float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Sum"]
    IntensitySum       =   IntenityList[DAPIChannel]
    return IntensityMean, IntensitySum
    
    
#For a surface containing several surface, this function selectes the most dense surface and creates a surface containing only this selected surface
def SelectSurface(vscene, surf, surfType, groupContainer):
    vTimeIndex          =   0
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    volID	            = 	[a for a, x in enumerate(vNames) if x == "Volume"]
    vol		            =	[vValues[x] for x in volID]
    if surfType         ==  "Nucleolus":
        areaID	        = 	[a for a, x in enumerate(vNames) if x == "Area"]
        area		    =	[vValues[x] for x in areaID]
        volumeToAreaRatio         =   [x/y for x,y in zip(vol, area)]
        SelectedID      =   max(xrange(len(volumeToAreaRatio)), key=volumeToAreaRatio.__getitem__)
    else:
        SelectedID      =   max(xrange(len(vol)), key=vol.__getitem__)
    IndexToRemove       =   [x for x in range(len(vol)) if x !=SelectedID]
#    for ind in IndexToRemove:
#        surf.RemoveSurface(ind)
    verts			=	surf.GetVertices(SelectedID)
    vNormals			=	surf.GetNormals(SelectedID)
    faces=surf.GetTriangles(SelectedID)
    vNucleiSurface		=	vImaris.GetFactory().CreateSurfaces()
    vNucleiSurface.AddSurface(verts, faces, vNormals, vTimeIndex)    
    vNucleiSurface.SetName(surfType)
    groupContainer.AddChild(vNucleiSurface, -1)
    vscene.AddChild(groupContainer, -1)
    return vNucleiSurface
    
#This function returns the volume value of surface
def GetVolume(surf):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    VolumeId	           = 	[a for a, x in enumerate(vNames) if x == "Volume"]
    Volume		           =	float(vValues[VolumeId[0]])
    return Volume


# ==============================================================================
# Functions required to ask user set immuno signal spot diameter for all channels
# ==============================================================================
def fetch(entries, root, ParameterOption):
    global SmothingFactor
    global aLocalContrastFilterWidth 
    text = entries.get()
    try:
        x = float(text)
    except ValueError:
        print "You must enter a number in field"
        return
    root.destroy()
    if  ParameterOption=="SmothingSurface":
        SmothingFactor = x  
    else :
        aLocalContrastFilterWidth = x


def makeform(root):
    field = "for nucleus segmentation"
    row = Frame(root)
    lab = Label(row, width=30, text=field, anchor='w')
    ent = Entry(row)
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent.pack(side=RIGHT, expand=YES, fill=X)
    return ent


def AskUserSmoothingFactor(ParameterOption):
    root = Tk()
    label_text = "Set smooth surface detail " if ParameterOption=="SmothingSurface" else "Set value of Local Contrast Filter Width "
    row = Frame(root)
    lab = Label(row, width=30, text=label_text, anchor='w')
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ents = makeform(root)
    root.bind('<Return>', (lambda event, e=ents: fetch(e, root, ParameterOption)))
    b1 = Button(root, fg="darkgreen", text='Submit',
                command=(lambda e=ents: fetch(e, root, ParameterOption)))
    b1.pack(side=LEFT, padx=5, pady=5)
    b2 = Button(root, fg="darkred", text='Quit', command=root.quit)
    b2.pack(side=LEFT, padx=5, pady=5)
    root.mainloop()



#==============================================================================
# Function is required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTNucleiPloidy_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)
        
#==============================================================================
# This function:
# Segments DAPI chanenl into nucleus surface 
# Saves the mean and sum DAPI intensity inside nucleus and the volume of the nucleus
#==============================================================================
def GetImageFeatures(FileIndex, Result_pathway, vFileName):
    global vImaris
    global DAPIChannel
    global SmothingFactor
    global aLocalContrastFilterWidth
    NucleusSurfaceVolume = None
    IntensityMean        = None
    IntensitySum         = None
    vImage               =   	vImaris.GetDataSet()
    vImaris.GetSurpassCamera().Fit() #scene fits the view of the surpass camera object 
    if vImage is not None:
        numberIndex        =    vImage.GetSizeC()
        if FileIndex       ==   1:
            DAPIChannel    =    Ask_user(numberIndex)
            AskUserSmoothingFactor("SmothingSurface")
            AskUserSmoothingFactor("LocalContrastFilterWidth")
        NucleusSurface     =    None
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        vImage                     =   vImaris.GetDataSet()
        logtime('Segmentation START - image_'+str(FileIndex))
        vScene                     =   vImaris.GetSurpassScene()
        vFactory		              =	vImaris.GetFactory().CreateFactory()
        GroupOfObjects	           =	vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
        try: #Select biggest surface
            SurfP=Segment_Surface(vImage, DAPIChannel, SmothingFactor, "Rougher Nucleus segmentation",aLocalContrastFilterWidth)
            if SurfP.GetNumberOfSurfaces()>0:
                NucleusSurface=SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
                logtime('Nucleus surface segmented - Image '+str(FileIndex))
                NucleusSurfaceVolume    =   GetVolume(NucleusSurface)
                IntensityMean, IntensitySum     =    GetNucleusIntensity(NucleusSurface)
#                vPathToSaveTables = os.path.join(Result_pathway, "Snapshot_"+vFileName+".tif")
#                vImaris.SaveSnapShot(vPathToSaveTables)
            else:
                print "No nucleus surface detected for image: "+vFileName
        except:#Select biggest surface
            logging.exception("Select biggest surface:")
            print "There was a problem a segmentation problem with "+vFileName +". Please check log file."
            pass
#        vPath = os.path.join(Result_pathway, vFileName+".ims")
#        vImaris.FileSave(vPath, "") #to save modified file
        vScene.RemoveChild(GroupOfObjects) #remove containter from scene
    return vImage is None, NucleusSurfaceVolume, IntensityMean, IntensitySum

#==============================================================================
# This function:
# Do regression fit of DAPI sum intensity vs nucleus volume to cluster nuclei into different ploidy level
# Check that the average intensity per cluster doubles in each higher ploidy
# Create boxplot of DAPI sum intensity as a function of cluster ID and save in tif file
# Create barplot of number of nuclei and average and save in tif file
# Saves sum and mean DAPI intensity, nucleus volume, cluster ID and filename in csv file
#==============================================================================
def ClusterNuclei(VolumeList,SumIntensityList, FileNameList,  MeanIntensityList, Result_pathway):
    IntensityVolumeArray=np.array([[x,y] for x, y in zip(SumIntensityList,VolumeList)])
    ClusterNotReached = True
    n_clusters=2 #starting number of cluster
    NumberOfRun=0 # the while loop will run at most 3 times
    while ClusterNotReached and NumberOfRun<=4: # run loop until average intensity between clusters are doubled
        kmeans                      =    KMeans(n_clusters=n_clusters, random_state=0).fit(IntensityVolumeArray) #I set the number of clusters to be created to 3 to have an intensity value that doubles in each cluster
        ClusterIndex                =    list(kmeans.labels_)
        UnikCluster                 =    set(ClusterIndex)
        IntensityGroup                        =    []
        for clusterID in UnikCluster:
            IntensityGroup.append(np.array([x for x,y in zip(SumIntensityList, ClusterIndex) if y==clusterID]))
        AverageIntensity=[np.mean(x) for x in IntensityGroup]
        AverageIntensity, UnikCluster= zip(*sorted(zip(AverageIntensity, UnikCluster)))
        oldNumberOfCluster=n_clusters
        for i in range(n_clusters-1):
            #if average intensity the next 2 clusters is less than the double of cluster x, lower the number of cluster
            if i+2<=n_clusters-1 and np.mean(AverageIntensity[i+1: i+3])<=2*AverageIntensity[i]:
                n_clusters=n_clusters-1 
                break  
            #if average intensity triples in next cluster, increase the number of cluster
            if AverageIntensity[i+1] >=3*AverageIntensity[i]:
                n_clusters=n_clusters+1 
                break
        if  n_clusters==oldNumberOfCluster: #non of the conditions above are true = the clustering is correct
            ClusterNotReached=False
        NumberOfRun=NumberOfRun+1
    #==============================================================================
    # Substitute cluster index with P0, P1, P2 etc
    #==============================================================================
    SubstitueClusterIndex=["P"+str(x) for x in range(n_clusters)]  
    for y,z in zip(SubstitueClusterIndex, UnikCluster):
        ClusterIndex=[y if x==z else x for x in ClusterIndex] 
    #==============================================================================
    # Create and save boxplot intensity of nuclei per cluster
    #==============================================================================
    IntensityGroup                        =    []
    for clusterID in SubstitueClusterIndex:
        IntensityGroup.append(np.array([x for x,y in zip(SumIntensityList, ClusterIndex) if y==clusterID]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot( IntensityGroup, notch=True, patch_artist=True)
    ax.set_xticklabels(SubstitueClusterIndex)
    ax.set_ylabel('Average DAPI intensity sum', color="darkblue")
    ax.set_xlabel('Ploidy Level', color="darkblue")
    ax.set_title('DAPI sum intensity by cluster of nuclei', color="darkblue")
    PlotFile = os.path.join(Result_pathway, "PloidyIntensityPlot.svg") 
    fig.savefig(PlotFile, format='svg', dpi=300)
    plt.close(fig)
    #==============================================================================
    # Create and save barplot number of nuclei per cluster
    #==============================================================================
    NumberOfNuclei=[len(x) for x in IntensityGroup]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(n_clusters), NumberOfNuclei, align='center')
    ax.set_xticks(range(n_clusters))
    ax.set_xticklabels(SubstitueClusterIndex)
    ax.set_ylabel('Number of nuclei', color="darkblue")
    ax.set_xlabel('Ploidy Level', color="darkblue")
    ax.set_title('Number of nuclei by cluster of nuclei', color="darkblue")
    PlotFile = os.path.join(Result_pathway, "PloidyNumberOfNucleiPlot.svg") 
    fig.savefig(PlotFile, format='svg', dpi=300)
    plt.close(fig)
    #==============================================================================
    # Save table col (filename, sum intensity,mean intensity, volume, clusterid), row(number of .ims file) in .csv file     
    #==============================================================================
    dataR       =       pd.DataFrame({"FileName": FileNameList, "Volume": VolumeList, "MeanIntensity": MeanIntensityList, "SumIntensity":  SumIntensityList, 'ClusterId': ClusterIndex})
    vPathToSaveTables = os.path.join(Result_pathway, "XTNucleiPloidy_Result.csv")
    dataR.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')


#==============================================================================
# Main function:
# Connects to Imaris 
# Open each image in folder set by user
# Segment nucleus and gets DAPI sum intensity and volume of nucleus
# Do a regression fit of DAPI sum intensity vs volume of nucleus  and cluster nuclei
# Create and save .csv table with columns= filename, DAPI sum intensity, volume, Cluster ID; rows= number of .ims files with DAPI channel inside folder
#==============================================================================
def XTNucleiPloidy(aImarisId):
    logging.basicConfig(level=logging.DEBUG, filename= "log[XTNucleiPloidy].log")
    #Declare global variables
    global gLasttime
    global vImaris
    gLasttime               =    None
    logtime('Extension XTNucleiPloidy START')
    FileNameList            =   []
    VolumeList               =   []
    MeanIntensityList               =   []
    SumIntensityList               =   []
    #==============================================================================
    #            STEP1 : CONNECT TO IMARIS
    #==============================================================================
    try: #Is there problem to connect to Imaris?
        vImarisLib			=	ImarisLib.ImarisLib()
        # Get an imaris object with id aImarisId
        vImaris             =   vImarisLib.GetApplication(aImarisId)
    except:#Is there problem to connect to Imaris?
        logging.exception("Imaris connection:")
        print "There was a problem to connect the to Imaris. Please check log file."
        pass
    if vImaris is not None : #Imaris is open?
        logtime('Connected to Imaris')
        vImaris.GetSurpassCamera().SetPerspective(0) #Set camera to orthographic view 
        FileIndex                       =       1 #This  variable is used to count the number of files analysed
        #==============================================================================
        # 		STEP 2 : ASK USER TO SET INPUT PARAMETERS
        #==============================================================================
        try: #Is there a problem to open the dialog window?  
            root1				       =	        Tk()
            Image_folder	       =	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
            root1.destroy()
        except:#Is there a problem to open the dialog window? 
            logging.exception("tkFileDialog issue:")
            print "There was a problem to open window to ask user to set input parameters. Please check log file."
            pass
        Result_pathway          =           os.path.join(Image_folder, "XTNucleiPloidy_Result")
        CreateDirectoryToSaveFiles(Result_pathway)
        AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
        logtime('Get all files')
        AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
        TotalNumberFile         =           len(AllFilesToBeProcessed)
        logtime('Start bach processing')
        NumberFISHList                   =       []
        #==============================================================================
        # 		STEP 3 : OPEN IMAGE FILES AND GET INTENSITY AND VOLUME OF NUCLEI
        #==============================================================================
        if TotalNumberFile  >   0: # Are there any .ims files?
            for vFileName in AllFilesToBeProcessed: #loop over all .ims or .ics files
                try: #Is there a problem inside one file (no image or segmentation issue)?
                    vFullFileName   =           os.path.join(Image_folder, vFileName)
                    vImaris.FileOpen(vFullFileName, "")
                    ImageIsEmpty,NucleusSurfaceVolume, IntensityMean, IntensitySum    =           GetImageFeatures(FileIndex, Result_pathway, vFileName)
                except: #Is there a problem inside one file (no image or segmentation issue)?
                    logging.exception("Image="+vFileName+" :")
                    print "There is an issue with file : "+ vFileName
                    continue #it will not treat this image and move on to the next one
                if not ImageIsEmpty and  NucleusSurfaceVolume is not None :
                    FileIndex	               +=  		1
                    FileNameList.append(vFileName)
                    VolumeList.append(NucleusSurfaceVolume)
                    MeanIntensityList.append(IntensityMean)
                    SumIntensityList.append(IntensitySum)
            #==============================================================================
            # 		STEP 4 : CLUSTER NUCLEI ACCORDING TO THEIR PLOIDY
            #==============================================================================
            try:    #Is there a clustering problem or plotting
                ClusterNuclei(VolumeList,SumIntensityList, FileNameList,  MeanIntensityList, Result_pathway)
            except:  #Is there a clustering problem
                logging.exception("Nuclei clustering issue:")
                print "There is a nuclei clustering issue. Please check log file."
                pass
        else: # Are there any .ims files?
            tkMessageBox.showinfo(title="Alert", message="There is no .ims or .ics file detected in the selected folder.")
            quit()
        logtime('Extension XTNucleiPloidy END')
        print "Resulting files are saved in the folder XTNucleiPloidy_Result"
        print "Log file is saved in the same folder containing this plugin"
        raw_input("Press Enter to terminate.")
#==============================================================================
# End of extension
#==============================================================================
