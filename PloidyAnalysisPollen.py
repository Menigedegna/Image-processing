#==============================================================================
# Objectives of this Python script :
#   Read output data from XTSegmentNucleus Imaris plugin,
#   Analyse ploidy
# Creator: Mariamawit S. Ashenafi, Celia Baroux, UZH
# Published on 23.01.2016
#==============================================================================
#python library
import numpy as np
import pandas as pd
import itertools
import scipy #TODO: edit out?
from sklearn import linear_model as lm
from sklearn.cluster import KMeans #both numpy and scipy needs to be imported before importing sklearns 
import matplotlib.pyplot as plt
import os
from sklearn.metrics import r2_score

def ScatterPlot(param1, param2, Result_pathway, spotType, UnikCluster, Param2, Param1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    col=["gold", "orangered", "olive", "mediumslateblue"]
    for i,x,y in zip(range(len(param1)), param1, param2):
        if Param1   != "Predicted" and Param1   != "Within-clusters sum of squares":
            x=[np.log(a) for a in x]
            y=[np.log(a) for a in y]
        if Param1   != "Within-clusters sum of squares" :            
            fit = np.polyfit(x, y, deg=1) #calculate regression
            p = np.poly1d(fit) #get r squared 
            coefficient_of_dermination = r2_score(y, p(x))
            fittedY     =   [fit[0] * a + fit[1] for a in x] 
            ax.plot(x, fittedY, color=col[i])
            Lab1        =   UnikCluster[i]+" : y= "+str(p)+", R2="+str(round(coefficient_of_dermination, 2))
        if Param1   == "Within-clusters sum of squares" : 
            Lab1        =  " "
            ax.plot(x, y, c=col[i], zorder=2)    
            ax.grid(True)
        ax.scatter(x, y, s=10, c=col[i], marker="X", label=Lab1)
    # axes and labels
    var1="" if Param1=="Predicted" or Param1   == "Within-clusters sum of squares" else ' (um3))'
    var2="" if Param1=="Predicted" or Param1   == "Within-clusters sum of squares" else 'log('
    var3="" if Param1=="Predicted" or Param1   == "Within-clusters sum of squares" else ')'
    ax.set_ylabel(var2+Param1+var3, color="white")
    ax.set_xlabel(var2+Param2+var1+var3, color="white")
    ax.tick_params(axis='y', colors='white')
    ax.tick_params(axis='x', colors='white')
    ax.set_title(Param1+" as a function of "+Param2+", "+spotType, color="white")
    # Put a legend below current axis
    lgd = ax.legend(loc='center left', fancybox=True, shadow=True, bbox_to_anchor=(1, 0.5), fontsize=8)
    PlotFile = os.path.join(Result_pathway, Param2+Param1+spotType+".png") 
    fig.savefig(PlotFile, bbox_extra_artists=(lgd,), bbox_inches='tight', transparent=True, dpi=300)
    plt.close(fig)
    
def BarplotFigure(ClusterIndex, SumIntensityList,Result_pathway):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ## the data
    UniqueValuesClusters        = list(set(ClusterIndex))
    NumberOfImages              = len(SumIntensityList)
    NumberOfNucleiPerCluster    = [(float(ClusterIndex.count(n)) / float(NumberOfImages)) * 100.0 for n in
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
    rects2 = ax.bar(ind + width, NumberOfNucleiPerCluster, width,
                    color='mediumslateblue')
    # axes and labels
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0.0, 100.0)
    ax.set_ylabel('Number of nuclei', color="mediumslateblue")
    ax.tick_params(axis='y', colors='mediumslateblue')
    ax.tick_params(axis='x', colors='white')
    ax2.set_ylabel('Average DAPI intensity sum', color="olive")
    ax2.tick_params(axis='y', colors='olive')
    ax.set_title('DAPI intensity sum and number of nuclei by cluster', color="white")
#    xTickMarks = ['Cluster' + str(i) for i in range(1, len(UniqueValuesClusters) + 1)]
    xTickMarks = UniqueValuesClusters
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ## add a legend
    ax.legend((rects1[0], rects2[0]), ('Intensity sum average', 'Number of nuclei'), fontsize=6)
#    lgd = ax.legend((rects1[0], rects2[0]), ('Intensity sum average', 'Number of nuclei'), loc='upper center', ncol=1, fancybox=True, shadow=True, bbox_to_anchor=(1, 0.5), fontsize=6)    
    PlotFile = os.path.join(Result_pathway, "PloidyPlot.png") 
    fig.savefig(PlotFile, transparent=True, dpi=300)
    plt.close(fig)
    dataR       =       pd.DataFrame({"Cluster": ind, "SumIntensity Average":  AverageIntensity, 'SumIntensity STD': IntensityStd, 'Number of Nuclei': NumberOfNucleiPerCluster})
    vPathToSaveTables = os.path.join(Result_pathway,"ClusterIntensityAverage.csv")
    dataR.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    return UniqueValuesClusters

def ReadTable(FileName, colTitle, nextCol, NumberCol):
    FileName0                                   =       os.path.join(r"D:\Workspace Aptana\ImarisPlugins\St0", FileName)
    FileName1                                   =       os.path.join(r"D:\Workspace Aptana\ImarisPlugins\St1", FileName) 
    Param                                       =       pd.read_csv(FileName0, sep=",", header='infer',decimal='.')
    Param0                                       =       list(Param[colTitle])
    Param2                                      =       pd.read_csv(FileName1, sep=",", header='infer',decimal='.')
    Param2                                      =       list(Param2[colTitle])[0:2]
    Param2.extend(Param0[22:30])
    Param2.extend(Param0[0:5])
    if NumberCol>1:
         Param1                                      =       list(Param[nextCol])
         Param2.extend(Param1[5:15])
    else: 
         Param2.extend(Param0[5:15])
    return Param2

def PloidyTest(Result_pathway):
    DapiSumIntensity                            =       ReadTable("DAPIntensitySum.csv", "0",None, 1)
    Volume                                      =       ReadTable("NucleusVolume.csv", "0", None, 1)
    #remove outliers, after running this function I found that there is one value with very low DAPI intensity
    MinDAPIiD                                   =       DapiSumIntensity.index(min(DapiSumIntensity))
    del DapiSumIntensity[MinDAPIiD]; del Volume[MinDAPIiD]
    #Ploidy analysis will be based on Sum of DAPI intensity and nucleus volume fit
    IntensityVolumeArray                        =       np.array([[x,y] for x, y in zip(DapiSumIntensity,Volume)])
    NumberCluster               =    10
    SDDList                     =       []
    for i in range(1,NumberCluster):
        kmeans                      =    KMeans(n_clusters=i, random_state=0).fit(IntensityVolumeArray) #I 
        SDD                         =    float(kmeans.inertia_)
        SDDList.append(SDD)
    ScatterPlot([range(1,NumberCluster)],[SDDList],  Result_pathway, "", [""],  "Number of Cluster", "Within-clusters sum of squares")  
    
def AnalysePloidy(Result_pathway, NumberCluster):
    DapiSumIntensity                            =       ReadTable("DAPIntensitySum.csv", "0",None, 1)
    Volume                                      =       ReadTable("NucleusVolume.csv", "0", None, 1)
    #remove outliers, after running this function I found that there is one value with very low DAPI intensity
    MinDAPIiD                                   =       DapiSumIntensity.index(min(DapiSumIntensity))
#    del DapiSumIntensity[MinDAPIiD]; del Volume[MinDAPIiD]
    #Ploidy analysis will be based on Sum of DAPI intensity and nucleus volume fit
    IntensityVolumeArray                        =       np.array([[x,y] for x, y in zip(DapiSumIntensity,Volume)])
    kmeans                      =    KMeans(n_clusters=NumberCluster, random_state=0).fit(IntensityVolumeArray) #I set the number of clusters to be created to 3 to have an intensity value that doubles in each cluster
    ClusterIndex                =    list(kmeans.labels_)    
    dataR       =       pd.DataFrame({"FileIndex": range(len(Volume)), "Volume": Volume, "SumIntensity":  DapiSumIntensity, 'ClusterId': ClusterIndex})
    vPathToSaveTables = os.path.join(Result_pathway, "XTNucleiPloidy_Result.csv")
    dataR.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    UniqueClusterID             =   BarplotFigure(ClusterIndex, DapiSumIntensity, Result_pathway)
    ScatterPlot([Volume], [DapiSumIntensity], Result_pathway, "", [""], "Volume", "DAPIIntensity")
    #Rename Cluster index
    for j in range(NumberCluster):
        ClusterIndex            =   [UniqueClusterID[j] if x==j else x for x in ClusterIndex]
    UniqueValuesClusters                      =        list(set(ClusterIndex))
    SelectedValuesNS            =   PloidySegregate(DapiSumIntensity, UniqueValuesClusters, ClusterIndex)
    SelectedValuesV            =   PloidySegregate(Volume, UniqueValuesClusters, ClusterIndex)
    ScatterPlot(SelectedValuesV, SelectedValuesNS, Result_pathway, "Ploidy", UniqueValuesClusters, "Volume", "DAPIIntensity")   
    return ClusterIndex, MinDAPIiD

def PloidySegregate(ListSel, UniqueValuesClusters, ClusterIndex):
    NumberOfImages              =   len(ListSel)
    SelectedValuesNS            =   []   
    Param                       =   ListSel
    for clusterID in UniqueValuesClusters:
        SelectedValuesInd = [Param[x] for x in range(NumberOfImages) if ClusterIndex[x] == clusterID]
        SelectedValuesNS.append(SelectedValuesInd)  
    return SelectedValuesNS

def BarplotFigureSimpler(Param1,Param2, ParamLabel,Result_pathway, param2):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ## the data
    AverageNumberSpotSer2           = [round(float(np.mean(n)), 2) for n in Param1]
    STDNumberSpotSer2               = [round(float(np.std(n)), 2) for n in Param1]
    AverageNumberSpotSer5           = [round(float(np.mean(n)), 2) for n in Param2]
    STDNumberSpotSer5               = [round(float(np.std(n)), 2) for n in Param2]
    ind = np.array(range(len(ParamLabel)))  # the x locations for the groups
    width = 0.35  # the width of the bars
    ## the bars
    rects1 = ax.bar(ind, AverageNumberSpotSer2, width,
                     color='olive', 
                     yerr=STDNumberSpotSer2,
                     error_kw=dict(elinewidth=2, ecolor='olive'))
    rects2 = ax.bar(ind + width, AverageNumberSpotSer5, width,
                     color='mediumslateblue',
                     yerr=STDNumberSpotSer5,
                     error_kw=dict(elinewidth=2, ecolor='mediumslateblue'))
    # axes and labels
    ax.set_xlim(-width, len(ind) + width)
#    ax.set_ylim(0.0, 100.0)
    ax.set_ylabel('Number of spots', color="white")
    ax.tick_params(axis='y', colors='white')
    ax.tick_params(axis='x', colors='white')
    ax.set_title('Number of spots by '+param2, color="white")
#    xTickMarks = ['Cluster' + str(i) for i in range(1, len(UniqueValuesClusters) + 1)]
    xTickMarks = ParamLabel
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    ax.legend((rects1[0], rects2[0]), ('Ser2P', 'Ser5P'), fontsize=8)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ## add a legend
    PlotFile = os.path.join(Result_pathway, "NumberSpotBy"+param2+".png") 
    fig.savefig(PlotFile, transparent=True)
    plt.close(fig)

def PredictValues(param1,param2, param3):
    RefFit=np.array([[np.log(u),np.log(v)] for u,v in zip(param1,param2)])
    clf = lm.LinearRegression(fit_intercept=False)
    clf.fit(RefFit, param3)
    predValue   =   list(np.dot(RefFit, clf.coef_)) 
    print clf.coef_
    return predValue

def PlotPrediction(param1,param2, param3, ClusterIndex, SelectedValuesNS, Result_pathway, LabS):
        UniqueValuesClusters                      =        list(set(ClusterIndex))
        predValue   =   PredictValues(param1,param2, param3) #Multivariate regression model
        ScatterPlot([param3], [predValue], Result_pathway, LabS, [""], "Observed", "Predicted")
        selPred     =     PloidySegregate(predValue, UniqueValuesClusters, ClusterIndex)
        ScatterPlot(SelectedValuesNS, selPred, Result_pathway, LabS+"_Ploidy", UniqueValuesClusters, "Observed", "Predicted") 
        Par=[np.log(x) for x in param2]
        SelectedDAPI     =     PloidySegregate(Par, UniqueValuesClusters, ClusterIndex )
        ScatterPlot(SelectedDAPI, selPred, Result_pathway, LabS+"_Ploidy", UniqueValuesClusters, "DAPIIntensity", "Predicted") 
      
def NumberOfSpotStudy(Result_pathway, ClusterIndex, RemoveID):
    Ser2PSpot                                 =       ReadTable("NumberOfSpots.csv" ,"0", "1", 2)
    Ser5PSpot                                 =       ReadTable("NumberOfSpots.csv" ,"1", "0", 2)
    #    del Ser2PSpot[RemoveID]; del Ser5PSpot[RemoveID]; del Volume[RemoveID]
    UniqueValuesClusters                      =        list(set(ClusterIndex))
    ListSpots                                 =     [Ser2PSpot, Ser5PSpot]
    LabSpots                                  =     ["Ser2P", "Ser5P"]
    res                                       =     []   
    res2                                       =     []   
    for i in range(len(ListSpots)):    
        SelectedValuesNS     =     PloidySegregate(ListSpots[i], UniqueValuesClusters, ClusterIndex)
        res.append(SelectedValuesNS)
        x = [ListSpots[i][0:11], ListSpots[i][11:-1]]
        res2.append(x)
        y=["3H1", "Wt"]
        a= ["NucleusVolume.csv", "DAPIntensitySum.csv"]
        b= ["Volume", "DAPIIntensity"]
        res3=[]
        for param2 in range(2):
            tab                                    =       ReadTable(a[param2], "0", None, 1)
#            del tab[RemoveID]
            res3.append(tab)
            SelectedValuesV     =     PloidySegregate(tab, UniqueValuesClusters, ClusterIndex )
            ScatterPlot([tab], [ListSpots[i]], Result_pathway, LabSpots[i], [""], b[param2], "NumberOfSpot")
            ScatterPlot(SelectedValuesV, SelectedValuesNS, Result_pathway, LabSpots[i]+"_Ploidy", UniqueValuesClusters, b[param2], "NumberOfSpot")            
        PlotPrediction(res3[0],res3[1], ListSpots[i], ClusterIndex,  SelectedValuesNS, Result_pathway, LabSpots[i])  
    #Bar plot number spot 
    BarplotFigureSimpler(res[0], res[1], UniqueValuesClusters,Result_pathway, "Ploidy")#by ploidy cluster
    BarplotFigureSimpler(res2[0], res2[1], y,Result_pathway, "PlantType")#by plant type cluster

           




Result_pathway                      =       r"D:\Workspace Aptana\ImarisPlugins\res"
PloidyTest(Result_pathway)
NumberCluster                       =       4
ClusterIndex, RemoveID                              =       AnalysePloidy(Result_pathway, NumberCluster)
NumberOfSpotStudy(Result_pathway, ClusterIndex, RemoveID)

#OtherCh=[1,0]
#selectedChannel=2
##def ClusterAnalyse(selectedChannel):
#for sp in range(selectedChannel):
#    tab0                               =       pd.read_csv("D:\Workspace Aptana\ImarisPlugins\St0\SpotsPositions_Ch"+str(OtherCh[sp])+".csv", sep=";", header='infer',decimal='.')
#    tab1                               =       pd.read_csv("D:\Workspace Aptana\ImarisPlugins\St0\SpotsPositions_Ch"+str(sp)+".csv", sep=";", header='infer',decimal='.')
#    tab2                               =       pd.read_csv("D:\Workspace Aptana\ImarisPlugins\St1\SpotsPositions_Ch"+str(sp)+".csv", sep=";", header='infer',decimal='.')
#    res             =                   pd.concat(tab0, tab1[5*3:15*3], tab1[21*3:90], tab2)
#    resList=[]
#    for x range(3):
#        resx            =           res[range(x,90,3)]
#        resList.append(resx)
#    for imageId in range(30):
#        NumberSpot=len(resx)
#        listProd        =   []
#        for x in range(3)
#            repx= resList[x][imageId]
#            repx=[list(repx)]*NumberSpot
#            repx0=pd.DataFrame(repx)
#            repx1=pd.DataFrame(repx).T
#            res=    repx1.subtract(repx0)
#            res=    np.square(res)
#            listProd.append(res)
#        SumProduct=listProd[0].add(listProd[1])
#        SumProduct=SumProduct.add(listProd[2])
#        R =SumProduct.apply(np.sqrt)
