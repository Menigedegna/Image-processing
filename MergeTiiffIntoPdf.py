# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 03:37:52 2017

@author: Pheonix
"""
import sys
import os
sys.path.insert(0, "H:\Python") #I need to add this to import futures on the VMs
from fpdf import FPDF
pdf = FPDF()
# imagelist is the list with all image filenames
x=y=0
w=h=2
count=0
FolderInitial="J:\06072017\ImsFiles\s16_Wt2\1\XTSegmentNuclei_Result"
AllFilesInDirectory     =           os.listdir(FolderInitial) #get all files in the Image_folder directory
imagelist   =           [i for i in AllFilesInDirectory if i.endswith('.tif')] #select files with .ims or .ics extensions
for image in imagelist:
    pdf.add_page()
    pdf.image(image,x,y,w,h)
    x+=w
    if count<6:
        count=count+1  
    else:
        count=0
        y+=h 
        x=0
OutputFile=os.path.join(FolderInitial,"SegmentSnapshots.pdf")
pdf.output(OutputFile, "F")