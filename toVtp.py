import vtkXMLFile
import numpy as np
import os
print("convert to .vtp format")
rltPath=os.path.split(os.path.realpath(__file__))[0]
insegment0=rltPath+"\\output\\tradition_centSize0.txt"
insegment1=rltPath+"\\output\\raytracing_centSize0.txt"
outsegment0=rltPath+"\\output\\tradition.vtp"
outsegment1=rltPath+"\\output\\raytracing.vtp"

data0 = np.loadtxt(open(insegment0, "rb"), delimiter=",", skiprows=0, ndmin=2)
data1 = np.loadtxt(open(insegment1, "rb"), delimiter=",", skiprows=0, ndmin=2)
vtkXMLFile.vtkXMLVoxel.to_vtp(data0,outsegment0,0,"segment")
vtkXMLFile.vtkXMLVoxel.to_vtp(data1,outsegment1,0,"segment")

voxcntFile=rltPath+"\\input\\voxelModel-CentSize.txt"
outvoxcntFile=rltPath+"\\output\\voxelModel.vtp"
vxlmdl = np.loadtxt(open(voxcntFile, "rb"), delimiter=",", skiprows=0, ndmin=2)
vtkXMLFile.vtkXMLVoxel.to_vtp(vxlmdl,outvoxcntFile,5,"voxelmodel")

coordFile=rltPath+"\\input\\rayCoordinates10.txt"
rayOutFile=rltPath+"\\output\\ray.vtp"
newdata=[]
raydata = np.loadtxt(open(coordFile, "rb"), delimiter=",", skiprows=0, ndmin=2)
for i in range(0,len(raydata)):
    a=list(raydata[i])
    a.append(0)
    newdata.append(a)
rayline=vtkXMLFile.vtkXMLLine(newdata)
rayline.to_vtk_xml(rayOutFile)
print("finshing conversion")
