from vtk import *
import numpy as np
class vtkXMLLine:
    def __init__(self,ptsVal,valName="val"):
        self.lines=ptsVal
        self.valName=valName
        pass

    def to_vtk_xml(self, outfile):
        if len(self.lines)<1:
            return
        # https://lorensen.github.io/VTKExamples/site/Python/UnstructuredGrid/UGrid/
        numOfComp=len(self.lines[0])-6
        # setup points and vertices
        Ray_pts = vtk.vtkPoints()
        Ray_lines = vtk.vtkCellArray()
        vals = vtk.vtkDoubleArray()
        vals.SetNumberOfComponents(numOfComp)
        vals.SetName(self.valName)
        for i in range(0,len(self.lines)):
            line=self.lines[i]
            # id = i  # validRayids[i]
            Ray_pts.InsertNextPoint(line[0], line[1], line[2])
            Ray_pts.InsertNextPoint(line[3], line[4], line[5])
            aline = vtk.vtkLine()
            aline.GetPointIds().SetId(0, i * 2)
            aline.GetPointIds().SetId(1, i * 2+1)
            Ray_lines.InsertNextCell(aline)
            for val in line[6:]:
                vals.InsertNextValue(val)

        Linedata = vtk.vtkPolyData()
        Linedata.SetPoints(Ray_pts)
        Linedata.SetPolys(Ray_lines)
        Linedata.GetCellData().SetScalars(vals)
        Linedata.Modified()
        writer1 = vtk.vtkXMLPolyDataWriter()
        writer1.SetFileName(outfile)
        writer1.SetDataModeToAscii()
        writer1.SetInputData(Linedata)
        # writer.SetInputConnection(ugReader.GetOutputPort())
        writer1.Write()
class vtkXMLVoxel:
    def __init__(self,valName="val"):
        '''
        @:param voxels, cent and size
        '''

        pass
    @classmethod
    def to_vtu(cls,lonCoords,latCoords,rCoords,vals,outfile,valName="val",beVector=False):
        '''
        :param lonCoords: list数组，经度方向坐标，弧度制
        :param latCoords: list数组，纬度方向坐标，弧度制
        :param rCoords: list数组，径向方向坐标，米
        :param vals: 4维的np数组，第0~3维依次为R,Lat,Lon,属性维
        :param outfile: 输出文件名，以vtu为后缀
        :param quality: 每个格网按2^quality进行加密
        :param valName: 输出变量的名称
        :return:
        '''
        ugrid=cls.__GenUnStructHexahedronCell(lonCoords,latCoords,rCoords,vals,valName,beVector)
        if ugrid!=None:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(outfile)
            writer.SetDataModeToAscii()

            if vtk.VTK_MAJOR_VERSION <= 5:
                writer.SetInput(ugrid)
            else:
                writer.SetInputData(ugrid)
            writer.Write()
            return True
        else:
            return False
    @classmethod
    def to_vtp(cls, centSizesVal,outfile, quality=1,valName="val"):
        voxs=centSizesVal
        if len(voxs)<1:
            return
        import math
        facepts = []
        faces = []
        for vox in voxs:
            minR = vox[2]-vox[5]/2
            maxR = vox[2]+vox[5]/2
            minLat = vox[1]-vox[4]/2
            maxLat = vox[1]+vox[4]/2
            minLon = vox[0]-vox[3]/2
            maxLon = vox[0]+vox[3]/2
            numOfSection = math.pow(2, quality)
            cls.__gen_llr_voxel_surface(minLon, maxLon, minLat, maxLat, minR, maxR, numOfSection, facepts,
                                                 faces, vox[6:])
                # setup points and vertices
        quadpts = vtk.vtkPoints()
        quadfaces = vtk.vtkCellArray()
        for pt in facepts:
            x, y, z = cls.__sphere_llr_geocent(pt[0], pt[1], pt[2])
            quadpts.InsertNextPoint(x, y, z)
        numOfComp=len(voxs[0])-6
        vals = vtk.vtkDoubleArray()
        vals.SetNumberOfComponents(numOfComp)
        vals.SetName(valName)
        for face in faces:
            aquad = vtk.vtkQuad()
            aquad.GetPointIds().SetId(0, face[0])
            aquad.GetPointIds().SetId(1, face[1])
            aquad.GetPointIds().SetId(2, face[2])
            aquad.GetPointIds().SetId(3, face[3])
            quadfaces.InsertNextCell(aquad)
            vals.InsertNextTypedTuple(face[4])



        polydata = vtk.vtkPolyData()
        polydata.SetPoints(quadpts)
        polydata.SetPolys(quadfaces)
        polydata.GetCellData().SetScalars(vals)
        polydata.Modified()
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(outfile)
        writer.SetDataModeToAscii()
        # writer.SetDataModeToBinary()
        writer.SetInputData(polydata)
        # writer.SetInputData(Linedata)
        # writer.SetInputConnection(ugReader.GetOutputPort())
        writer.Write()

    @classmethod
    def __GenStructedCells(cls, lonCoords, latCoords, rCoords, vals):
        '''
        :param lonCoords: list对象
        :param latCoords:list对象
        :param rCoords:
        :param vals: list对象，每个元素是一个truple
        :param cells:
        :param pts:
        :param scals:
        :return:
        '''
        pass
        '''
        import vtk
        lonNum = np.shape(lonCoords)[0]
        latNum = np.shape(latCoords)[0]
        rNum = np.shape(rCoords)[0]
        valNum = len(vals)
        if valNum != 1 and valNum != (lonNum - 1) * (latNum - 1) * (rNum - 1):
            return 0
        pts = vtk.vtkPoints()
        sgrid = vtk.vtkStructuredGrid()
        sgrid.SetDimensions()
        scals = vtk.vtkFloatArray()
        for r in rCoords:
            for lat in latCoords:
                for lon in lonCoords:
                    pt = cls.__sphere_llr_geocent(lon, lat, r)
                    pts.InsertNextPoint(pt)
        latLonNum = lonNum * latNum
        count = 0
        for iR in range(0, rNum - 1):
            for iLat in range(0, latNum - 1):
                for iLon in range(0, lonNum - 1):
                    hex = vtk.vtkHexahedron()
                    hex.GetPointIds().SetId(0, iR * latLonNum + iLat * lonNum + iLon + ptstartid)
                    hex.GetPointIds().SetId(1, iR * latLonNum + iLat * lonNum + iLon + 1 + ptstartid)
                    hex.GetPointIds().SetId(2, iR * latLonNum + (iLat + 1) * lonNum + iLon + 1 + ptstartid)
                    hex.GetPointIds().SetId(3, iR * latLonNum + (iLat + 1) * lonNum + iLon + ptstartid)
                    hex.GetPointIds().SetId(4, (iR + 1) * latLonNum + iLat * lonNum + iLon + ptstartid)
                    hex.GetPointIds().SetId(5, (iR + 1) * latLonNum + iLat * lonNum + iLon + 1 + ptstartid)
                    hex.GetPointIds().SetId(6, (iR + 1) * latLonNum + (iLat + 1) * lonNum + iLon + 1 + ptstartid)
                    hex.GetPointIds().SetId(7, (iR + 1) * latLonNum + (iLat + 1) * lonNum + iLon + ptstartid)
                    sgrid.InsertNextCell(hex)
                    if valNum > 1:
                        scals.InsertNextTypedTuple(vals[count])
                    else:
                        scals.InsertNextTypedTuple(vals[0])
                    count = count + 1
        return count
        '''
    @classmethod
    def __GenUnStructHexahedronCell(cls,lonCoords,latCoords,rCoords,vals,valName="val",beVector=False):
        latNum=len(latCoords)
        lonNum=len(lonCoords)
        rNum=len(rCoords)
        valShp=np.shape(vals)
        if(rNum!=valShp[0]+1 or latNum!=valShp[1]+1 or lonNum!=valShp[2]+1):
            #if latDim * lonDim * rDim < 1 or valDim!=latDim*lonDim*rDim:
            print ("Voxel dim is not match with data dim")
            return None

        pts = vtk.vtkPoints()
        cells = vtk.vtkCellArray()
        scals = vtk.vtkFloatArray()
        numOfComp=valShp[3]
        scals.SetNumberOfComponents(numOfComp)
        scals.SetName(valName)
        #cellNum = cls.__GenUnStructHexahedronCells(lons, lats, rs, vals, cells, points, sals)
        for k in range(0, rNum):  # vtk将属性数据定义在顶点上，为此去掉了第一个顶点，以保持顶点的数量跟属性数据的数量一致
            for j in range(0, latNum):
                for i in range(0, lonNum):
                    pt = cls.__sphere_llr_geocent(lonCoords[i], latCoords[j], rCoords[k])
                    pts.InsertNextPoint(pt)
                    #if numOfComp > 1:
                    #    scals.InsertNextTypedTuple(vals[k, j, i])
                    #else:
                    #    scals.InsertNextTypedTuple(vals[k, j, i])
        # lonNum=lonNum-1#每个维度去掉了一个顶点，格网数量也应相应减少1个
        # latNum=latNum-1
        # rNum=rNum-1
        latLonNum = lonNum * latNum
        count = 0
        for iR in range(0, rNum-1):
            for iLat in range(0, latNum-1):
                for iLon in range(0, lonNum-1):
                    hex = vtk.vtkHexahedron()
                    hex.GetPointIds().SetId(0, iR * latLonNum + iLat * lonNum + iLon )
                    hex.GetPointIds().SetId(1, iR * latLonNum + iLat * lonNum + iLon + 1 )
                    hex.GetPointIds().SetId(2, iR * latLonNum + (iLat + 1) * lonNum + iLon + 1 )
                    hex.GetPointIds().SetId(3, iR * latLonNum + (iLat + 1) * lonNum + iLon )
                    hex.GetPointIds().SetId(4, (iR + 1) * latLonNum + iLat * lonNum + iLon )
                    hex.GetPointIds().SetId(5, (iR + 1) * latLonNum + iLat * lonNum + iLon + 1 )
                    hex.GetPointIds().SetId(6, (iR + 1) * latLonNum + (iLat + 1) * lonNum + iLon + 1 )
                    hex.GetPointIds().SetId(7, (iR + 1) * latLonNum + (iLat + 1) * lonNum + iLon )
                    cells.InsertNextCell(hex)
                    scals.InsertNextTypedTuple(vals[iR, iLat, iLon])


                    count = count + 1



        '''
        if quality==0:
            
        else:

            numOfSection = int(math.pow(2, quality))
            latCoords=[]
            lonCoords=[]
            for iLat in range(0,latDim):
                latRes = (lats[iLat+1] - lats[iLat]) / numOfSection
                for i in range(0,numOfSection):
                    lat=lats[iLat]+i*latRes
                    latCoords.append(lat)
            latCoords.append(lats[latDim])
            for iLon in range(0, lonDim):
                lonRes = (lons[iLon + 1] - lons[iLon]) / numOfSection
                for i in range(0, numOfSection):
                    lon = lons[iLon] + i * lonRes
                    lonCoords.append(lon)
            lonCoords.append(lons[lonDim])
            rCoords=rs
            count=0
            cellNum=0
            for iR in range(0,rDim):
                for iLat in range(0,latDim):
                    for iLon in range(0,lonDim):
                        cellNum=cellNum+cls.__GenUnStructHexahedronCells(lonCoords[iLon*numOfSection:(iLon+1)*numOfSection+1],
                                                         latCoords[iLat*numOfSection:(iLat+1)*numOfSection+1],
                                                         rCoords[iR*numOfSection:(iR+1)*numOfSection+1],
                                                         [vals[count]],cells,points,sals)
                        count=count+1
        '''


        if count>0:
            uGrid = vtk.vtkUnstructuredGrid()
            uGrid.SetPoints(pts)
            uGrid.SetCells(vtk.VTK_HEXAHEDRON, cells)

            if beVector:
                uGrid.GetCellData().SetVectors(scals)
            else:
                uGrid.GetCellData().SetScalars(scals)

            uGrid.Modified()
            if vtk.VTK_MAJOR_VERSION <= 5:
                uGrid.Update()
            return uGrid
        else:
            return None

    @classmethod
    def __GenUnStructHexahedronCells1(cls, lonCoords, latCoords, rCoords,vals,cells,pts,scals):
        '''
        :param lonCoords: list对象
        :param latCoords:list对象
        :param rCoords:
        :param vals: list对象，每个元素是一个truple
        :param cells:
        :param pts:
        :param scals:
        :return:
        '''
        pass
        '''
        lonNum=np.shape(lonCoords)[0]-1
        latNum=np.shape(latCoords)[0]-1
        rNum=np.shape(rCoords)[0]-1
        valNum = np.shape(vals)[3]
        ptstartid=pts.GetNumberOfPoints()
        count=0
        for k in range(0,rNum):#vtk将属性数据定义在顶点上，为此去掉了第一个顶点，以保持顶点的数量跟属性数据的数量一致
            for j in range(0,latNum):
                for i in range(0,lonNum):
                    pt=cls.__sphere_llr_geocent(rCoords[k],latCoords[j],lonCoords[i])
                    pts.InsertNextPoint(pt)
                    if valNum>1:
                        scals.InsertNextTypedTuple(vals[k,j,i])
                    else:
                        scals.InsertNextTypedTuple(vals[0])
                    count=count+1
        #lonNum=lonNum-1#每个维度去掉了一个顶点，格网数量也应相应减少1个
        #latNum=latNum-1
        #rNum=rNum-1
        latLonNum=lonNum*latNum
        count=0
        for iR in range(0,rNum-1):
            for iLat in range(0,latNum-1):
                for iLon in range(0,lonNum-1):
                    hex = vtk.vtkHexahedron()
                    hex.GetPointIds().SetId(0, iR*latLonNum+iLat*lonNum+iLon+ptstartid)
                    hex.GetPointIds().SetId(1, iR*latLonNum+iLat*lonNum+iLon+1+ptstartid)
                    hex.GetPointIds().SetId(2, iR*latLonNum+(iLat+1)*lonNum+iLon+1+ptstartid)
                    hex.GetPointIds().SetId(3, iR*latLonNum+(iLat+1)*lonNum+iLon+ptstartid)
                    hex.GetPointIds().SetId(4, (iR+1)*latLonNum+iLat*lonNum+iLon+ptstartid)
                    hex.GetPointIds().SetId(5, (iR+1)*latLonNum+iLat*lonNum+iLon+1+ptstartid)
                    hex.GetPointIds().SetId(6, (iR+1)*latLonNum+(iLat+1)*lonNum+iLon+1+ptstartid)
                    hex.GetPointIds().SetId(7, (iR+1)*latLonNum+(iLat+1)*lonNum+iLon+ptstartid)
                    cells.InsertNextCell(hex)

                    count=count+1
        return count
        '''

    @classmethod
    def __genQuadsPtsForSurface(cls, x0, x1, dx, y0, y1, dy, cz, xyzDimOrder, pts, faces, val):
        import copy
        nnx = (x1 - x0) / dx
        nny = (y1 - y0) / dy
        if nnx - int(nnx) > 0:
            nx = int(nnx) + 1
        else:
            nx = int(nnx)
        if nny - int(nny) > 0:
            ny = int(nny) + 1
        else:
            ny = int(nny)
        xcoords = []
        for i in range(0, nx):
            xcoords.append(x0 + dx * i)
        xcoords.append(x1)
        ycoords = []
        for i in range(0, ny):
            ycoords.append(y0 + dy * i)
        ycoords.append(y1)

        newpt = [0, 0, 0]
        newpt[xyzDimOrder[2]] = cz
        ptIndex = len(pts)
        for j in range(0, ny + 1):
            newpt[xyzDimOrder[1]] = ycoords[j]
            for i in range(0, nx + 1):
                newpt[xyzDimOrder[0]] = xcoords[i]
                pts.append(copy.deepcopy(newpt))

        for j in range(0, ny):
            for i in range(0, nx):
                faces.append(
                    [ptIndex + j * (nx + 1) + i, ptIndex + j * (nx + 1) + i + 1, ptIndex + (j + 1) * (nx + 1) + i + 1,
                     ptIndex + (j + 1) * (nx + 1) + i, val])

    @classmethod
    def __sphere_llr_geocent(cls, phi, lamda, r):
        import math
        x = r * math.cos(lamda) * math.cos(phi)
        y = r * math.cos(lamda) * math.sin(phi)
        z = r * math.sin(lamda)
        return [x, y, z]
    @classmethod
    def __gen_llr_voxel_surface(cls, minLon, maxLon, minLat, maxLat, minR, maxR, numOfSection, face_pts, faces,
                                usr_data):
        lonRes = (maxLon - minLon) / numOfSection
        latRes = (maxLat - minLat) / numOfSection
        cls.__genQuadsPtsForSurface(minLon, maxLon, lonRes, minLat, maxLat, latRes, minR, [0, 1, 2], face_pts, faces,
                                     usr_data)  # 底球面
        cls.__genQuadsPtsForSurface(minLon, maxLon, lonRes, minLat, maxLat, latRes, maxR, [0, 1, 2], face_pts, faces,
                                     usr_data)  # 顶球面
        cls.__genQuadsPtsForSurface(minLat, maxLat, latRes, minR, maxR, maxR - minR, minLon, [1, 2, 0], face_pts,
                                     faces, usr_data)  # 西侧面
        cls.__genQuadsPtsForSurface(minLat, maxLat, latRes, minR, maxR, maxR - minR, maxLon, [1, 2, 0], face_pts,
                                     faces, usr_data)  # 东侧面
        cls.__genQuadsPtsForSurface(minLon, maxLon, lonRes, minR, maxR, maxR - minR, minLat, [0, 2, 1], face_pts,
                                     faces, usr_data)  # 南侧面
        cls.__genQuadsPtsForSurface(minLon, maxLon, lonRes, minR, maxR, maxR - minR, maxLat, [0, 2, 1], face_pts,
                                     faces, usr_data)  # 北侧面



