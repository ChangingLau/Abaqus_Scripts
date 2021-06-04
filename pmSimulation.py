# -----------------------------------------------------------------------------------
# import modules
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import regionToolset
import visualization
import odbAccess
from datetime import *
import csv
import copy
import math


# define a class named as Simulation, with necessary input:
#   pmName = name of pavement;
#   dimension = [horizontal, vertical] dimension of pavement;
#   pmLayers = {layeri: [], ...} layers of pavement;
#   n_Layers = number of layers;
#   seedInfo = information of seed, [global seed size, max local seed size, min local seed size]

class Simulation:
    def __init__(self, pm, dim, seedInfo, ):
        self.pms = pm  # overall pavement information
        self.dims = dim  # overall pavement dimensions
        self.seedInfos = seedInfo  # overall model seed information
        self.modelSize = None  # size of model
        self.modelName = None  # name of model
        self.modelLayers = None  # layers dictionary of model
        self.nLayers = None  # number of layers
        self.seed = None  # seed information for seeding
        self.path = {}  # path information for XY data
        self.modelFace = {}  # model face dictionary of model
        self.modelEdge = {}  # model edge dictionary of model
        self.surfacePartitionEdge = {}  # rectangular partition edge points at model surface
        self.surfacePartitionVertex = {}  # rectangular partition vertices at model surface
        self.loadFaceEdge = None  # edge at load area
        self.loadFaceVertex = None  # vertex at load area
        self.loadFace = None  # Face of load area
        self.pmModel = None  # Model
        self.pmPart = None  # Part
        self.pmDatum = None  # Datums
        self.pmAssembly = None  # Assembly
        self.pmInstance = None  # Instance
        self.pmJob = None  # Job
        self.startTime = None  # start time of each model
        self.duration = 0  # calculation duration of each model

    def log(self, work, judge = False, mode = "a"):
        # writing log
        if not judge:
            endTime = datetime.now()
            duration = (endTime - self.startTime).total_seconds() 
            self.duration += duration
            log = [work, self.startTime.time(), endTime.time(), duration, ]
        else:
            log = [work, judge, "", "", ]
        with open("pmLog.csv", mode, ) as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            spamwriter.writerow(log)

    def initialVariables(self):
        # initial some varibales
        self.modelSize = None  # size of model
        self.modelName = None  # name of model
        self.modelLayers = None  # layers dictionary of model
        self.nLayers = None  # number of layers
        self.seed = None  # seed information for seeding
        self.path = {}  # path information for XY data
        self.modelFace = {}  # model face dictionary of model
        self.modelEdge = {}  # model edge dictionary of model
        self.surfacePartitionEdge = {}  # rectangular partition edge points at model surface
        self.surfacePartitionVertex = {}  # rectangular partition vertices at model surface
        self.loadFace = None  # vertex at load area

    # Diagram of pavement, and coding of faces and edges
    #                              / Y
    #                             /
    #                            /
    #             +-------------+                     ___________________________
    #            /|            /|                     ___________________________ HMA 1 (NO.: 0)
    #           / |           / |
    #          /  |   3      /  |                     ___________________________ HMA 2 (NO.: 1)
    #  X      /   |         /   |
    #  <-----+-------------+    |
    #        |    |        |    |                     ___________________________ Base (or Functional Layer)
    #        |    |        |  1 |                                                 (NO.: 2)
    #        |    =        |    =
    #        |             |
    #        |      2      |
    #        |             |                          ___________________________ Subbase (NO.: 3)
    #        |             |
    #        =             =                                                      Subgrade (NO.: 4)
    #
    #     No of Faces of a layer    / Y
    #                              /
    #                             /
    #              +-------------+
    #             /|            /|
    #            / |           / |
    #           /  |          1  |
    #  X       /   |         /   5
    #  <------+------2------+    |
    #         |    |        |    |
    #         |    |        |    |
    #         |    =        |    =
    #         3             4
    #         |             |
    #         |             |
    #         |             |
    #         =             =
    #
    #     No of Edges of a layer
    #                             ^  Y                                          * * + 1
    #                             |                                        * *      |
    #      +----------------------+                                    * *          |
    #      |   4                  |                                3                |
    #      |  +-------------------+                            *                    |
    #      |  |                   |                          *       ^              1
    #      |  |   3               |                         *        |              |
    #      |  |      +------------+      partition        *          arc            |
    #      |  |      |            |    +------------>  *                            |
    #      |  |      |   2  +-----+  _/               *                             |
    #      |  |      |      |  1  |                   +---------------2-------------+ 0
    #  <---+--+------+------+-----+                   2
    #  X
    #     No of circular partition at surface
    #

    # ====================================================================================  
    def generateGeoInfo(self):
        print("{}Here We go{}".format(2*"\n", "\n"))
        # -----------------------------------------------------------------------------------
        # determine thickness of subgrade layer
        # subtract depth of subgrade from total model hight
        self.modelLayers[self.nLayers - 1][2] = self.modelSize[1] - self.modelLayers[self.nLayers - 1][1]
        # -----------------------------------------------------------------------------------
        # generate model faces and edges dictionary
        hCoord = self.modelSize[0] / 2  # horizontal (x and y) coordinate
        for layer_i in self.modelLayers:
            # generate model faces points
            temp = {}
            vCoord = self.modelLayers[layer_i][1] + \
                self.modelLayers[layer_i][2] / 2  # vertical (z) coordinate
            temp[1] = [0, hCoord, vCoord]
            temp[2] = [hCoord, 0, vCoord]
            temp[3] = [hCoord, hCoord, self.modelLayers[layer_i][1]]
            self.modelFace[layer_i] = temp
            # modelFace = {
            #   layer_i: {
            #       1(X=0):[x,y,z],
            #       2(Y=0):[x,y,z],
            #       3(Z=0):[x,y,z],
            #   },
            #   ...
            # }
            # generate model edges points
            temp = {}
            temp[1] = [0, hCoord, self.modelLayers[layer_i][1]]
            temp[2] = [hCoord, 0, self.modelLayers[layer_i][1]]
            temp[3] = [hCoord * 2, 0, vCoord]
            temp[4] = [0, 0, vCoord]
            temp[5] = [0, hCoord * 2, vCoord]
            self.modelEdge[layer_i] = temp
            # modelEdge = {
            #   layer_i: {
            #       1: [x,y,z],
            #       2: [x,y,z],
            #       3: [x,y,z],
            #       4: [x,y,z],
            #       5: [x,y,z],
            #   }
            #   ...
            # }
        # print(self.modelFace)
        # print(self.modelEdge)
        # -----------------------------------------------------------------------------------
        # generate load area vertex, edge and face
        #
        #             Load area diagram
        #                                             =
        #                +------2-----+               |
        #               /              \              |
        #             /                  \            |
        #           /    half circle       \          |
        #         /                          \        |
        #        |                             \      |
        # -||----+--------------0---------------+--1--+
        #        2               0              1    
        # 
        self.loadFaceVertex = {
            0: [319.5 / 2, 0, 0],
            1: [53.25, 0, 0],
            2: [266.25, 0, 0], 
        }
        r = self.loadFaceVertex[2][0] - self.loadFaceVertex[0][0]  # initial circular partition radius
        self.loadFaceEdge = {
            0: self.loadFaceVertex[0],
            1: [self.loadFaceVertex[1][0] / 2, 0, 0],
            # 2: [self.loadFaceVertex[0][0], r, 0],
        }
        self.loadFace = [self.loadFaceVertex[0][0], r / 2, 0, ]
        # -----------------------------------------------------------------------------------
        # generate surface circular partition vertex dictionary
        cir_i = 1  # initial NO. of partition
        d = 106.5
        cirR = self.loadFaceVertex[2][0] + self.loadFaceVertex[1][0]
        # generate partition vertex
        #   while radius of circular partition is less than model horizontal size
        while cirR < self.modelSize[0] * 0.6:
            if cir_i != 1:
                cirR = 1.1 ** (cir_i - 1) * d + cirR  # radius determination 
            xCoord = yCoord = cirR
            # generate circular partition vertices
            temp = {}
            temp[1] = [0, yCoord, 0]
            temp[2] = [xCoord, 0, 0]
            self.surfacePartitionVertex[cir_i] = temp
            temp = {}
            temp[1] = [0, yCoord - 0.01, 0]
            temp[2] = [xCoord - 0.01, 0, 0]
            temp[3] = [xCoord*math.sqrt(2)/2, yCoord*math.sqrt(2)/2, 0]
            self.surfacePartitionEdge[cir_i] = temp
            cir_i += 1
        # redefine the 1st edge of 1st partition
        self.surfacePartitionEdge[1][2][0] = \
            (self.surfacePartitionEdge[1][2][0] + self.loadFaceVertex[2][0]) / 2
        # print(self.loadFace)
        # print(self.loadFaceVertex)
        # print(self.loadFaceEdge)
        # print(self.surfacePartitionVertex)
        # print(self.surfacePartitionEdge)
        # -----------------------------------------------------------------------------------
        # generate pathes
        depths = [0, 100, 300, 600, 1000, 5000, ]
        tempRange = [float(i) for i in range(0, int(self.modelSize[0]), 25)]
        xs = [0, ]
        xs[1: ] = tempRange[:]
        self.path = {
            # path on X axis at different depth
            # "depth": [],
            # "xCoordinate": [],
            "depth": depths,
            "xCoordinate": xs,
        }
    
    # ====================================================================================  
    def createModel(self):
        # -----------------------------------------------------------------------------------
        # create model
        self.pmModel = mdb.Model(self.modelName)
    
    # ====================================================================================  
    def createPart(self):
        # -----------------------------------------------------------------------------------
        # create part
        partSketch = self.pmModel.ConstrainedSketch(
            name='__profile__', sheetSize=8000)
        partSketch.rectangle(
            point1=(0.0, 0.0), point2=(self.modelSize[0], self.modelSize[0]))
        # create and assign the part to partPart
        self.pmPart = self.pmModel.Part(dimensionality=THREE_D, name=self.modelName,
                                        type=DEFORMABLE_BODY)
        self.pmPart.BaseSolidExtrude(
            depth=self.modelSize[1], sketch=partSketch)
        del partSketch
        session.viewports['Viewport: 1'].setValues(displayedObject=self.pmPart)
    
    # ====================================================================================  
    def createDatum(self):
        # -----------------------------------------------------------------------------------
        # create datum planes
        for layer_i in range(1, self.nLayers):
            plane_pt = self.modelFace[0][3]
            planeOrig = self.pmPart.faces.findAt((plane_pt,), )[0]  # [0]
            self.pmPart.DatumPlaneByOffset(flip=SIDE2,
                                           offset=self.modelLayers[layer_i][1], plane=planeOrig)
        self.pmDatum = self.pmPart.datums
    
    # ====================================================================================  
    def createPartition(self):
        # -----------------------------------------------------------------------------------
        # create layer partitions
        for layer_j in range(1, self.nLayers):
            cell = self.pmPart.cells.findAt(
                (self.modelFace[self.nLayers - 1][3],))[0]
            # create cell partition
            self.pmPart.PartitionCellByDatumPlane(cells=cell,
                                                  datumPlane=self.pmDatum[layer_j + 1])
        # -----------------------------------------------------------------------------------
        # create load and seed partition
        sketchFace = self.pmPart.faces.findAt((self.modelFace[0][3],))[0]
        sketch_UpEdge = self.pmPart.edges.findAt((self.modelEdge[0][1],))[0]
        partitionSketch = self.pmModel.ConstrainedSketch(gridSpacing=1, name='__profile__',
                                                        sheetSize=8000, \
                                                        transform=self.pmPart.MakeSketchTransform(
                                                             sketchPlane=sketchFace,
                                                             sketchPlaneSide=SIDE1,
                                                             sketchUpEdge=sketch_UpEdge,
                                                             sketchOrientation=RIGHT, 
                                                             origin=(0.0, 0.0, 0.0)))
        self.pmPart.projectReferencesOntoSketch(
            filter=COPLANAR_EDGES, sketch=partitionSketch)
        # -----------------------------------------------------------------------------------
        # sketch load partition
        pt0 = [a*-1 for a in self.loadFaceVertex[0][:2]]
        pt1 = [b*-1 for b in self.loadFaceVertex[1][:2]]
        pt2 = [c*-1 for c in self.loadFaceVertex[2][:2]]
        partitionSketch.ArcByCenterEnds(center=pt0, direction=COUNTERCLOCKWISE, \
            point1=pt1, point2=pt2)
        # -----------------------------------------------------------------------------------
        # sketch circular partition
        for cir_j in self.surfacePartitionVertex:
            pt1 = self.surfacePartitionVertex[cir_j][1][:2]
            pt2 = [j*-1 for j in self.surfacePartitionVertex[cir_j][2][:2]]
            partitionSketch.ArcByCenterEnds(center=(0, 0), direction=COUNTERCLOCKWISE,
                                            point1=pt1, point2=pt2)
        # create face partition
        self.pmPart.PartitionFaceBySketch(
            faces=sketchFace, sketch=partitionSketch, sketchUpEdge=sketch_UpEdge)
        del partitionSketch
    
    # ====================================================================================  
    def createMaterial(self):
        # -----------------------------------------------------------------------------------
        # create materials
        for layer_k in range(self.nLayers):
            materialName = self.modelLayers[layer_k][0]
            self.pmModel.Material(name=materialName)
            self.pmModel.materials[materialName].Elastic(table=((self.modelLayers[layer_k][3],
                                                                 self.modelLayers[layer_k][4]),))
    
    # ====================================================================================  
    def create_assignSection(self):
        # -----------------------------------------------------------------------------------
        # create and sections
        for layer_l in range(self.nLayers):
            sectionName = self.modelLayers[layer_l][0]
            self.pmModel.HomogeneousSolidSection(material=sectionName,
                                                 name=sectionName, thickness=None)
            # assign sections
            sectionRegion = (self.pmPart.cells.findAt(
                (self.modelFace[layer_l][1],))[0],)
            self.pmPart.SectionAssignment(offset=0.0, offsetField='',
                                          offsetType=MIDDLE_SURFACE, region=sectionRegion,
                                          sectionName=sectionName,
                                          thicknessAssignment=FROM_SECTION)
    
    # ====================================================================================  
    def createInstance(self):
        # -----------------------------------------------------------------------------------
        # create instance
        self.pmAssembly = self.pmModel.rootAssembly.DatumCsysByDefault(
            CARTESIAN)
        self.pmInstance = self.pmModel.rootAssembly.Instance(dependent=ON,
                                                             name=self.modelName, part=self.pmPart)
    
    # ====================================================================================  
    def createStep(self):
        # -----------------------------------------------------------------------------------
        # create step
        self.pmModel.StaticStep(name='ApplyingLoad', previous='Initial')
    
    # ====================================================================================  
    def createLoad(self):
        # -----------------------------------------------------------------------------------
        # create load
        loadRegion = ((self.pmInstance.faces.findAt((self.loadFace,)), SIDE1),)
        self.pmModel.Pressure(amplitude=UNSET, createStepName='ApplyingLoad',
                              distributionType=UNIFORM, field='', magnitude=0.7, name='VehicleLoad',
                              region=loadRegion)
    
    # ====================================================================================  
    def createBC(self):
        # -----------------------------------------------------------------------------------
        # create boundary condition
        # create bootom boundary condition
        bottomPt = self.modelFace[self.nLayers - 1][3]
        bottomPt[2] = self.modelSize[1]
        bottomBcRegion = (self.pmInstance.faces.findAt((bottomPt,))[0],)
        self.pmModel.DisplacementBC(amplitude=UNSET, createStepName='ApplyingLoad',
                                    distributionType=UNIFORM, fieldName='', fixed=OFF,
                                    localCsys=None, name='BottomBC', region=bottomBcRegion,
                                    u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        # define X=0 and Y=0 regions
        x0BcRegion = []  # X=0 region
        y0BcRegion = []  # Y=0 region
        # append X=0 and Y=0 regions using find
        for layer_m in self.modelFace:
            x0BcRegion.append(
                self.pmInstance.faces.findAt((self.modelFace[layer_m][1],)))
            y0BcRegion.append(
                self.pmInstance.faces.findAt((self.modelFace[layer_m][2],)))
        # -----------------------------------------------------------------------------------
        # create X=0 boundary condition
        self.pmModel.DisplacementBC(amplitude=UNSET, createStepName='ApplyingLoad',
                                    distributionType=UNIFORM, fieldName='', fixed=OFF,
                                    localCsys=None, name='X=0BC', region=x0BcRegion,
                                    u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        # -----------------------------------------------------------------------------------
        # create Y=0 boundary condition
        self.pmModel.DisplacementBC(amplitude=UNSET, createStepName='ApplyingLoad', distributionType=UNIFORM,
                                    fieldName='', fixed=OFF,
                                    localCsys=None, name='Y=0BC', region=y0BcRegion,
                                    u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    # ====================================================================================  
    def createMesh(self):
        # -----------------------------------------------------------------------------------
        # create mesh
        # define mesh control
        overallCellRegion = []
        for layer_n in self.modelFace:
            overallCellRegion.append(
                self.pmPart.cells.findAt((self.modelFace[layer_n][1],))[0])
        self.pmPart.setMeshControls(algorithm=MEDIAL_AXIS,
                                    regions=overallCellRegion, technique=SWEEP)
        self.pmPart.setElementType(elemTypes=(ElemType(elemCode=C3D20,
                                                       elemLibrary=STANDARD),
                                              ElemType(elemCode=C3D15,
                                                       elemLibrary=STANDARD),
                                              ElemType(elemCode=C3D10, elemLibrary=STANDARD)),
                                   regions=overallCellRegion)
        # -----------------------------------------------------------------------------------
        # seed the part
        # -----------------------------------------------------------------------------------
        # seed the part with global size
        self.pmPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1,
                             size=self.seed[0])
        # define an empty list to contain edges need to be assign constant seed size 1
        localConstantSizeEdge = []
        # add edges of layers except the bottom layer
        for layer_o in range(self.nLayers - 1):
            edgeList = self.modelEdge[layer_o]
            for edge_i in range(3, len(edgeList)+1):
                localConstantSizeEdge.append(
                    self.pmPart.edges.findAt((edgeList[edge_i],)[0])
                )
        # -----------------------------------------------------------------------------------
        # seed the part with constant size
        self.pmPart.seedEdgeBySize(
            constraint=FINER, deviationFactor=0.1,
            edges=localConstantSizeEdge,
            minSizeFactor=0.1, size=self.seed[2])
        # -----------------------------------------------------------------------------------
        # seed circular partition at surface
        # define an empty list to contain edges need to be assign constant seed size 1 to n
        seedMiniSize = self.seed[2]
        error = (self.seed[0] - seedMiniSize) / len(self.surfacePartitionEdge)
        for circle_k in self.surfacePartitionEdge:
            seedMiniSize = seedMiniSize + (circle_k - 1) * error
            edgeList = self.surfacePartitionEdge[circle_k]
            # define an empty list to contain edges need to be assign constant seed size 1 to n
            tempEdge = []
            # add load area when seeding the first circular partition
            if circle_k == 1:
                for edge_c in self.loadFaceEdge:
                    tempEdge.append(self.pmPart.edges.findAt((self.loadFaceEdge[edge_c],))[0])
                for edge_d in range(1, 3):
                    tempEdge.append(self.pmPart.edges.findAt((edgeList[edge_d],))[0])
            tempEdge.append(self.pmPart.edges.getClosest((edgeList[3],))[0][0])
            # for edge_j in edgeList:
            #     tempEdge.append(
            #         self.pmPart.edges.getClosest((edgeList[edge_j],))[0][0]
            #     )
            self.pmPart.seedEdgeBySize(
                constraint=FINER, deviationFactor=0.1,
                edges=tempEdge,
                minSizeFactor=0.1, size=seedMiniSize)
        # seed the part with gradational size
        # define an empty list to contain edges need to be assign gradational seed size
        localGradationalSizeEdge = []
        # add vertical edges of bottom layer
        edgeList = self.modelEdge[self.nLayers - 1]
        for edge_j in range(3, len(edgeList) + 1):
            localGradationalSizeEdge.append(
                self.pmPart.edges.findAt((edgeList[edge_j],))[0]
            )
        self.pmPart.seedEdgeByBias(
            biasMethod=SINGLE, constraint=FINER,
            end1Edges=localGradationalSizeEdge,
            maxSize=self.seed[1], minSize=self.seed[2])
        # generate mesh
        self.pmPart.generateMesh()
        self.pmModel.rootAssembly.regenerate()
        session.viewports['Viewport: 1'].setValues(displayedObject=self.pmPart)
    
    # ====================================================================================  
    def create_submitJob(self):
        # -----------------------------------------------------------------------------------
        # create job
        self.pmJob = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model=self.modelName, modelPrint=OFF,
                multiprocessingMode=THREADS, name=self.modelName,
                nodalOutputPrecision=SINGLE, numCpus=28, numDomains=28, numGPUs=2, \
                queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='',
                waitHours=0, waitMinutes=0)
        # -----------------------------------------------------------------------------------
        # submit job
        self.pmJob.submit(consistencyChecking=OFF)
        # Do not return control till job is finished running
        self.pmJob.waitForCompletion()
    
    # ====================================================================================  
    def accessDatabase(self):
        # -----------------------------------------------------------------------------------
        # open database file
        session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
            referenceRepresentation=ON)
        o1 = session.openOdb(
            name='{}.odb'.format(self.modelName), readOnly=False)
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)
        # -----------------------------------------------------------------------------------
        # create pathes
        for path_i in self.path["depth"]:
            path_name = "{}_Depth".format(path_i/1000.0)
            pathPointList = []
            for i in self.path["xCoordinate"]:
                temp = [0, 0, path_i]
                temp[0] = i
                pathPointList.append(temp)
            session.Path(name=path_name, type=POINT_LIST,
                         expression=(pathPointList))
        # -----------------------------------------------------------------------------------
        # set field output as S33
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S33'))
        # define S33 according to pathes defined before
        for path_j in self.path["depth"]:
            pth = session.paths["{}_Depth".format(path_j/1000.0)]
            S33_name = "{}{}_S33".format(self.modelName[0], path_j/1000.0)
            session.XYDataFromPath(name=S33_name, path=pth, includeIntersections=True,
                                   projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10,
                                   projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE,
                                   removeDuplicateXYPairs=True, includeAllElements=False)
            # output S33 report as "Point_S33"
            session.writeXYReport(fileName='{}_{}.rpt'.format(self.modelName, S33_name), appendMode=OFF,
                                  xyData=session.xyDataObjects[S33_name])
        # -----------------------------------------------------------------------------------
        # set field output as E33
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='E', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'E33'))
        # define E33 according to pathes defined before
        for path_k in self.path["depth"]:
            pth = session.paths["{}_Depth".format(path_k/1000.0)]
            E33_name = "{}{}_E33".format(self.modelName[0], path_k/1000.0)
            session.XYDataFromPath(name=E33_name, path=pth, includeIntersections=True,
                                   projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10,
                                   projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE,
                                   removeDuplicateXYPairs=True, includeAllElements=False)
            # output S33 report as "Point_E33"
            session.writeXYReport(fileName='{}_{}.rpt'.format(self.modelName, E33_name), appendMode=OFF,
                                  xyData=session.xyDataObjects[E33_name])
    
    # ====================================================================================  
    def main(self):
        self.startTime = datetime.now()
        self.log("Start...", self.startTime, "w", )
        for pm_i in self.pms:
            self.initialVariables()
            self.modelLayers = self.pms[pm_i]
            self.nLayers = len(self.modelLayers)
            for dim_i in self.dims:
                self.modelSize = [i * 106.5 for i in dim_i]
                for seed_i in self.seedInfos:
                    self.startTime = datetime.now()
                    self.seed = self.seedInfos[seed_i]
                    # model name = [Pavement Name][Horizontal dimension factor]_[Global seed size]
                    self.modelName = "{}{}_{}".format(pm_i, \
                        int(dim_i[0]), int(self.seed[0] / 100))
                    self.generateGeoInfo()
                    self.createModel()
                    self.createPart()
                    self.createDatum()
                    self.createPartition()
                    self.createMaterial()
                    self.create_assignSection()
                    self.createInstance()
                    self.createStep()
                    self.createLoad()
                    self.createBC()
                    self.createMesh()
                    self.create_submitJob()
                    self.accessDatabase()
                    self.log(self.modelName)
        del mdb.models['Model-1']
        self.log("Finished!", self.duration, )


# pavement information
pm = {
    "ReversedBasePM": {
        # No: [name, depth, thickness, E, PR] # PR is Poisson's Ratio
        0: ["AC16", 0, 50, 9500, 0.25],
        1: ["AC25", 50, 70, 8500, 0.25],
        2: ["UGM", 120, 200, 350, 0.35],
        3: ["CTM", 320, 360, 12000, 0.35],
        4: ["Soil", 680, None, 70, 0.4],
    },
    "SemiRigidBasePM": {
        # No: [name, depth, thickness, E, PR] # PR is Poisson's Ratio
        0: ["AC16", 0, 50, 9500, 0.25],
        1: ["AC25", 50, 70, 8500, 0.25],
        2: ["CTM", 120, 360, 12000, 0.35],
        3: ["UGM", 480, 20, 350, 0.35],
        4: ["Soil", 680, None, 70, 0.4],
    },
    "UGMBasePM": {
        # No: [name, depth, thickness, E, PR] # PR is Poisson's Ratio
        0: ["AC16", 0, 50, 9500, 0.25],
        1: ["AC25", 50, 70, 8500, 0.25],
        2: ["UGM", 120, 560, 350, 0.35],
        3: ["Soil", 680, None, 70, 0.4],
    },
}

# model size, [horizontal times, vertical times]
h_v = [
    [10.0, 100], [15.0, 120], [20.0, 140], 
    [25.0, 160], [30.0, 180], [35.0, 200],
    [40.0, 220], [45.0, 240],
]
# seed size information
seedInfo = {
    # [globalSize, maximumLocalSize, minmumLocalSize]
    0: [500, 1000, 50, ],
    1: [400, 800, 40, ],
    2: [300, 600, 30, ],
    3: [200, 400, 20, ],
}


# define a function to start simulation
demo = Simulation(pm, h_v, seedInfo)
demo.main()
