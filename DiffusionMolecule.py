import mdtraj as md
from math import sqrt

def selectRefAtom(atomList):
    check = False
    while not check:
        print "Select a reference atom :"
        atom = raw_input(" >>  ")
        if atom in atomList:
            check = True
        else:
            print("The selected atom is not available in atom list, you have "
                  "to select another atom")
    return atom

def selectSlice():
    check = False
    while not check:
        print "Select a range of time to do clusters of distances (min = 0.01):"
        slice = raw_input(" >>  ")
        if slice >= 0.01:
            check = True
        else:
            print("The selected value is too small, choose another value")
    return slice

def buildAtomList(trajectory):
    print("Building atom list in progress ...")
    atomList = []
    for atom in trajectory.topology.atoms:
        if atom.name not in atomList:
            atomList.append(str(atom.name))
    print("Building atom list : ok")
    return atomList

def buildResiduesList(trajectory,refAtom):
    print("Building residues list in progress ...")
    residuesList = []
    for atom in trajectory.topology.atoms:
        if atom.name == refAtom:
            residuesList.append(atom.residue)
    print("Building residues list : ok")
    return residuesList


def getCoordinate(trajectory, refAtom, frame):
    coordinate = {}
    for atom in trajectory.topology.atoms:
        if (atom.name == refAtom):
            coordinate[atom.residue] = list(trajectory.xyz[frame, atom.index])[
                                       0:2]
    return coordinate


def calculateDistances(time, trajectory, residueList, refAtom):
    print("Distance calculs in progress ...")
    distance = {}
    for timeI in time:
        distance[timeI] = {}
        for residu in residueList:
            distance[timeI][residu] = []
        tStart = 0
        tEnd = timeI * 10
        while tEnd <= 2000 :
            coordinate1 = getCoordinate(trajectory, refAtom, tStart)
            coordinate2 = getCoordinate(trajectory, refAtom, tEnd)
            for residu in residueList:
                distResI = sqrt((coordinate1[residu][0] -
                                 coordinate2[residu][0])**2 +
                                (coordinate1[residu][1] -
                                 coordinate2[residu][0])**2)
                distance[timeI][residu].append(distResI)
            tStart += 10
            tEnd += 10
    print("Distance calculs : ok")
    return distance

def minMaxDistances(distances, residuesList):
    print("Min/Max calcul in progress ...")
    minimalValues = []
    maximalValues = []
    for timeI in distances.keys():
        for residu in residuesList:
            print timeI
            print residu
            print distances[timeI][residu]
            minimalValues.append(min(distances[timeI][residu]))
            maximalValues.append(max(distances[timeI][residu]))
    minValue = min(minimalValues)
    maxValue = max(maximalValues)
    print("Min/Max calcul : ok")
    return minValue,maxValue

def distanceClusters(slice,distances, residuesList):
    min,max = minMaxDistances(distances, residuesList)
    print("Building clusters in progress ...")
    clusters = {}
    for timeI in distances.keys():
        clusters[timeI] = {}
        # parcours de 0 to slice (ex: slice = 0.5 : [0,0.5[, [0.5,1[, [1,
        # 1.5[...
        cluster = 0
        while cluster <= round(max,2):
            clusters[timeI][cluster] = 0
            cluster += slice
        for residu in residuesList:
            for value in distances[timeI][residu]:
                clusters[timeI][round(value-value%slice,2)] += 1
    print("Building clusters : ok")
    return(clusters)

def exportClustersCSV(data_path,clusters):
    print("Export CSV file in progress ...")
    csvFile = open(data_path+'\Results.csv','w')
    headline = 'Frame_t;Cluster;Effectif\n'
    csvFile.write(headline)
    for timeI in sorted(clusters.keys()):
        for cluster in sorted(clusters[clusters.keys()[0]]):
            csvFile.write(str(timeI)+';'+str(cluster)+';'+
                          str(clusters[timeI][cluster])+'\n')
    csvFile.close()
    print("Export CSV file : ok")



# Path Paul :
data_path= 'D:\Cours\Master 2\Modelisation Bioinformatique\p3_p4_p5_p8'
# Path Delphine :
# data_path = 'C:\Users\Asus\Documents\Master 2 GPhy\Semestre 1\Modelisation'
print("Load data in progress ...")
trajectory = md.load(data_path+'\md_200ns_OK.xtc',
                     top=data_path+'\start.pdb')
print("Load data : ok")
atomList = buildAtomList(trajectory)
#refAtom = selectRefAtom(atomList)
refAtom = 'P1'
residueList = buildResiduesList(trajectory, refAtom)
time = [1, 5, 10, 20, 30, 50, 70, 100, 200]
distances = calculateDistances(time, trajectory, residueList, refAtom)
#slice = selectSlice()
slice = 0.1
clusters = distanceClusters(slice, distances, residueList)
exportClustersCSV(data_path,clusters)