import mdtraj as md
import os
import sys
import argparse
from math import sqrt
from time import localtime
from subprocess import call

def checkArgs(pdb,xtc,rpath):
    if not(os.path.exists(pdb)):
        print("The selected PDB file does not exist")
    if not(os.path.exists(xtc)):
        print("The selected XTC file does not exist")
    if not(os.path.exists(rpath)):
        print("Thr selected Rscript PATH does not exist (check the double "
              "backslash '\\\\')")
    if not(os.path.exists(pdb)) or not(os.path.exists(xtc)) or not(os.path.exists(rpath)):
        sys.exit()
    return(pdb,xtc,rpath)

def selectRefAtom(atomList):
    check = False
    while not check:
        print "Select a reference atom (for the test we used P1):"
        atom = raw_input(" >>  ")
        if atom in atomList:
            check = True
        else:
            print("The selected atom is not available in atom list, you have "
                  "to select another atom")
    return atom

def selectTimes():
    time =[]
    print("\nYou have to select a list of time ranges will be used for "
          "calculs")
    print("By default, the time ranges selected are : 1, 5, 10, 20, 30, 50, "
          "70, 100 and 200")
    while time == []:
        print("Do you want keep these values or change them to choose your own "
              "values ? (Keep/Change)")
        choice = raw_input(" >>  ")
        if choice.lower() == "k" or choice.lower() == "keep":
            time = [1, 5, 10, 20, 30, 50, 70, 100, 200]
        elif choice.lower() == "c" or choice.lower() == "change":
            check = False
            while not check:
                print("\nEnter all the values (integers) you want at once, "
                      "separeted by a space (ex: '1 2 3 10 50')")
                values = raw_input(" >>  ")
                try:
                    time = sorted(list(set([int(v) for v in values.split(' '
                                                                         '')])))
                    check = True
                except:
                   print("\nYou have to select only integers")
        else:
            print("Invalide choice, you have to enter 'keep' or 'change'")
    print time
    return time

def selectSlice():
    check = False
    while not check:
        print "Select a range of distances to do the clusters (min = 0.01, " \
              "max = 10), for the tests we used 0.1:"
        slice = raw_input(" >>  ")
        try:
            if 0.01 <= float(slice) <= 10:
                check = True
            else:
                error = "small" if float(slice) < 0.01 else "large"
                print("The selected value is too "+error+", choose another "
                                                         "value\n")
        except:
            print("You have to select a value, not a character\n")
    return float(slice)

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
        print("For the range of time = "+str(timeI)+" ...")
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
            clusters[timeI][round(cluster,2)] = 0
            cluster += slice
        for residu in residuesList:
            for value in distances[timeI][residu]:
                clusters[timeI][round(value-value%slice,2)] += 1
    print("Building clusters : ok")
    return(clusters)

def exportClustersCSV(data_path,clusters):
    print("Export CSV file in progress ...")
    timeResult = str(localtime().tm_year)+'-'+str(localtime().tm_mon)+'-'\
                 +str(localtime().tm_mday)+'_'+str(localtime().tm_hour)\
                 +'h'+str(localtime().tm_min)+'min'+str(localtime().tm_sec)
    csvName = os.path.join(data_path,'Results_'+timeResult+'.csv')
    csvFile = open(csvName,'w')
    headline = 'Frame_t;Cluster;Effective\n'
    csvFile.write(headline)
    for timeI in sorted(clusters.keys()):
        for cluster in sorted(clusters[clusters.keys()[0]]):
            csvFile.write(str(timeI)+';'+str(cluster)+';'+
                          str(clusters[timeI][cluster])+'\n')
    csvFile.close()
    print("Export CSV file : ok")
    return csvName



if __name__=="__main__":
    parser=argparse.ArgumentParser(description='This program allows '
                                               'calculating and represent the '
                                               'distribution of distances '
                                               'traveled by molecules '
                                               'through a system')
    parser.add_argument('pdb', type=str, help='PDB input file')
    parser.add_argument('xtc', type=str, help='XTC input file')
    parser.add_argument('RscriptPath', type=str, help='Rscrip PATH')
    args=parser.parse_args()

    pdbfile, xtcfile, rscriptpath = checkArgs(args.pdb, args.xtc,
                                           args.RscriptPath)
    resultdir = os.path.commonprefix([xtcfile, pdbfile])
    maindir = os.path.dirname(__file__)
    scriptR = os.path.join(maindir,'ScriptR.R')
    print("Load data in progress ...")
    trajectory = md.load(xtcfile,top=pdbfile)
    print("Load data : ok")
    atomList = buildAtomList(trajectory)
    refAtom = selectRefAtom(atomList)
    #refAtom = 'P1'
    residueList = buildResiduesList(trajectory, refAtom)
    #time = [1, 5, 10, 20, 30, 50, 70, 100, 200]
    time = selectTimes()
    distances = calculateDistances(time, trajectory, residueList, refAtom)
    slice = selectSlice()
    #slice = 0.1
    clusters = distanceClusters(slice, distances, residueList)
    csvName = exportClustersCSV(resultdir,clusters)
    call([rscriptpath,scriptR,csvName])