#!/usr/bin/env python

import mdtraj as md
import os
import sys
import argparse
from math import sqrt
from time import localtime
from subprocess import call

#
# Check function
#


def check_args(pdb, xtc, rpath):
    """
    Check the existence of given file/application in parameters.

    :param pdb: The first argument is the PDB file PATH
    :param xtc: The second argument is the XTC file PATH
    :param rpath: The third argument is the Rscript application PATH
    :return pdb: The variable PDB file PATH
    :return xtc: The variable XTC PATH
    :return rpath: The variable Rscript PATH
    """

    if not(os.path.exists(pdb)):
        print("The selected PDB file does not exist")
    if not(os.path.exists(xtc)):
        print("The selected XTC file does not exist")
    if not(os.path.exists(rpath)):
        print("Thr selected Rscript PATH does not exist (check the double "
              "backslash '\\\\')")
    if not(os.path.exists(pdb)) or not(os.path.exists(xtc)) \
            or not(os.path.exists(rpath)):
        sys.exit()  # If one of file does not exist then the program closes
    return pdb, xtc, rpath


#
# Input functions
#


def select_ref_atom(atom_list):
    """
    Allows the user to choose the reference atom for the calculations.

    :param atom_list: The list of all available atoms in the PDB file
    :return atom: The selected atom
    """
    atom = ""
    check = False
    while not check:
        print "Select a reference atom (for the test we used P1):"
        atom = raw_input(" >>  ")
        if atom in atom_list:
            check = True
        else:
            print("The selected atom is not available in atom list, you have "
                  "to select another atom")
    return atom


def select_times():
    """
    Allows the user to choose the ranges of time will be used for calculations.
    The user can be choose default values : 1, 5, 10, 20, 30, 50, 70, 100, 200.
    This function is contains secure user entries.

    :return time_list: The list of selected ranges of time
    """
    time_list = []
    print("\nYou have to select a list of time ranges will be used for "
          "calculations")
    print("By default, the time ranges selected are: 1, 5, 10, 20, 30, 50, "
          "70, 100 and 200")
    while not time_list:  # While the list is empty
        print("Do you want keep these values or change them to choose "
              "your own values ? (Keep/Change)")
        choice = raw_input(" >>  ")
        if choice.lower() == "k" or choice.lower() == "keep":
            time_list = [1, 5, 10, 20, 30, 50, 70, 100, 200]
        elif choice.lower() == "c" or choice.lower() == "change":
            check = False
            while not check:
                print("\nEnter all the values (integers) you want at once, "
                      "separated by a space (ex: '1 2 3 10 50')")
                values = raw_input(" >>  ")
                try:  # The entry list is sorted and the duplicates are
                    # deleted.
                    time_list = sorted(list(set(
                        [int(v) for v in values.split(' ')])))
                    check = True
                except:
                    print("\nYou have to select only integers")
        else:
            print("Invalid choice, you have to enter 'keep' or 'change'")
    return time_list


def select_slice():
    """
    Allows
    :return:
    """
    check = False
    selected_slice = ""
    while not check:
        print "Select a range of distances to do the clusters (min = 0.01, " \
              "max = 10), for the tests we used 0.1:"
        selected_slice = raw_input(" >>  ")
        try:
            if 0.01 <= float(selected_slice) <= 10:
                check = True
            else:
                error = "small" if float(selected_slice) < 0.01 else "large"
                print("The selected value is too " + error + ", choose "
                                                             "another value\n")
        except:
            print("You have to select a value, not a character\n")
    return float(selected_slice)


#
# Help functions
#


def build_atom_list(trajectory_data):
    """

    :param trajectory_data:
    :return atom_list:
    """
    print("Building atom list in progress ...")
    atom_list = []
    for atom in trajectory_data.topology.atoms:
        if atom.name not in atom_list:
            atom_list.append(str(atom.name))
    print("Building atom list: ok")
    return atom_list


def build_residues_list(trajectory_data, ref_atom):
    """

    :param trajectory_data:
    :param ref_atom:
    :return residues_list:
    """
    print("Building residues list in progress ...")
    residues_list = []
    for atom in trajectory_data.topology.atoms:
        if atom.name == ref_atom:
            residues_list.append(atom.residue)
    print("Building residues list: ok")
    return residues_list


#
# Calculation functions
#


def get_coordinate(trajectory_data, ref_atom, frame):
    """

    :param trajectory_data:
    :param ref_atom:
    :param frame:
    :return coordinate:
    """
    coordinate = {}
    for atom in trajectory_data.topology.atoms:
        if atom.name == ref_atom:
            coordinate[atom.residue] = \
                list(trajectory_data.xyz[frame, atom.index])[0:2]
    return coordinate


def calculate_distances(time_list, trajectory_data, residue_list, ref_atom):
    """

    :param time_list:
    :param trajectory_data:
    :param residue_list:
    :param ref_atom:
    :return distance:
    """
    print("Calculating distances in progress ...")
    distance = {}
    for time_i in time_list:
        print("For the range of time = " + str(time_i) + " ...")
        distance[time_i] = {}
        for residue in residue_list:
            distance[time_i][residue] = []
        t_start = 0
        t_end = time_i * 10
        while t_end <= 2000:
            coordinate_1 = get_coordinate(trajectory_data, ref_atom, t_start)
            coordinate_2 = get_coordinate(trajectory_data, ref_atom, t_end)
            for residue in residue_list:
                dist_residue_i = sqrt((coordinate_1[residue][0] -
                                       coordinate_2[residue][0])**2 + (
                    coordinate_1[residue][1] - coordinate_2[residue][0])**2)
                distance[time_i][residue].append(dist_residue_i)
            t_start += 10
            t_end += 10
    print("Calculating distances: ok")
    return distance


def min_max_distances(dist, residues_list):
    """

    :param dist:
    :param residues_list:
    :return min_value:
    :return max_value:
    """
    print("Calculating min/max in progress ...")
    minimal_values = []
    maximal_values = []
    for time_i in dist.keys():
        for residue in residues_list:
            minimal_values.append(min(dist[time_i][residue]))
            maximal_values.append(max(dist[time_i][residue]))
    min_value = min(minimal_values)
    max_value = max(maximal_values)
    print("Calculating min/max: ok")
    return min_value, max_value


def distance_clusters(slice_used, dist, residues_list):
    """

    :param slice_used:
    :param dist:
    :param residues_list:
    :return new_clusters:
    """
    mini, maxi = min_max_distances(distances, residues_list)
    print("Building clusters in progress ...")
    new_clusters = {}
    for time_i in dist.keys():
        new_clusters[time_i] = {}
        cluster = 0
        while cluster <= round(maxi, 2):
            new_clusters[time_i][round(cluster, 2)] = 0
            cluster += slice_used
        for residue in residues_list:
            for value in distances[time_i][residue]:
                new_clusters[time_i][round(value - value % slice_used, 2)] += 1
    print("Building clusters: ok")
    return new_clusters


#
#
#


def export_clusters_csv(data_path, cluster_list):
    """

    :param data_path:
    :param cluster_list:
    :return csv_file_name:
    """
    print("Export CSV file in progress ...")
    time_result = str(localtime().tm_year)\
                  + '-' + str(localtime().tm_mon) \
                  + '-' + str(localtime().tm_mday) \
                  + '_'+str(localtime().tm_hour) \
                  + 'h' + str(localtime().tm_min) \
                  + 'min' + str(localtime().tm_sec)
    csv_file_name = os.path.join(data_path, 'Results_' + time_result + '.csv')
    csv_file = open(csv_file_name, 'w')
    head_line = 'Frame_t;Cluster;Effective\n'
    csv_file.write(head_line)
    for time_i in sorted(cluster_list.keys()):
        for cluster in sorted(cluster_list[cluster_list.keys()[0]]):
            csv_file.write(str(time_i) + ';' + str(cluster) + ';' + str(
                cluster_list[time_i][cluster]) + '\n')
        csv_file.close()
    print("Export CSV file: ok")
    return csv_file_name


#
# Main program
#


if __name__ == "__main__":
    """
    
    """
    parser = argparse.ArgumentParser(description='This program allows '
                                                 'calculating and represent '
                                                 'the distribution of '
                                                 'distances traveled by '
                                                 'molecules through a system')
    parser.add_argument('pdb', type=str, help='PDB input file')
    parser.add_argument('xtc', type=str, help='XTC input file')
    parser.add_argument('rscript_path', type=str, help='Rscript PATH')
    args = parser.parse_args()

    pdb_file, xtc_file, rscript_path = check_args(args.pdb, args.xtc,
                                                  args.rscript_path)
    result_dir = os.path.commonprefix([xtc_file, pdb_file])
    main_dir = os.path.dirname(__file__)
    script_r = os.path.join(main_dir, 'ScriptR.R')
    print("Load data in progress ...")
    trajectory = md.load(xtc_file, top=pdb_file)
    print("Load data: ok")
    list_of_atoms = build_atom_list(trajectory)
    reference_atom = select_ref_atom(list_of_atoms)
    list_of_residues = build_residues_list(trajectory, reference_atom)
    time = select_times()
    distances = calculate_distances(time, trajectory, list_of_residues,
                                    reference_atom)
    slice = select_slice()
    clusters = distance_clusters(slice, distances, list_of_residues)
    csv_name = export_clusters_csv(result_dir, clusters)
    call([rscript_path, script_r, csv_name])
