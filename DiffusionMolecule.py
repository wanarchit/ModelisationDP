#!/usr/bin/env python

#
# Authors : Delphine Rousse, Paul Gand
# Date of last modification : 2017-10-16
#
# This program allows calculating the traveled distances by a selected
# reference atom for each selected range of time. The distances are
# clustered by a slice of distance chosen by the user.
#
# Read the README.txt file for more information.
#

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
    Allows the user to choose a slice value of distance to do the clusters.
    This function is contains secure user entries : if the user enters
    something other than an integer or a float then he must start over.

    :return selected_slice: The selected slice value transformed into float.
    """
    check = False
    selected_slice = ""
    while not check:
        print "Select a range of distances to do the clusters (min = 0.01, " \
              "max = 10), for the tests we used 0.1:"
        selected_slice = raw_input(" >>  ")
        try:  # The min and max are fixed arbitrarily to have enough clusters.
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
    This function allows building the list of all atoms.
    This list is necessary to the function 'select_ref_atom' to check if the
    user selects an existing reference atom.

    :param trajectory_data: Trajectory object from PDB and XTC files
    obtained by function load of MDtraj.
    :return atom_list: The list of all atoms (no duplicates)
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
    This function allows building the list of all residues which contain the
    selected reference atom.

    :param trajectory_data: Trajectory object from PDB and XTC files
    obtained by function load of MDtraj.
    :param ref_atom: The reference atom chosen by the user.
    :return residues_list: The list of all residues which contain the
    selected reference atom.
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
    Allows recovering the coordinates X and Y for a given time frame.

    :param trajectory_data: Trajectory object from PDB and XTC files
    obtained by function load of MDtraj.
    :param ref_atom: The reference atom chosen by the user.
    :param frame: The frame of time
    :return coordinate: A list of two values : X and Y coordinates.
    """
    coordinate = {}
    for atom in trajectory_data.topology.atoms:
        if atom.name == ref_atom:
            coordinate[atom.residue] = \
                list(trajectory_data.xyz[frame, atom.index])[0:2]
    return coordinate


def calculate_distances(time_list, trajectory_data, residue_list, ref_atom):
    """
    Allows calculating all distances between the position of reference atom of
    each residue and the position of this same atom at greater time interval.
    For example: the first given time is equal to 5 then the function
    calculate the distances between 0 and 5, 1 and 6, 2 and 7, ...

    :param time_list: The list of selected ranges of time
    :param trajectory_data: Trajectory object from PDB and XTC files
    obtained by function load of MDtraj.
    :param residue_list: The list of all residues which contain the
    selected reference atom.
    :param ref_atom: The reference atom chosen by the user.
    :return distance: Dictionary of each distance for each time for each
    residue.
    """
    print("Calculating distances in progress ...")
    distance = {}
    for time_i in time_list:
        print("For the range of time = " + str(time_i) + " ...")
        distance[time_i] = {}
        for residue in residue_list:
            distance[time_i][residue] = []
        t_start = 0
        t_end = time_i * 10  # 1 frame = 200 picosecond --> x10 to have 1
        # nanosecond
        while t_end <= 2000:
            coordinate_1 = get_coordinate(trajectory_data, ref_atom, t_start)
            # Coordinates X and Y for the start frame
            coordinate_2 = get_coordinate(trajectory_data, ref_atom, t_end)
            # Coordinates X and Y for the last frame of the time range
            for residue in residue_list:
                dist_residue_i = sqrt((coordinate_1[residue][0] -
                                       coordinate_2[residue][0])**2 + (
                    coordinate_1[residue][1] - coordinate_2[residue][0])**2)
                # Calculation of the distance : d = square_root((x1-x2)^2 +
                # (y1-y2)^2)
                distance[time_i][residue].append(dist_residue_i)
            t_start += 10
            t_end += 10
    print("Calculating distances: ok")
    return distance


def min_max_distances(dist, residues_list):
    """
    Allows recovering the minimal and maximal distance values among the
    list of all distances for each residue.

    :param dist: The list of each distance for each time for each residue.
    :param residues_list: The list of all residues which contain the
    selected reference atom.
    :return min_value: The minimal value of distance
    :return max_value: The maximal value of distance any residues included
    """
    print("Calculating min/max in progress ...")
    minimal_values = []
    maximal_values = []
    for time_i in dist.keys():  # For each range time we add the largest /
        # smallest distance value for each residue.
        for residue in residues_list:
            minimal_values.append(min(dist[time_i][residue]))
            maximal_values.append(max(dist[time_i][residue]))
    min_value = min(minimal_values)  # We keep only the minimal and maximal
    # distance values among all the time.
    max_value = max(maximal_values)
    print("Calculating min/max: ok")
    return min_value, max_value


def distance_clusters(slice_used, dist, residues_list):
    """
    Allows sorting all the distance values. The clusters are built according
    to the selected slice of distance value.

    For example: if the selected slice of distance is equal to 0.5 Angstrom
    then we will have a first cluster for all distances between 0 and 0.5 (
    non include), a second cluster for all distances between 0.5 and 1.0, ...

    For each range of distance we count the number of residues for each
    frame of time that traveled this distance.

    :param slice_used: The selected slice of distance value.
    :param dist: the list of all distances for each time and each residue.
    :param residues_list: The list of all residues which contain the
    selected reference atom.
    :return new_clusters: The table of counting distances traveled by all
    residues and for all ranges of time.
    """
    mini, maxi = min_max_distances(distances, residues_list)
    print("Building clusters in progress ...")
    new_clusters = {}
    for time_i in dist.keys():
        new_clusters[time_i] = {}
        cluster = 0
        while cluster <= round(maxi, 2):  # We create all cluster possibilities
            new_clusters[time_i][round(cluster, 2)] = 0
            cluster += slice_used
        for residue in residues_list:
            for value in distances[time_i][residue]:
                new_clusters[time_i][round(value - value % slice_used, 2)] += 1
                # For each distance, we calculating to which cluster it belong.
    print("Building clusters: ok")
    return new_clusters


#
# Output function
#


def export_clusters_csv(data_path, cluster_list):
    """
    Allows writing into a csv file the counting table.
    The output file name is generate according to the date and the hour of
    execution of the program. This allows it to be unique.

    :param data_path: The PATH of output directory
    :param cluster_list: The table of counting distances traveled by all
    residues and for all ranges of time.
    :return csv_file_name: The name of output CSV file.
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
    """ The core of program """
    # First we retrieve the arguments given as parameter.
    parser = argparse.ArgumentParser(description='This program allows '
                                                 'calculating and represent '
                                                 'the distribution of '
                                                 'distances traveled by '
                                                 'molecules through a system')
    parser.add_argument('pdb', type=str, help='PDB input file')
    parser.add_argument('xtc', type=str, help='XTC input file')
    parser.add_argument('rscript_path', type=str, help='Rscript PATH')
    args = parser.parse_args()

    # Next, we check if the arguments are in conformity
    pdb_file, xtc_file, rscript_path = check_args(args.pdb, args.xtc,
                                                  args.rscript_path)
    result_dir = os.path.commonprefix([xtc_file, pdb_file])
    main_dir = os.path.dirname(__file__)
    script_r = os.path.join(main_dir, 'ScriptR.R')

    # We load the data using the module MDtraj from the given PDB and XTC
    # files in parameters. Recovery of trajectory object.
    print("Load data in progress ...")
    trajectory = md.load(xtc_file, top=pdb_file)
    print("Load data: ok")

    # We build and retrieve all the data necessary for the calculations:
    # - The list of all available atoms
    list_of_atoms = build_atom_list(trajectory)
    # - The atom chosen by the user
    reference_atom = select_ref_atom(list_of_atoms)
    # - The list of all residues with this reference atom
    list_of_residues = build_residues_list(trajectory, reference_atom)
    # - The range of time between two frames for calculations.
    time = select_times()

    # We calculate each distance for each residue and for each time range
    distances = calculate_distances(time, trajectory, list_of_residues,
                                    reference_atom)

    # We calculate clusters of distances according to the slice of
    # distance chosen by the user.
    slice = select_slice()
    clusters = distance_clusters(slice, distances, list_of_residues)

    # We export the results (counting table) into csv file
    csv_name = export_clusters_csv(result_dir, clusters)

    # We execute the R script to analyze the results.
    call([rscript_path, script_r, csv_name, reference_atom])
    print("\nThe program has finished, the results are saved in the "
          "directory : "+result_dir)