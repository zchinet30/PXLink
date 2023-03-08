import subprocess as sp
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import logging
import re
import os
import numpy as np
# import math
# import random
# from sympy import Plane, Point3D


class GromacsSys:
    """This class is used to describe a Gromacs system.
    It includes all topology data read from .top file and
    coordinate data from .gro file, as well as other required
    data.

    Class includes these variables:
    topfile {string} -- [path of .top file]
    grofile {string} -- [path of .gro file]
    logfile {string} -- [path of data file]
    ndxfile {string} -- [path of index file for [frozen] group]
    frozen_resname {list} -- [list of names to be frozen]
    atoms {pandas dataframe} -- [atom data .top [atom]]
    coordinate {pandas dataframe} -- [coordinate data from .gro]
    time {float} -- [time, read from .gro file]
    velocity {bool} -- [whether .gro file contains velocity]
    box {list of float[3]} -- [box size vector from .gro]
    atomnum {int} -- [atom numbers in the system]
    bonds {list} -- [[bonds] from .top file]
    pairs {list} -- [[pairs] from .top file]
    angles {list} -- [[angles] from .top file]
    dihedrals {list} -- [(proper) [dihedrals] from .top file]
    impropers {list} -- [(improper) [dihedrals] from .top file
    NOTE: code will try to read and write parameters directly in .top file]
    carboxyl {list} -- [list of carboxyl groups, used in polymerization]
    amine {list} -- [list of amine groups, used in polymerization]
    charge_shift {list} -- [list of atoms whose charges are changed during
    polymerization. In case of PA crosslinking this is a list of
    [carboxyl C/ amine N, amide O/H, alpha aromatic carbon, beta atomatic carbons]]
    cl_num {int} -- [number of crosslinking reactions occured]
    rings {pandas dataframe} -- [benzene ring planes]
    """

    def __init__(self,
                 topfile,
                 grofile,
                 logfile,
                 tprfile='',
                 ndxfile='',
                 frozen_resname=[],
                 residue_connectivity=True,
                 xlink_list_file='residue_connectivity.txt',
                 old_gmx=True):
        """Initiates the class. Read data from .top and .gro file.
        Create and write in log file.

        Arguments:
            topfile {string} -- [path of .top file]
            grofile {string} -- [path of .gro file]
            logfile {string} -- [path of data file]
            tprfile {string} -- [path of .tpr file (optional).
            Will be renewed when calling gmx_run]
            ndxfile {string} -- [path of.ndx file (optional)]
            frozen_resname {list} -- [list of names to be frozen]
                The two arguments above are used when e.g. freezing atoms. Use one of them.
            residue_connectivity {bool} -- [If true a network will be created to record connectivity
            between residues] Default to true.
            xlink_list_file {string} -- [path of a file in which bonds between RESIDUES are recorded.
            Only used to record residue_connectivity]
            old_gmx {bool} -- [True if you are using Gromacs 4.x; False if higher versions.]
        """
        self.topfile = topfile
        self.grofile = grofile
        self.tprfile = tprfile
        self.time = 0
        self.old_gmx = old_gmx
        # Creating logging settings:
        self.logfile = logfile
        logging.basicConfig(filename=self.logfile,
                            filemode='a',
                            level="INFO",
                            format='%(asctime)s: %(levelname)s: %(message)s')
        logging.info(
            f"Created log file {self.logfile}, initiating GromacsSys class.")
        # Define the charge changes during amide bond formation.
        # delta_charge_C: charge change on carboxyl side, in atoms as follows in order:
        # amide C, amide O, alpha aromatic C, 2 beta aromatic C, 2 gamma aromatic C
        # delta_charge_N: charge change on amine side, in atoms as follows in order:
        # amide N, amide H, alpha aromatic C, 2 beta aromatic C, 2 gamma aromatic C
        # In further work, we need to make this easy to adjust by the user, so that
        # this program can be used for more types of monomers.
        self.delta_charge_C = [
            0.143, -0.043, -0.204, 0.018, 0.018, 0.001, 0.001
        ]
        self.delta_charge_N = [
            0.329, -0.043, 0.089, -0.006, -0.006, -0.002, -0.002
        ]

        self.read_top()
        self.read_gro()
        self.carboxyl = []
        self.amine = []
        self.cl_num = 0
        self.rings = pd.DataFrame(
            columns=['resnr', 'resname', 'ring_plane', 'center_point'])
        logging.info("Read topology data from " + topfile +
                     " and coordinate data from " + grofile + ".\n")
        # Read ndx file. The file should include a [ frozen ] group for freezing atoms.
        self.freeze_group = []
        if ndxfile:
            self.ndxfile = ndxfile
            with open(ndxfile) as f:
                for line in f:
                    if '[ frozen ]' in line:
                        break
                for line in f:
                    self.freeze_group = self.freeze_group + [
                        int(num) for num in line.split()
                    ]
                    if '[' in line:
                        break
            logging.info(
                f"Read frozen atom numbers from {ndxfile}, a total of {len(self.freeze_group)} atoms."
            )
            self.freeze_res = set(self.atoms.loc[self.atoms['nr'].isin(
                self.freeze_group)]['resnr'].values)
        else:
            self.ndxfile = ''
        if frozen_resname:
            for rn in frozen_resname:
                self.freeze_group = self.freeze_group + self.atoms.loc[
                    self.atoms['resname'] == rn].values
        # Might need more functions on this one.
        # Create network of residue connectivity
        self.residue_connectivity = residue_connectivity
        if residue_connectivity:
            self.Residuenet = nx.Graph()
            self.Residuenet.add_nodes_from(set(self.atoms['resnr'].values))
            logging.info(
                f"Created network of {len(set(self.atoms['resnr'].values))} residues."
            )
            if self.freeze_group:
                self.Residuenet.remove_nodes_from(self.freeze_res)
                logging.info(
                    f"Removed {len(self.freeze_res)} frozen residues.")
            # Keep residue types: 'T' for 'TMA' and 'M' for 'MPD'
            dic_restype = {
                i: self.atoms.loc[self.atoms['resnr'] == i]
                ['resname'].values[0][0]
                for i in set(self.atoms['resnr'].values)
            }
            nx.set_node_attributes(self.Residuenet, dic_restype,
                                   'residue_type')
            # Set the file to keep crosslink bonds between residues as edges
            self.xlink_list_file = xlink_list_file
            exist_edges = []
            # Read exist edges
            if os.path.isfile(self.xlink_list_file):
                with open(self.xlink_list_file, mode='r') as f:
                    for line in f:
                        if not line.strip():
                            continue
                        line_contents = line.strip().split()
                        exist_edges += [[
                            int(line_contents[0]),
                            int(line_contents[2])
                        ]]
            self.Residuenet.add_edges_from(exist_edges)

    def read_top(self):
        """Read topology data from file self.topology
        """
        logging.info(f"Running read_top to read topology file {self.topfile}")
        with open(self.topfile) as topf:
            # These lists stores information read from topology file.
            atoms = []
            bonds = []
            pairs = []
            angles = []
            dihedrals = []
            impropers = []
            cmap = []

            # These flags marks which part of the topology file is being read.
            parse_atoms = False
            parse_bonds = False
            parse_pairs = False
            parse_angles = False
            parse_dihedrals = False
            parse_impropers = False
            parse_cmap = False

            # Read topology line by line:
            for line in topf:
                # Change parse flag if reading specific lines:
                if line.strip() == "[ atoms ]":
                    parse_atoms = True
                    continue
                elif line.strip() == "[ bonds ]":
                    parse_atoms = False
                    parse_bonds = True
                    continue
                elif line.strip() == "[ pairs ]":
                    parse_bonds = False
                    parse_pairs = True
                    continue
                elif line.strip() == "[ angles ]":
                    parse_pairs = False
                    parse_angles = True
                    continue
                elif line.strip() == "[ dihedrals ]" and not parse_dihedrals:
                    parse_angles = False
                    parse_dihedrals = True
                    continue
                elif line.strip() == "[ dihedrals ]" and parse_dihedrals:
                    parse_dihedrals = False
                    parse_impropers = True
                    continue
                elif line.strip() == "[ cmap ]":
                    parse_dihedrals = False
                    parse_impropers = False
                    parse_cmap = True
                    continue
                elif line.startswith("#") and (parse_impropers or
                                               parse_dihedrals or parse_cmap):
                    parse_dihedrals = False
                    parse_impropers = False
                    parse_cmap = False
                    break
                # Read data:
                if not line.strip() or re.match(r'^\s*;', line):
                    # Skip empty/annotation lines.
                    continue
                if parse_atoms:
                    atoms.append(line.split()[0:8])
                    continue
                if parse_bonds:
                    bonds.append(line.split())
                    continue
                if parse_pairs:
                    pairs.append(line.split())
                    continue
                if parse_angles:
                    angles.append(line.split())
                    continue
                if parse_dihedrals:
                    dihedrals.append(line.split())
                    continue
                if parse_impropers:
                    impropers.append(line.split())
                    continue
                if parse_cmap:
                    cmap.append(line.split())
                    continue
            # Sorting out data, turn atoms list into pandas Dataframe and write column names.
            self.atoms = pd.DataFrame(atoms)
            self.atoms.columns = [
                'nr', 'type', 'resnr', 'resname', 'atom', 'cgnr', 'charge',
                'mass'
            ]
            col_int = ['nr', 'resnr', 'cgnr']
            col_float = ['charge', 'mass']
            self.atoms[col_int] = self.atoms[col_int].apply(
                pd.to_numeric, downcast='unsigned')
            self.atoms[col_float] = self.atoms[col_float].apply(
                pd.to_numeric, downcast='float')
            self.atomnum = self.atoms.shape[0]
            self.resnum = max(self.atoms['resnr'])
            # Write original atom numbers  and crosslinking marks:
            self.atoms['sep'] = ';'
            self.atoms['original_nr'] = self.atoms['nr']
            self.atoms['marks'] = ''
            # Convert strings in topology entries (bonds, angles, etc.) to numbers:
            self.bonds = [[int(i[0]), int(i[1]), int(i[2])] for i in bonds]
            self.pairs = [[int(i[0]), int(i[1]), int(i[2])] for i in pairs]
            self.angles = [[int(i[0]),
                            int(i[1]),
                            int(i[2]),
                            int(i[3])] for i in angles]
            self.dihedrals = [[
                int(i[0]),
                int(i[1]),
                int(i[2]),
                int(i[3]),
                int(i[4])
            ] for i in dihedrals]
            if impropers:
                if len(impropers[0]) == 7:
                    self.impropers = [[
                        int(i[0]),
                        int(i[1]),
                        int(i[2]),
                        int(i[3]),
                        int(i[4]),
                        int(i[5]),
                        float(i[6])
                    ] for i in impropers]
                else:
                    self.impropers = [[
                        int(i[0]),
                        int(i[1]),
                        int(i[2]),
                        int(i[3]),
                        int(i[4]), 0, 443.504
                    ] for i in impropers]
            else:
                self.impropers = []
            if cmap:
                self.cmap = [[
                    int(i[0]),
                    int(i[1]),
                    int(i[2]),
                    int(i[3]),
                    int(i[4]),
                    int(i[5])
                ] for i in cmap]
            else:
                self.cmap = []
            # Make a graph for molecular topology, where nodes are atoms and edges are bonds:
            self.Graph = nx.Graph()
            self.Graph.add_nodes_from(range(1, self.atomnum + 1))
            self.Graph.add_edges_from([row[0:2] for row in self.bonds])
            self.charge_shift = []
            logging.info("Loaded topology data from " + self.topfile + ".")

    def read_gro(self):
        """Reads the gro file designated by self.gro.
        """
        # Read .gro file
        logging.info(f"Running read_gro, reading from {self.grofile}")
        with open(self.grofile) as grof:
            # Read time data
            lineone = grof.readline().split()
            if len(lineone) >= 3:
                if lineone[-2] == 't=':
                    self.time = lineone[-1]
            # Read coordinate data
            grof.readline()
            linethree = grof.readline().split()
            if len(linethree) == 6:
                self.velocity = False
                cols = ['resnr', 'resname', 'atom', 'nr', 'x', 'y', 'z']
            elif len(linethree) == 9:
                self.velocity = True
                cols = [
                    'resnr', 'resname', 'atom', 'nr', 'x', 'y', 'z', 'vx',
                    'vy', 'vz'
                ]
            else:
                logging.error(
                    "Something is wrong with the .gro input file format.")
                raise ValueError(
                    "Something is wrong with the .gro input file format.")
        if self.velocity:
            self.coordinate = pd.read_fwf(self.grofile,
                                          colspecs=[(0, 5), (5, 9), (9, 15),
                                                    (15, 20), (20, 28), (28,
                                                                         36),
                                                    (36, 44), (44, 52),
                                                    (52, 60), (60, 68)],
                                          dtype={
                                              3: 'str',
                                              4: 'str',
                                              5: 'str',
                                              6: 'str',
                                              7: 'str',
                                              8: 'str'
                                          },
                                          header=None,
                                          skiprows=2,
                                          skipfooter=1,
                                          engine='python')
            self.coordinate[[3, 4, 5, 6, 7,
                             8]] = self.coordinate[[3, 4, 5, 6, 7,
                                                    8]].astype('float')
        else:
            self.coordinate = pd.read_fwf(self.grofile,
                                          colspecs=[(0, 5), (5, 9), (9, 15),
                                                    (15, 20), (20, 28),
                                                    (28, 36), (36, 44)],
                                          dtype={
                                              3: 'str',
                                              4: 'str',
                                              5: 'str'
                                          },
                                          header=None,
                                          skiprows=2,
                                          skipfooter=1,
                                          engine='python')
            self.coordinate[[3, 4, 5]] = self.coordinate[[3, 4,
                                                          5]].astype('float')
        self.coordinate.columns = cols
        self.coordinate['nr'] = self.coordinate['nr'].astype(int)
        self.coordinate[['resnr', 'nr']].apply(pd.to_numeric,
                                               downcast='unsigned')
        self.coordinate[['x', 'y', 'z']].apply(pd.to_numeric, downcast='float')
        if self.velocity:
            self.coordinate[['vx', 'vy', 'vz']].apply(pd.to_numeric,
                                                      downcast='float')
        # Read box size from the last line
        with open(self.grofile, 'rb') as grof:
            grof.seek(-2, os.SEEK_END)
            while grof.read(1) != b'\n':
                grof.seek(-2, os.SEEK_CUR)
            last_line = grof.readline().decode()
        self.box = [float(i) for i in last_line.split()]
        logging.info("Loaded coordinate data from" + self.grofile + ".")

    def mark_atoms(self, nofreeze=False):
        """
        Mark all free carboxyl and amine groups with the following procedures:
        1. Find all hydroxyl oxygen (OG311) and amine nitrogen (NG2S3)
        2. Find carboxyl and amine groups by looking for atoms nearby these atoms.
        3. Write list of carboxyl and amine groups in the given order:
        Carboxyl list self.carboxyl: each element of this list is a sub-list consisting of
        indices of atoms in one carboxyl group, in the following order: [carboxyl C, hydroxyl O, hydroxyl H]
        Amine list self.amine: each element of this list is a sub-list consisting of indices of
        atoms in one amine group, in the following order:
        [amine N, amine H (deleted when forming bond), amine H (remains when forming bond)]
        These two lists are used to modify topology when adding bond. Each bond is formed between atoms indicated
        by the first atoms in a self.carboxyl sub-list and a self.amine sub-list.
        The method also writes another list self.charge_shift.
        For carboxyl carbon, a sub-list consistes of indices in this order:
        [carboxyl C, carbonyl O, alpha aromatic C, beta aromatic C]
        For amine carbon, a sub-list consistes of indices in this order:
        [amine N, amine H (remains when forming bond), alpha aromatic C, beta aromatic C]
        This list marks the atoms whose charge are changed when forming bond, and is used when modifying charges.
        NOTE: This method needs to be modified if force field other then CHARMM is used.
        """
        logging.info("Running mark_atoms\n")
        # Find hydroxyl oxygen OG311 and amine nitrogen NG2S3
        O_atoms = self.atoms.loc[self.atoms['type'] == 'OG311']['nr'].tolist()
        N_atoms = self.atoms.loc[self.atoms['type'] == 'NG2S3']['nr'].tolist()
        if nofreeze and self.freeze_group:
            O_atoms = [i for i in O_atoms if i not in self.freeze_group]
            N_atoms = [i for i in N_atoms if i not in self.freeze_group]
        # Mark the carboxyl groups
        for i in O_atoms:
            hydroxyl_O_neighbor = list(nx.neighbors(self.Graph, i))
            # Find the two atoms connected with OG311 (which are carboxyl C and H),
            # and record them in order.
            if self.atoms.loc[self.atoms['nr'] == hydroxyl_O_neighbor[0]][
                    'type'].tolist()[0] == 'CG2O2':
                carboxyl_C = hydroxyl_O_neighbor[0]
                self.carboxyl.append([hydroxyl_O_neighbor[0]] + [i] +
                                     [hydroxyl_O_neighbor[1]])
            elif self.atoms.loc[self.atoms['nr'] == hydroxyl_O_neighbor[0]][
                    'type'].tolist()[0] == 'HGP1':
                carboxyl_C = hydroxyl_O_neighbor[1]
                self.carboxyl.append([hydroxyl_O_neighbor[1]] + [i] +
                                     [hydroxyl_O_neighbor[0]])
            # Find atoms around the carboxyl group and write the charge_shift list
            hydroxyl_C_neighbor = list(nx.neighbors(self.Graph, carboxyl_C))
            alpha_C = self.atoms.loc[
                (self.atoms['nr'].isin(hydroxyl_C_neighbor))
                & (self.atoms['type'] == 'CG2R61')]['nr'].tolist()[0]
            charge_list = [
                carboxyl_C,
                self.atoms.loc[(self.atoms['nr'].isin(hydroxyl_C_neighbor)) & (
                    self.atoms['type'] == 'OG2D1')]['nr'].tolist()[0], alpha_C
            ]
            beta_C_list = self.atoms.loc[
                (self.atoms['nr'].isin(list(nx.neighbors(self.Graph, alpha_C)))
                 ) & (self.atoms['type'] == 'CG2R61')]['nr'].tolist()
            beta_C_neighbors = []
            for j in beta_C_list:
                beta_C_neighbors = beta_C_neighbors + self.atoms.loc[
                    (self.atoms['nr'].isin(list(nx.neighbors(self.Graph, j))))
                    & (self.atoms['type'] == 'CG2R61')]['nr'].tolist()
            gamma_C_list = [
                k for k in beta_C_neighbors if k not in charge_list
            ]
            charge_list = charge_list + beta_C_list + gamma_C_list
            self.charge_shift.append(charge_list)

        # Mark the amine groups
        for i in N_atoms:
            amine_N_neighbor = list(nx.neighbors(self.Graph, i))
            amine_H_list = self.atoms.loc[
                (self.atoms['nr'].isin(amine_N_neighbor))
                & (self.atoms['type'] == 'HGP4')]['nr'].tolist()
            self.amine.append([i] + amine_H_list)
            # Find atoms around the amine group and write the charge_shift list
            alpha_C = self.atoms.loc[
                (self.atoms['nr'].isin(amine_N_neighbor))
                & (self.atoms['type'] == 'CG2R61')]['nr'].tolist()[0]
            beta_C_list = self.atoms.loc[
                (self.atoms['nr'].isin(list(nx.neighbors(self.Graph, alpha_C)))
                 ) & (self.atoms['type'] == 'CG2R61')]['nr'].tolist()
            charge_list = [i, amine_H_list[1], alpha_C]
            beta_C_neighbors = []
            for j in beta_C_list:
                beta_C_neighbors = beta_C_neighbors + self.atoms.loc[
                    (self.atoms['nr'].isin(list(nx.neighbors(self.Graph, j))))
                    & (self.atoms['type'] == 'CG2R61')]['nr'].tolist()
            gamma_C_list = [
                k for k in beta_C_neighbors if k not in charge_list
            ]
            charge_list = charge_list + beta_C_list + gamma_C_list
            self.charge_shift.append(charge_list)
        # Write log
        logging.info(f"Found {len(self.carboxyl)} carboxyl carbon atoms.\n")
        logging.info(f"Found {len(self.amine)} amine nitrogen atoms.\n")
        # logging.info(f"Carboxyl carbon: {' '.join(str(i[0]) for i in self.carboxyl)}")
        # logging.info(f"Amine nitrogen: {' '.join(str(i[0]) for i in self.amine)}\n")

    def remove_atoms(self, remove_list):
        """
        Remove all atoms whose indices are given in argument remove_list.
        Also removes all bond, pair, angle, dihedral and improper dihedral entries that
        contains theses atoms.
        Arguments:
            remove_list {list} -- [list of atoms to be removed]
        """
        logging.info(
            f"Running remove_atoms({' '.join(str(x) for x in remove_list)})")
        self.atoms = self.atoms[~self.atoms['nr'].isin(remove_list)]
        self.coordinate = self.coordinate[~self.coordinate['nr'].
                                          isin(remove_list)]
        bonds_r = [
            i for i in self.bonds if any([j in remove_list for j in i[0:2]])
        ]
        pairs_r = [
            i for i in self.pairs if any([j in remove_list for j in i[0:2]])
        ]
        angles_r = [
            i for i in self.angles if any([j in remove_list for j in i[0:3]])
        ]
        dihedrals_r = [
            i for i in self.dihedrals
            if any([j in remove_list for j in i[0:4]])
        ]
        impropers_r = [
            i for i in self.impropers
            if any([j in remove_list for j in i[0:4]])
        ]
        self.bonds = [i for i in self.bonds if i not in bonds_r]
        self.pairs = [i for i in self.pairs if i not in pairs_r]
        self.angles = [i for i in self.angles if i not in angles_r]
        self.dihedrals = [i for i in self.dihedrals if i not in dihedrals_r]
        self.impropers = [i for i in self.impropers if i not in impropers_r]
        self.amine = [
            i for i in self.amine if not (any([j in remove_list for j in i]))
        ]
        self.carboxyl = [
            i for i in self.carboxyl
            if not (any([j in remove_list for j in i]))
        ]
        self.Graph.remove_nodes_from(remove_list)
        self.atomnum = self.atomnum - len(remove_list)
        self.reset_index()
        # Write log:
        # logging.info(f"Removed atoms: {' '.join(str(x) for x in remove_list)}")
        # logging.info(f"Removed bonds: {' '.join(str(x) for x in bonds_r)}")
        # logging.info(f"Removed pairs: {' '.join(str(x) for x in pairs_r)}")
        # logging.info(f"Removed angles: {' '.join(str(x) for x in angles_r)}")
        # logging.info(f"Removed dihedrals: {' '.join(str(x) for x in dihedrals_r)}")
        # logging.info(f"Removed improper dihedrals: {' '.join(str(x) for x in impropers_r)}\n")

    def reset_index(self, startn=1):
        """
        Sets atom indices to consecutive numbers, beginning with startn,
        ending with startn + atomnum - 1.
        Also changes atom numbers in bond, pair, angle etc. entries.

        Arguments:
            startn {integer} -- [starting number of atom indices.] Default 1
        """
        logging.info("Running reset_index to reset atom indices.")
        index = dict(
            zip(self.atoms['nr'], range(startn, startn + self.atomnum)))
        self.atoms = self.atoms.replace({'nr': index})
        self.coordinate = self.coordinate.replace({'nr': index})
        self.bonds = GromacsSys.reindex_list(self.bonds, index, 2)
        self.pairs = GromacsSys.reindex_list(self.pairs, index, 2)
        self.angles = GromacsSys.reindex_list(self.angles, index, 3)
        self.dihedrals = GromacsSys.reindex_list(self.dihedrals, index, 4)
        self.impropers = GromacsSys.reindex_list(self.impropers, index, 4)
        self.amine = GromacsSys.reindex_list(self.amine, index, 3)
        self.carboxyl = GromacsSys.reindex_list(self.carboxyl, index, 3)
        self.charge_shift = GromacsSys.reindex_list(self.charge_shift, index,
                                                    7)
        self.Graph = nx.relabel_nodes(self.Graph, index)
        logging.info("Atom indices reset.")

    def output_contents(self,
                        new_top='',
                        new_gro='',
                        file=True,
                        graph=False,
                        change_files=False):
        """Outputs the contents of the system.

        Keyword Arguments:
            new_top {string} -- [Filename. If given, outputs topology data to new_top.
            Otherwise outputs to self.topfile.] (default: {''})
            new_gro {string} -- [Filename. If given, outputs topology data to new_gro.
            Otherwise outputs to self.grofile.] (default: {''})
            file {bool} -- [If True, outputs topology data and coordinate data to self.topfile and
            self.grofile respectively. If false, output to terminal.] (default: {True})
            graph {bool} -- [If True, draws the molecule topology graph with networkx.draw.
            The graph will be saved to graph_[num].png if file is True, or shown directly if otherwise.]
            (default: {False})
            change_files (bool): [If both this and file are true, use the newly output files as system files.]
            Defaults to False.
        """
        logging.info("Running output_contents")
        if file:
            # Decide output file based on whether arguments new_top and
            # new_gro are given.
            if new_top:
                out_topfile = new_top
            else:
                out_topfile = self.topfile
            if new_gro:
                out_grofile = new_gro
            else:
                out_grofile = self.grofile
            # Check if topfile and grofile exists. If they do, change their names:
            if os.path.isfile(out_topfile):
                logging.info(
                    "Topology file with the same name already exists.")
                n = 1
                while True:
                    fn = out_topfile + '_' + str(n)
                    if not os.path.isfile(fn):
                        break
                    else:
                        n = n + 1
                os.rename(out_topfile, fn)
                logging.info(f"Renamed old topology file to {fn}")
            if os.path.isfile(out_grofile):
                logging.info(
                    "Coordinate file with the same name already exists.")
                n = 1
                while True:
                    fn = out_grofile + '_' + str(n)
                    if not os.path.isfile(fn):
                        break
                    else:
                        n = n + 1
                os.rename(out_grofile, fn)
                logging.info(f"Renamed old coordinate file to {fn}")
            # Write topology file:
            with open(out_topfile, 'w+') as topf:
                print(
                    ";\n;Topology created by atuomatic script\ngrafting TMAO onto TMC.\n",
                    "\n; Include forcefield parameters\n#include \"./charmm36-mar2019.ff/forcefield.itp\"\n\n",
                    "[ moleculetype ]\nOther               3\n[ atoms ]\n",
                    file=topf)
                for i in self.atoms.values.tolist():
                    print(
                        "{:d}\t{}\t{:d}\t{}\t{}\t{:d}\t{:>6.3f}\t{:>6.3f} {}{:d}\t{}"
                        .format(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7],
                                i[8], i[9], i[10]),
                        file=topf)
                print('\n[ bonds ]\n', file=topf)
                for i in self.bonds:
                    print(*i, sep='    ', file=topf)
                print('\n[ pairs ]\n', file=topf)
                for i in self.pairs:
                    print(*i, sep='    ', file=topf)
                print('\n[ angles ]\n', file=topf)
                for i in self.angles:
                    print(*i, sep='    ', file=topf)
                print('\n[ dihedrals ]\n', file=topf)
                for i in self.dihedrals:
                    print(*i, sep='    ', file=topf)
                if self.impropers:
                    print('\n[ dihedrals ]\n', file=topf)
                    for i in self.impropers:
                        print(*i, sep='    ', file=topf)
                if self.cmap:
                    print('\n[ cmap ]\n', file=topf)
                    for i in self.cmap:
                        print(*i, sep='    ', file=topf)
                print(
                    "\n; Include Position restraint file\n#ifdef POSRES\n#include \"posre.itp\"\n#endif\n\n",
                    "; Include water topology\n#include \"./charmm36-mar2019.ff/tip3p.itp\"\n\n#ifdef POSRES_WATER\n",
                    "; Position restraint for each water oxygen\n[ position_restraints ]\n",
                    ";  i funct       fcx        fcy        fcz\n    1    1       1000       1000       1000\n#endif",
                    "\n\n[ system ]\nTestCrosslinkSys\n\n[ molecules ]\nOther               1\n",
                    file=topf)
            logging.info(f"Written topology file {out_topfile}")
            # Write coordinate file:
            with open(out_grofile, 'w+') as grof:
                print('Coordinate file created by script.', file=grof)
                print(' ', self.atomnum, file=grof)
                if self.velocity:
                    for i in self.coordinate.values.tolist():
                        print(
                            "{:>5}{:5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.4f}{:>8.4f}{:>8.4f}"
                            .format(i[0], i[1], i[2], i[3], i[4], i[5], i[6],
                                    i[7], i[8], i[9]),
                            file=grof)
                else:
                    for i in self.coordinate.values.tolist():
                        print("{:>5}{:5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}".
                              format(i[0], i[1], i[2], i[3], i[4], i[5], i[6]),
                              file=grof)
                print('   {:.5f}  {:.5f}   {:.5f}'.format(
                    self.box[0], self.box[1], self.box[2]),
                      file=grof)
            if change_files:
                self.topfile = out_topfile
                self.grofile = out_grofile
            if graph and file:
                logging.info(f"Written coordinate file {out_grofile}")
            else:
                logging.info(f"Written coordinate file {out_grofile}\n")
        else:
            # Write to terminal:
            print('[ atoms ]\n')
            print(self.atoms.to_string(index=False, header=False))
            print('\n[ bonds]\n')
            for i in self.bonds:
                print(*i, sep='    ')
            print('\n[ pairs ]\n')
            for i in self.pairs:
                print(*i, sep='    ')
            print('\n[ angles ]\n')
            for i in self.angles:
                print(*i, sep='    ')
            print('\n[ dihedrals ]\n')
            for i in self.dihedrals:
                print(*i, sep='    ')
            print('\n[ dihedrals ]\n')
            for i in self.impropers:
                print(*i, sep='    ')
            print('Current time t= ', self.time, '\n')
            print(self.atomnum, '\n')
            print(self.coordinate.to_string(index=False, header=False))
            print(*self.box, sep=' ')
        # Draw graph:
        if graph:
            if file:
                nx.draw(self.Graph,
                        with_labels=True,
                        node_color='cyan',
                        node_size=60,
                        font_size=8)
                fn = 'graph_1.png'
                if os.path.isfile(fn):
                    n = 2
                    while True:
                        fn = 'grpah_' + str(n) + '.png'
                        if not os.path.isfile(fn):
                            break
                        else:
                            n = n + 1
                plt.savefig(fn, format='PNG')
                logging.info(f"Output graph file {fn}\n")
            else:
                nx.draw(self.Graph,
                        with_labels=True,
                        node_color='cyan',
                        node_size=60,
                        font_size=8)
                plt.show()

    @staticmethod
    def reindex_list(list, index_dic, d):
        """
        Find and replace the indices in a 2d list.
        When atom indices are changed, this method is called to update the indices
        in the lists of bonds, pairs, angles etc.

        Arguments:
            list {list} -- [The 2d list to be worked on.]
            index_dic {dictionary} -- [A dictionary of {"Find this value" : "Replace with this value"}]
            d {integer} -- [For each sub-list, replacce values in [0:d] while ignoring [d:None].
            This is used so that only the atom indices in each entries are changed, and the force field
            parameters coming after them are left unchanged.
            For bonds and pairs, d=2; for angels, d=3; for dihedrals and improper dihedrals, d=4.]

        Returns:
            [list] -- [list after replacement]
        """
        new_list = []
        for i in list:
            new_list.append(
                [index_dic[n] if n in index_dic.keys() else n
                 for n in i[0:d]] + i[d:None])
        # print(list)
        return new_list

    def add_bond(self, newbond):
        """
        Add an amide bond in the two atoms in argument newbond, removing one hydrogen
        from the amine group and one hydroxyl group from the carboxyl group.
        Then, the script makes changes to atom indices, atom types and charges.
        The script also finds all the pairs, angles and dihedrals that are formed by the
        creation of new bond from the networkx graph of the system, then adds these entries
        to the corresponding lists.

        Arguments:
            newbond {list} -- [Indices of atoms to form bond between:
            [Atom A (carboxyl carbon), Atom B (amine nitroge)] Note the order.
        """
        logging.info(f'Running add_bond {newbond}')
        # Check atom types.
        atomA = newbond[0]
        atomB = newbond[1]
        if (self.atoms.loc[self.atoms['nr'] == atomA]['type'].tolist()[0] != 'CG2O2') or\
                (self.atoms.loc[self.atoms['nr'] == atomB]['type'].tolist()[0] != 'NG2S3'):
            logging.error(
                "Argument newbond should [carboxyl C, nitrogen N] list.")
            raise ValueError(
                "Argument newbond should [carboxyl C, nitrogen N] list.")
        # Change residue names according to the number of amide bonds on them.
        nr_A = self.atoms.loc[self.atoms['nr'] == atomA]['resnr'].tolist()[0]
        nr_B = self.atoms.loc[self.atoms['nr'] == atomB]['resnr'].tolist()[0]
        if self.atoms.loc[self.atoms['nr'] ==
                          atomA]['resname'].tolist()[0] == 'TMA':
            restype_a = 'TMAm'
        elif self.atoms.loc[self.atoms['nr'] ==
                            atomA]['resname'].tolist()[0] == 'TMAm':
            restype_a = 'TMAd'
        elif self.atoms.loc[self.atoms['nr'] ==
                            atomA]['resname'].tolist()[0] == 'TMAd':
            restype_a = 'TMAt'
        self.atoms.loc[self.atoms['resnr'] == nr_A, 'resname'] = restype_a
        self.coordinate.loc[self.coordinate['resnr'] == nr_A,
                            'resname'] = restype_a
        if self.atoms.loc[self.atoms['nr'] ==
                          atomB]['resname'].tolist()[0] == 'MPD':
            restype_b = 'MPDm'
        elif self.atoms.loc[self.atoms['nr'] ==
                            atomB]['resname'].tolist()[0] == 'MPDm':
            restype_b = 'MPDd'
        self.atoms.loc[self.atoms['resnr'] == nr_B, 'resname'] = restype_b
        self.coordinate.loc[self.coordinate['resnr'] == nr_B,
                            'resname'] = restype_b
        # Change atom types of the amide bond C and N.
        self.atoms.loc[self.atoms['nr'] == atomA, 'type'] = 'CG2O1'
        self.atoms.loc[self.atoms['nr'] == atomB, 'type'] = 'NG2S1'
        self.atoms.loc[self.atoms['nr'] == atomA, 'marks'] = 'amide_C'
        self.atoms.loc[self.atoms['nr'] == atomB, 'marks'] = 'amide_N'
        H_change = [i[2] for i in self.amine if i[0] == atomB][0]
        logging.info(f'Amine hydrogen {H_change} is kept.')
        self.atoms.loc[self.atoms['nr'] == H_change, 'type'] = 'HGP1'
        # Change atomic charges.
        charge_group_A = [i for i in self.charge_shift if i[0] == atomA][0]
        charge_group_B = [i for i in self.charge_shift if i[0] == atomB][0]
        for i in range(len(charge_group_A)):
            self.atoms.loc[self.atoms['nr'] == charge_group_A[i],
                           'charge'] += self.delta_charge_C[i]
        for i in range(len(charge_group_B)):
            self.atoms.loc[self.atoms['nr'] == charge_group_B[i],
                           'charge'] += self.delta_charge_N[i]
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[0], 'charge'] = 0.596
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[1],
        #                'charge'] = -0.471
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[2],
        #                'charge'] = -0.134
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[3],
        #                'charge'] += 0.018
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[4],
        #                'charge'] += 0.018
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[5],
        #                'charge'] += 0.001
        # self.atoms.loc[self.atoms['nr'] == charge_group_A[6],
        #                'charge'] += 0.001
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[0],
        #                'charge'] = -0.508
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[1], 'charge'] = 0.339
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[2], 'charge'] = 0.151
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[3],
        #                'charge'] -= 0.006
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[4],
        #                'charge'] -= 0.006
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[5],
        #                'charge'] -= 0.002
        # self.atoms.loc[self.atoms['nr'] == charge_group_B[6],
        #                'charge'] -= 0.002
        # Delete the amine H and carboxyl OH removed when forming amide bond.
        # Atom indices are changed here.
        remove_list = [i[1:3] for i in self.carboxyl if i[0] == atomA
                       ][0] + [i[1] for i in self.amine if i[0] == atomB]
        self.remove_atoms(remove_list)
        # Changes variable atomA and atomB to the new indices of the corresponding atoms.
        orig_atomA = atomA
        orig_atomB = atomB
        atomA = atomA - len([i for i in remove_list if i < atomA])
        atomB = atomB - len([i for i in remove_list if i < atomB])
        # Change the indices in self.freeze_group.
        logging.info("Editing atom indices in group [ frozen ].")
        if self.freeze_group:
            self.freeze_group = [
                i - sum(1 for j in remove_list if j < i)
                for i in self.freeze_group
            ]
        # Add the new bond in self.bonds and self.Graph.
        self.bonds.append([atomA, atomB, 1])
        self.Graph.add_edge(atomA, atomB)
        # Make a subgraph U consisting of the atoms (nodes) nearby the new bond.
        # All nodes in U are no more then 2 edges away from either A or B.
        U = self.Graph.subgraph(
            set(nx.single_source_shortest_path(
                self.Graph, atomA, 2).keys()).union(
                    set(
                        nx.single_source_shortest_path(self.Graph, atomB,
                                                       2).keys())))
        # Find all paths consisting of less then 4 nodes (3 edges) in the graph U. If the beginning
        # and ending nodes of a path are from different residues, then it is a path
        # newly created after adding the bond. Note that this might not be true if you are not
        # crosslinking aromatic polymers, but we don't consider about those situations here.
        paths = []
        new_angles = []
        new_dihedrals = []
        new_pairs = []
        nodes_unvisited = list(U.nodes())
        # Iterate all the starting node - ending node combinations.
        for i in list(U.nodes()):
            nodes_unvisited.remove(i)
            for j in nodes_unvisited:
                # If starting node i and ending node j are from different residues, then
                # paths connecting them are related with the new bond. Find all paths connecting
                # them with length less then 3 edges.
                if self.atoms.loc[self.atoms['nr'] == i]['resname'].tolist()[0] is not\
                        self.atoms.loc[self.atoms['nr'] == j]['resname'].tolist()[0]:
                    paths = paths + list(nx.all_simple_paths(U, i, j, 3))
        # Among all the paths found here, those with a length of 2 edges indicate angles;
        # those with a length of 3 edges indicate dihedrals; the starting and ending nodes of
        # dihedrals make up pairs.
        for p in paths:
            if len(p) == 3:
                new_angles.append(p)
            elif len(p) == 4:
                new_dihedrals.append(p)
                new_pairs.append([p[0], p[-1]])
        # Add funct parameter to pair, angle and dihedral entries.
        self.pairs = self.pairs + [x + [1] for x in new_pairs]
        self.angles = self.angles + [x + [5] for x in new_angles]
        self.dihedrals = self.dihedrals + [x + [9] for x in new_dihedrals]
        # Add improper dihedral entry.
        impr_bcd = list(nx.neighbors(U, atomA))
        # Arrange the latter 3 atoms in the improper dihedral in NG2S1 - CG2R61 - OG2D1 order
        order = {"NG2S1": 1, "CG2R61": 2, "OG2D1": 3}
        impr_types = list(
            self.atoms.loc[self.atoms["nr"].isin(impr_bcd)]["type"].values)
        impr_nco = [
            x for y, x in sorted(zip([order[z] for z in impr_types], impr_bcd))
        ]
        new_improper = [atomA] + impr_nco + [2, 0, 1004.160]
        self.impropers.append(new_improper)
        self.cl_num = self.cl_num + 1
        # Write log.
        logging.info(
            f'Bond created between atoms {atomA} and {atomB} (formerly {orig_atomA} and {orig_atomB})'
        )
        logging.info(f'Residues {nr_A} and {nr_B} are connected.')
        if self.residue_connectivity:
            self.Residuenet.add_edge(nr_A, nr_B)
            if self.xlink_list_file:
                with open(self.xlink_list_file, mode='a') as f:
                    f.write(f"{nr_A} {restype_a} {nr_B} {restype_b}\n")
        # logging.info(f'Found new pairs: {" ".join(str(x) for x in new_pairs)}')
        # logging.info(f'Found new angles: {" ".join(str(x) for x in new_angles)}')
        # logging.info(f'Found new dihedrals: {" ".join(str(x) for x in new_dihedrals)}')
        # logging.info(f'Found new improper dihedrals: {" ".join(str(x) for x in new_improper)}\n')

    def gmx_run(self,
                operation,
                mdpfile,
                tpr_name='',
                use_ndx=True,
                verbose=False):
        """ Call Gromacs to make an MD run, using self.topfile and self.grofile as input files.
            Performs energy minimization when operation = opt, with the following command:
            grompp -f [mdpfile -c [grofile] -p [topfile] -o [tpr_name.tpr]
            mdrun -deffnm [tpr_name]
            Performs MD run when operation = md, with the following command:
            grompp -f [mdpfile] -c [grofile] -r [grofile] -p [topfile] -o [tpr_name.tpr]
            mdrun -deffnm [tpr_name]
            If tpr_name is given, this would be the name of files output by grompp and mdrun.
            Otherwise, it would use self.tprfile.

            The script would check whether files with the same name exist, and if they do,
            the script would add number postfix to avoid this. Note that as a result of this,
            all files that uses tpr_name in their names would have their names changes.

            NOTE: This is wrote for older versions (4.6.5) of Gromacs. To use this script with
            newer versions of Gromacs, make changes to the commands.
        Arguments:
            operation  {string} -- [Commands to be run 'opt' - energy minimization, 'md' - MD run]
            mdpfile {string} -- [Output path and filename of mdp file to be]
            tpr_name {string} -- [Name of output files, without extension (passed to argument -deffnm)]
            verbose {bool, optional} -- [If True, run mdrun verbosely (with -v argument)] Default False
        """
        logging.info(f"Running gmx_run to do {operation}")
        # Check duplicate file names
        if not tpr_name:
            tpr_name = self.tprfile
        if os.path.isfile(tpr_name + '.tpr'):
            n = 1
            tpr_pfx = tpr_name
            while True:
                tpr_name = tpr_pfx + '_' + str(n)
                if not os.path.isfile(tpr_name + '.tpr'):
                    break
                else:
                    n = n + 1
            logging.info(
                f"Tpr file already existed. Changed tpr_name to {tpr_name}")
        # Run the commands:
        if operation == 'opt':
            cmd = f"grompp -f {mdpfile} -c {self.grofile} -p {self.topfile} "
            if self.ndxfile and use_ndx:
                cmd = cmd + f"-n {self.ndxfile} -o {tpr_name}.tpr"
            else:
                cmd = cmd + f"-o {tpr_name}.tpr"
            if not self.old_gmx:
                cmd = "gmx " + cmd
            logging.info("Running geometry optimization: " + cmd)
            sp.check_call(cmd, shell=True)
            if verbose:
                if not self.old_gmx:
                    logging.info(
                        f"Running geometry optimization: gmx mdrun -v -deffnm {tpr_name}"
                    )
                    sp.check_call('gmx mdrun -v -deffnm ' + tpr_name,
                                  shell=True)
                else:
                    logging.info(
                        f"Running geometry optimization: mdrun -v -deffnm {tpr_name}"
                    )
                    sp.check_call('mdrun -v -deffnm ' + tpr_name, shell=True)
            else:
                if not self.old_gmx:
                    logging.info(
                        f"Running geometry optimization: gmx mdrun -deffnm {tpr_name}"
                    )
                    sp.check_call('gmx mdrun -deffnm ' + tpr_name, shell=True)
                else:
                    logging.info(
                        f"Running geometry optimization: mdrun -deffnm {tpr_name}"
                    )
                    sp.check_call('mdrun -deffnm ' + tpr_name, shell=True)
            self.grofile = tpr_name + '.gro'
            logging.info(f"New gro file: {self.grofile}")
            self.tprfile = tpr_name + '.tpr'
            logging.info(f"New tpr file: {self.tprfile}")
            self.read_gro()
        elif operation == 'md':
            cmd = f"grompp -f {mdpfile} -c {self.grofile} -p {self.topfile} "
            if self.ndxfile and use_ndx:
                cmd = cmd + f"-n {self.ndxfile} -o {tpr_name}.tpr"
            else:
                cmd = cmd + f"-o {tpr_name}.tpr"
            if not self.old_gmx:
                cmd = "gmx " + cmd
            logging.info("Running MD simulation: " + cmd)
            sp.check_call(cmd, shell=True)
            if verbose:
                if not self.old_gmx:
                    logging.info(
                        f"Running MD simulation: gmx mdrun -v -deffnm {tpr_name}"
                    )
                    sp.check_call('gmx mdrun -v -deffnm ' + tpr_name,
                                  shell=True)
                else:
                    logging.info(
                        f"Running MD simulation: mdrun -v -deffnm {tpr_name}")
                    sp.check_call('mdrun -v -deffnm ' + tpr_name, shell=True)
            else:
                if not self.old_gmx:
                    logging.info(
                        f"Running MD simulation: gmx mdrun -deffnm {tpr_name}")
                    sp.check_call('gmx mdrun -deffnm ' + tpr_name, shell=True)
                else:
                    logging.info(
                        f"Running MD simulation: mdrun -deffnm {tpr_name}")
                    sp.check_call('mdrun -deffnm ' + tpr_name, shell=True)
            self.grofile = tpr_name + '.gro'
            logging.info(f"New gro file: {self.grofile}")
            self.tprfile = tpr_name + '.tpr'
            logging.info(f"New tpr file: {self.tprfile}")
            self.read_gro()
        else:
            logging.error(("Argument Operation must be either 'opt' or 'md'."))
            raise Exception("Argument Operation must be either 'opt' or 'md'.")

    def remove_residue(self, res_to_remove):
        """Find residues with names matching res_to_remove and remove them.

        Args:
            res_to_remove ([String]): [Name of residues to be removed.]
        """
        logging.info(
            f"\nRunning remove_residue to delete {res_to_remove} residues")
        res_atoms = self.atoms.loc[self.atoms['resname'] ==
                                   res_to_remove]['nr'].values.tolist()
        res_residues = list(
            set(self.atoms.loc[self.atoms['resname'] == res_to_remove]
                ['resnr']))
        logging.info(f"Found {len(res_residues)} {res_to_remove} residues.")
        # logging.info(f"Residue indices: {res_residues}")
        # logging.info(f"These residues are included: {res_atoms}")
        self.remove_atoms(res_atoms)
        # Change residue numbers
        self.resnum = self.resnum - len(res_residues)
        resnr_old = list(set(self.atoms['resnr']))
        resnr_new = list(range(1, 1 + self.resnum))
        resnr_index = zip(resnr_old, resnr_new)
        self.atoms = self.atoms.replace({'resnr': resnr_index})
        self.coordinate = self.coordinate.replace({'resnr': resnr_index})
        logging.info("Residue numbers changed.")
        self.Residuenet.remove_nodes_from(res_residues)
        return len(res_residues)

    def atom_count(self,
                   dpc_bonds=True,
                   dpc_bonds_appr=True,
                   dpc_NO_ratio=True,
                   dpc_OO_ratio=True,
                   dpc_carboxyl_ratio=True,
                   dpc_TMAt=True,
                   exclude_freeze=True,
                   log_details=False):
        """Count the number of atoms and residues in the system, and calculate DPC.
        There are many different ways to calculate DPC in a polyamide system.
        The script can calculate them all, but by default,
        only dpc_TMAt is used for tracking crosslinking process.
        DPC(bonds) = (Number of amide bonds formed) / (Number of max possible amide bonds)
         = (Total nitrogen - amine nitrogen)/(Total nitrogen + carboxyl carbon)
        DPC(bonds, appropriate) = (Total nitrogen) / (Total nitrogen + carboxyl carbon)
        DPC(N/O ratio) = (4 * (Nitrogen/Oxygen ratio) - 2) / (1 + (Nitrogen/Oxygen ratio))
        DPC(O/O ratio) = (Number of hydroxyl oxygen) / (Number of carbonyl oxygen)
        DPC(Carboxyl_ratio) = (Number of amide bonds) / (Number of amide bonds and carboxyl groups)
        DPC(TMAt) = (Number of TMA linked with 3 MDP residues) / (Number of TMA residues)
        Args:
            exclude_freeze {bool}: [Wheter to ignore [ frozen ] atoms] Default True.
            log_details {bool}: [Whether to output atom numbers, residue numbers and DPC values to log file.]
                Default False.

        Returns:
            [Pandas Series]: [Number of each types of atoms]
            [Pandas Series]: [Number of each types of residues]
            [Dictionary]: [Dictionary of DPC values: {DPC type: value}]
        """
        dpc = {}
        if exclude_freeze:
            count_atom = self.atoms[~self.atoms['nr'].isin(self.freeze_group)][
                'type'].value_counts()
        else:
            count_atom = self.atoms['type'].value_counts()
        for name in ['CG2O1', 'NG2S1', 'CG2O2', 'NG2S3', 'OG311']:
            if name not in count_atom:
                count_atom[name] = 0
        if exclude_freeze:
            count_res = self.atoms[~self.atoms['nr'].isin(self.freeze_group)][
                'resname'].value_counts()
        else:
            count_res = self.atoms['resname'].value_counts()
        for name in ['TMA', 'TMAm', 'TMAd', 'TMAt', 'MPD', 'MPDm', 'MPDd']:
            if name not in count_res:
                count_res[name] = 0
        res_atomn = pd.Series(
            [16, 15, 14, 21, 19, 17, 15],
            index=['MPD', 'MPDm', 'MPDd', 'TMA', 'TMAm', 'TMAd', 'TMAt'])
        count_res = (count_res.div(res_atomn)).astype(int)
        if dpc_bonds:
            dpc.update({
                "DPC_bonds": (count_atom['NG2S1'] /
                              (count_atom['NG2S1'] + count_atom['NG2S3'] +
                               count_atom['CG2O2']))
            })
        if dpc_bonds_appr:
            dpc.update({
                "DPC_bonds_approximated":
                (count_atom['NG2S1'] + count_atom['NG2S3']) /
                (count_atom['NG2S1'] + count_atom['NG2S3'] +
                 count_atom['CG2O2'])
            })
        if dpc_NO_ratio:
            no_ratio = (count_atom['NG2S1'] + count_atom['NG2S3']) / (
                count_atom['OG2D1'] + count_atom['OG311'])
            dpc.update({
                "DPC_crosslinked_portion(NO_ratio)":
                (4 * no_ratio - 2) / (1 + no_ratio)
            })
        if dpc_OO_ratio:
            oo_ratio = count_atom['OG311'] / count_atom['OG2D1']
            dpc.update(
                {"DPC_crosslinked_portion(Oxygen_ratio)": 1 - 3 * oo_ratio})
        if dpc_carboxyl_ratio:
            cx_ratio = count_atom['CG2O1'] / (count_atom['CG2O2'] +
                                              count_atom['CG2O1'])
            dpc.update({"DPC_carboxyl_ratio": cx_ratio})
        if dpc_TMAt:
            dpc.update({
                "DPC_tri-linked-TMA_portion":
                count_res['TMAt'] /
                (count_res['TMAm'] + count_res['TMAd'] + count_res['TMAt'])
            })
        if log_details:
            for name in ['MPD', 'MPDm', 'MPDd', 'TMA', 'TMAm', 'TMAd', 'TMAt']:
                if name not in count_res:
                    count_res[name] = 0

            total_C = count_atom['CG2R61'] + count_atom['CG2O2'] + count_atom[
                'CG2O1']
            total_N = count_atom['NG2S1'] + count_atom['NG2S3']
            total_H = count_atom['HGP1'] + count_atom['HGP4'] + count_atom[
                'HGR61']
            total_O = count_atom['OG2D1'] + count_atom['OG311']
            mass_atoms = [
                total_C * 12.011, total_N * 14.007, total_H * 1.008,
                total_O * 15.999
            ]
            mass_fraction = [round(i / sum(mass_atoms), 3) for i in mass_atoms]

            logging.info(
                f"The numbers of C, N, H, O atoms are {total_C}, {total_N}, {total_H}, {total_O}\n"
            )
            logging.info(f"Their mass fractions are {mass_fraction}\n")
            logging.info(f"Atom numbers:\n{count_atom}")
            logging.info(f"\nNumbers of residues:\n{count_res}")
            logging.info(f"\nDegree of polymer crosslinking:\n{dpc}\n")
            return dpc, count_atom, count_res
        else:
            return dpc

    def CN_dist(self,
                r,
                use_zlim=False,
                zlim=[0, 0],
                use_xymargin=False,
                xymargin=0,
                nofreeze=True,
                wrap_PBC=True,
                output=''):
        """Running CN_dist to find carboxyl group (according to the coordinate of carboxyl C) and amine group
        (according to the coordinate of amine N) whose distance is less than a certain value of  {r}.
        Principle: Calculate their distance based on the coordinate of carboxyl C and amine N in the gro file.
        Steps as follows:
          1. Filter out the numbers and coordinates of the carboxyl C and amine N atoms
          and put them in the C_xyz and N_xyz lists respectively;
          2. For periodic reasons, the atoms that not in box are adjusted back into a box.
          3. According to the value of the x and y coordinates of each carboxyl C in the C_xyz list,
             divide the region into 9 regions.
          4. For the list of carboxyl C atoms in each region, find the list of amine N {N_points}
           whose x, y, z coordinates are at the distance r from the carboxyl C atom,
           then calculate the distance between the carboxyl C and  {N_points} respectively,
           the distances less than {r} are saved in the {list_dist} list.
        Arguments zlim and xymargin can be used to limit crosslinking in a certain range, that is,
        atoms outsides of this range is ignored during crosslinking.

        Arguments:
            r-- {float} -- the max distance between a C-N group.
            use_zlim (bool, optional): Whether to limit crosslinking in range zlim[0]<z<zlim[0]. Defaults to False.
            zlim (list, optional): Z range of crosslinking, used only when use_zlim=True. Defaults to [0, box[2]].
            use_xymargin (bool, optional): Whether to limit crosslinking in range xymargin < x < box[0] - xymargin
                and xymargin < y < box[1] - xymargin. Defaults to False.
            xymargin (int, optional): XY margin of crosslink-allowed region. Defaults to 0.
            nofreeze (bool, optional): If true, exclude atoms in self.freeze_group. Defaults to True.
            wrap_PBC (bool, optional): If true, wrap all atoms inside main PBC cell in XY directions.
                Defaults to True. NOTE: does not wrap in Z direction. Also wrapping is applied before xymargin.
            output (str, optional): A file name can be given here. List of found pairs are output to file.
                Defaults to ''.

        Returns:
            If C-N pairs with distance less then dist are found, returns a list of {list_dist}.
            Otherwise return an empty list.

            list_dist -- {list} -- the list of C-N  pairs and their distance.
            [[carboxyl C, amine N, distance of C-N]]
        """
        logging.info(f"Running CN_dist with a minimal distance of {r}.")
        # list_xyz = self.coordinate[['nr', 'x', 'y', 'z']].values.tolist()
        # carboxyl_C = [i[0] for i in self.carboxyl]
        # amine_N = [i[0] for i in self.amine]

        # Filter out the numbers and coordinates containing carboxyl C and amine N atoms
        # and put them in the C_xyz and N_xyz lists respectively
        # C_xyz -- {list} -- List of carboxyl C coordinates， [[Number of carboxyl C, x, y ,z ]]
        # N_xyz -- {list} -- List of amine N coordinates,  [[Number of amine N, x, y ,z ]]
        # C_xyz = []
        # N_xyz = []
        C_xyz = self.coordinate[self.coordinate['nr'].isin(
            [i[0] for i in self.carboxyl])][['nr', 'x', 'y',
                                             'z']].values.tolist()
        N_xyz = self.coordinate[self.coordinate['nr'].isin(
            [i[0] for i in self.amine])][['nr', 'x', 'y',
                                          'z']].values.tolist()
        if use_zlim:
            C_xyz = [i for i in C_xyz if (zlim[0] < i[3] < zlim[1])]
            N_xyz = [i for i in N_xyz if (zlim[0] < i[3] < zlim[1])]
        # for j in list_xyz:
        #     if j[0] in carboxyl_C:
        #         C_xyz.append(j)
        #     if j[0] in amine_N:
        #         N_xyz.append(j)
        if nofreeze and self.freeze_group:
            C_xyz = [i for i in C_xyz if i[0] not in self.freeze_group]
            N_xyz = [i for i in N_xyz if i[0] not in self.freeze_group]
        # For periodic reasons, the atoms that not in box are adjusted back into a box.
        if wrap_PBC:
            C_xyz = [[a, b % self.box[0], c % self.box[1], d]
                     for a, b, c, d in C_xyz]
            N_xyz = [[a, b % self.box[0], c % self.box[1], d]
                     for a, b, c, d in N_xyz]
        # TODO: warp PBC and cut xy margin: which one comes first?
        if use_xymargin:
            C_xyz = [
                i for i in C_xyz
                if (xymargin < i[1] < self.box[0] - xymargin) and (
                    xymargin < i[2] < self.box[1] - xymargin)
            ]
            N_xyz = [
                i for i in N_xyz
                if (xymargin < i[1] < self.box[0] - xymargin) and (
                    xymargin < i[2] < self.box[1] - xymargin)
            ]
        # According to the value of the x and y coordinates of each carboxyl C in the C_xyz list,
        # divide the region into 9 regions.
        # C_Center，X_Left，X_Right，Y_Left，Y_Right
        # are the center area and the four sides of the xy plane, respectively.
        C_Center = []
        X_Left = []
        X_Right = []
        Y_Left = []
        Y_Right = []
        # XX_Left，XY_Left，XY_Right，YY_Right are the four vertex positions in the xy plane.
        XX_Left = []
        XY_Left = []
        XY_Right = []
        YY_Right = []
        for i in C_xyz:
            if r <= i[1] <= self.box[0] - r and r <= i[2] <= self.box[1] - r:
                C_Center.append(i)
            elif i[1] <= r and r <= i[2] <= self.box[1] - r:
                X_Left.append(i)
            elif i[1] >= self.box[0] - r and r <= i[2] <= self.box[1] - r:
                X_Right.append(i)
            elif i[2] <= r and r <= i[1] <= self.box[0] - r:
                Y_Left.append(i)
            elif i[2] >= self.box[1] - r and r <= i[1] <= self.box[0] - r:
                Y_Right.append(i)
            if i[1] < r and i[2] < r:
                XX_Left.append(i)
            elif i[1] < r and i[2] > self.box[1] - r:
                XY_Left.append(i)
            elif i[1] > self.box[0] - r and i[2] < r:
                XY_Right.append(i)
            elif i[1] > self.box[0] - r and i[2] > self.box[1] - r:
                YY_Right.append(i)
        Dict_C_list = {
            'C_Center': C_Center,
            'X_Left': X_Left,
            'X_Right': X_Right,
            'Y_Left': Y_Left,
            'Y_Right': Y_Right,
            'XX_Left': XX_Left,
            'XY_Left': XY_Left,
            'XY_Right': XY_Right,
            'YY_Right': YY_Right
        }
        list_dist = []
        for key, value in Dict_C_list.items():
            for i in value:
                C_point = np.array([i[1], i[2], i[3]])
                '''
                    N_points -- {list} --list of amine whose x, y, z coordinates
                    are at the distance r from the carboxyl C ,
                    [[Number of amine N, x, y ,z ]]
                    N_nr -- {list} --[Number of amine N]
                    Dict_C_list-- {dict} --Index of seven region and the carboxy C list of the corresponding region.
                '''
                N_points = []
                N_nr = []
                # Consider according to the key in the dictionary.
                if Dict_C_list[key] == C_Center:
                    for j in N_xyz:
                        if abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == X_Left:
                    for j in N_xyz:
                        # The box x coordinate size should be greater than 3r
                        # j[1]   the x coordinate of amine N
                        if j[1] >= self.box[0] - r:
                            if abs(j[1] - self.box[0] - i[1]) <= r and abs(
                                    j[2] - i[2]) <= r and abs(j[3] -
                                                              i[3]) <= r:
                                N_points.append(
                                    [j[1] - self.box[0], j[2], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == X_Right:
                    for j in N_xyz:
                        if j[1] <= r:
                            if abs(j[1] + self.box[0] - i[1]) <= r and abs(
                                    j[2] - i[2]) <= r and abs(j[3] -
                                                              i[3]) <= r:
                                N_points.append(
                                    [j[1] + self.box[0], j[2], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == Y_Left:
                    for j in N_xyz:
                        # The box y  coordinate size should be greater than 3r
                        # j[2]   the y coordinate of amine N
                        if j[2] >= self.box[1] - r:
                            if abs(j[1] - i[1]) <= r and abs(
                                    j[2] - self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append(
                                    [j[1], j[2] - self.box[1], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == Y_Right:
                    for j in N_xyz:
                        if j[2] <= r:
                            if abs(j[1] - i[1]) <= r and abs(
                                    j[2] + self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append(
                                    [j[1], j[2] + self.box[1], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == XX_Left:
                    for j in N_xyz:
                        if j[1] >= self.box[0] - r and j[2] >= self.box[1] - r:
                            if abs(j[1] - self.box[0] - i[1]) <= r and abs(
                                    j[2] - self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append([
                                    j[1] - self.box[0], j[2] - self.box[1],
                                    j[3]
                                ])
                                N_nr.append(j[0])
                        elif j[1] >= self.box[0] - r and j[2] <= 2 * r:
                            if abs(j[1] - self.box[0] - i[1]) <= r and abs(
                                    j[2] - i[2]) <= r and abs(j[3] -
                                                              i[3]) <= r:
                                N_points.append(
                                    [j[1] - self.box[0], j[2], j[3]])
                                N_nr.append(j[0])
                        elif j[1] <= 2 * r and j[2] >= self.box[1] - r:
                            if abs(j[1] - i[1]) <= r and abs(
                                    j[2] - self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append(
                                    [j[1], j[2] - self.box[1], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])

                elif Dict_C_list[key] == XY_Left:
                    for j in N_xyz:
                        if j[1] >= self.box[0] - r and j[2] <= r:
                            if abs(j[1] - self.box[0] - i[1]) <= r and abs(
                                    j[2] + self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append([
                                    j[1] - self.box[0], j[2] + self.box[1],
                                    j[3]
                                ])
                                N_nr.append(j[0])
                        elif j[1] <= 2 * r and j[2] <= r:
                            if abs(j[1] - i[1]) <= r and abs(
                                    j[2] + self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append(
                                    [j[1], j[2] + self.box[1], j[3]])
                                N_nr.append(j[0])
                        elif j[1] >= self.box[0] - r and j[
                                2] >= self.box[0] - 2 * r:
                            if abs(j[1] - self.box[0] - i[1]) <= r and abs(
                                    j[2] - i[2]) <= r and abs(j[3] -
                                                              i[3]) <= r:
                                N_points.append(
                                    [j[1] - self.box[0], j[2], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == XY_Right:
                    for j in N_xyz:
                        if j[1] < r and j[2] >= self.box[0] - r:
                            if abs(j[1] + self.box[0] - i[1]) <= r and abs(
                                    j[2] - self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append([
                                    j[1] + self.box[0], j[2] - self.box[1],
                                    j[3]
                                ])
                                N_nr.append(j[0])
                        elif j[1] <= r and j[2] <= 2 * r:
                            if abs(j[1] + self.box[0] - i[1]) <= r and abs(
                                    j[2] - i[2]) <= r and abs(j[3] -
                                                              i[3]) <= r:
                                N_points.append(
                                    [j[1] + self.box[0], j[2], j[3]])
                                N_nr.append(j[0])
                        elif j[1] >= self.box[0] - 2 * r and j[
                                2] >= self.box[0] - r:
                            if abs(j[1] - i[1]) <= r and abs(
                                    j[2] - self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append(
                                    [j[1], j[2] - self.box[1], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                elif Dict_C_list[key] == YY_Right:
                    for j in N_xyz:
                        if j[1] <= r and j[2] <= r:
                            if abs(j[1] + self.box[0] - i[1]) <= r and abs(
                                    j[2] + self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append([
                                    j[1] + self.box[0], j[2] + self.box[1],
                                    j[3]
                                ])
                                N_nr.append(j[0])
                        elif j[1] >= self.box[0] - 2 * r and j[2] <= r:
                            if abs(j[1] - i[1]) <= r and abs(
                                    j[2] + self.box[1] -
                                    i[2]) <= r and abs(j[3] - i[3]) <= r:
                                N_points.append(
                                    [j[1], j[2] + self.box[1], j[3]])
                                N_nr.append(j[0])
                        elif j[1] <= r and j[2] >= self.box[1] - 2 * r:
                            if abs(j[1] + self.box[0] - i[1]) <= r and abs(
                                    j[2] - i[2]) <= r and abs(j[3] -
                                                              i[3]) <= r:
                                N_points.append(
                                    [j[1] + self.box[0], j[2], j[3]])
                                N_nr.append(j[0])
                        elif abs(j[1] - i[1]) <= r and abs(
                                j[2] - i[2]) <= r and abs(j[3] - i[3]) <= r:
                            N_points.append([j[1], j[2], j[3]])
                            N_nr.append(j[0])
                if len(N_points) > 0:
                    # Calculate the Euclidean distance from one point to multiple points
                    distances = np.sqrt(
                        np.sum(np.asarray(C_point - N_points)**2, axis=1))
                    for k in range(len(distances)):
                        if distances[k] <= r:
                            list_dist.append([
                                int(i[0]),
                                int(N_nr[k]),
                                round(distances[k], 4)
                            ])
        # Delete the items that may appear repeatedly in the list of list_dist
        list_dist = [list(i) for i in set(tuple(j) for j in list_dist)]
        # In case break_coord is used, delete the items in which two atoms are on different sides
        # of the "barrier"
        if list_dist:
            # Sort according to the third element (C-N distance).
            list_dist_sorted = sorted(list_dist, key=lambda x: x[2])
            # dist_min = min([i[2] for i in list_dist])
            # atom_mindist = [i[0:2] for i in list_dist if i[2] == dist_min][0]
            if output:
                with open(output, 'w') as f:
                    for i in list_dist_sorted:
                        f.write(f'{i}\n')
            logging.info(
                f"Found {len(list_dist)} C-N pairs that meet crosslinking conditions."
            )
            logging.info(
                f"Closest pair being {list_dist_sorted[0][0:2]} with a distance of {list_dist_sorted[0][2]}\n"
            )
            # Check if this is a cross link betweeen different periodic images.
            cross_PBC = []
            coord_a = self.coordinate.loc[self.coordinate['nr'] ==
                                          list_dist_sorted[0][0]]
            coord_b = self.coordinate.loc[self.coordinate['nr'] ==
                                          list_dist_sorted[0][1]]
            if (abs(coord_a['x'].values[0] - coord_b['x'].values[0]) >
                    self.box[0] / 2):
                cross_PBC += [True]
            else:
                cross_PBC += [False]
            if (abs(coord_a['y'].values[0] - coord_b['y'].values[0]) >
                    self.box[1] / 2):
                cross_PBC += [True]
            else:
                cross_PBC += [False]
            return list_dist_sorted[0][0:2], round(list_dist_sorted[0][2],
                                                   3), cross_PBC
        else:
            logging.warning(
                "No C-N pair that meet crosslinking conditions found.\n")
            return [], 0, False

    def output_ndx(self, ndxname, usetpr=False):
        """
        Outputs ndx file corresponding to self.tprfile or self.grofile.
        First, call make_ndx to output an ndx file. Then, if self.freeze_group is not empty,
        then add a group [ frozen ] to the ndx file and write the indices there. Finally, self.ndxfile
        becomes the newly created ndx file.
        This is only used when:
        1) You are keeping a frozen "wall" to make sure that the polymer does not crosslink
        across a period boundary.
        2) The frozen "wall" is designated by giving their indices in the ndx file.
        (You can also give another name to the frozen residues in the top and gro files.
        In that case, you don't need the ndx file.)

        NOTE: This is wrote for older versions (4.6.5) of Gromacs. To use this script with
        newer versions of Gromacs, make changes to the commands. Also, using ndx files with
        other groups manually defined is not covered in the method.
        Args:
            ndxname (string): [Name of the ndx file output.]
            usetpr (bool, optional): [If true, use self.tprfile as input file (-f argument)
            of make_ndx; otherwise self.grofile.] Default True
        """
        logging.info(f"Calling output_ndx to generate ndx file{ndxname}")
        # Check for duplicate names.
        if os.path.isfile(ndxname):
            if os.path.isfile(ndxname + '_old'):
                logging.info('Delete old ndx file.')
                os.remove(ndxname + '_old')
            logging.info(f"Former index file renamed to {ndxname}_old")
            os.rename(ndxname, ndxname + '_old')
        # Call make_ndx to make an ndx file.
        if usetpr:
            inputfile = self.tprfile
        else:
            inputfile = self.grofile
        cmd = f"""make_ndx -f {inputfile} -o {ndxname} << EOF
        q
        EOF
        """
        if not self.old_gmx:
            cmd = "gmx " + cmd
        sp.check_call(cmd, shell=True)
        # If self.freeze_group is not empty, write its contents under a [ frozen ] group
        # in the ndx file. 15 indices each row.
        if self.freeze_group:
            logging.info("Writting [ frozen ]")
            with open(ndxname, mode='a') as f:
                f.write('\n[ frozen ]\n')
                for i in range(0, len(self.freeze_group) + 1, 15):
                    f.write(' '.join(
                        str(n) for n in self.freeze_group[i:(i + 15)]) + '\n')
        logging.info(
            f"ndx file {ndxname} written and will be used in subsequent run.")
        self.ndxfile = ndxname

    def coord_wrap(self):
        """Wrap all atom coordinates back inside PBC boundary.
        """
        self.coordinate['x'] = self.coordinate['x'] % self.box[0]
        self.coordinate['y'] = self.coordinate['y'] % self.box[1]
        self.coordinate['z'] = self.coordinate['z'] % self.box[2]

    def stop_using_ndx(self, remove_res=False):
        """Stops using frozen walls.

        Args:
            remove_res (bool, optional): If True remove all frozen residues.
            If False frozen residues will be kept and allowed to move or react. Defaults to False.
        """
        logging.info(
            "stop_using_ndx called. Index file with [ frozen ] will not be used."
        )
        self.ndxfile = ''
        if remove_res:
            logging.info("Frozen monomers will be removed.")
            self.remove_atoms(self.freeze_group)
            # NOTE: While this removes a whole set of residues and leave disconnected residue indices.
            # This seems to be OK
        else:
            logging.info("Frozen monomers are free to move and react.")
            # Add the residue nodes back
            logging.info("Mark atoms again:")
            self.mark_atoms()
            if self.residue_connectivity:
                self.Residuenet.add_nodes_from(self.freeze_res)

    def count_clusters(self):
        """Counts "clusters" of atoms connected by covalent bonds in the system, using netwokx's
        number_connected_components function. In case some monomers are frozen, the frozen monomers
        are not count.

        Returns:
            [int]: number of connected components.
        """
        if self.residue_connectivity:
            return nx.number_connected_components(self.Residuenet)
        else:
            G_counted = self.Graph.copy()
            if self.freeze_group:
                G_counted.remove_nodes_from(self.freeze_group)
            return nx.number_connected_components(G_counted)

    def output_residue_network(self, res_file):
        """Outputs residue network to res_file.

        Args:
            res_file (string): Path of output file
        """
        if not self.Residuenet:
            logging.error("Residue network information not kept!")
            raise ValueError("Residue network information not kept!")
        nx.write_gexf(self.Residuenet, res_file)

    def output_gro(self, groname):
        """Outputs a gro file only.
        The class will use the new gro file as self.grofile, while self.topfile remains the same.

        Args:
            groname (str): path of the new gro file.
        """
        with open(groname, 'w+') as grof:
            print('Coordinate file created by script.', file=grof)
            print(' ', self.atomnum, file=grof)
            if self.velocity:
                for i in self.coordinate.values.tolist():
                    print(
                        "{:>5}{:5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.4f}{:>8.4f}{:>8.4f}"
                        .format(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7],
                                i[8], i[9]),
                        file=grof)
            else:
                for i in self.coordinate.values.tolist():
                    print("{:>5}{:5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}".format(
                        i[0], i[1], i[2], i[3], i[4], i[5], i[6]),
                          file=grof)
            print('   {:.5f}  {:.5f}   {:.5f}'.format(self.box[0], self.box[1],
                                                      self.box[2]),
                  file=grof)
        self.grofile = groname

    def adjust_z(self,
                 zl,
                 zrange,
                 nwater,
                 filename,
                 mdp_min,
                 mdp_nvt,
                 dzstep=0.2):
        """Adjusts Z length of the box after liberating 1 water molecule,
        so that system density is stable.
        Specifically, this method calculates the total mass in the box before and after losing
        1 water molecule, and then reduces self.box[2] (Z length of the periodic box) so that system
        density stays the same. After reducing Z length, it will run an energy minimization and an NVT
        run (both supposed to be without freezeing residues), to keep the system equilibrated. If the
        method is required to reduce Z length by a large amount (i.e. nwater is big), it will do the
        reduction over multiple steps.
        Finally, itoutputs a new coordinate file for the adjusted system.
        Frozen walls are ignored during calculations of both Z length and tiotal mass.

        Args:
            zl (folat): Z length of the periodic box, minus thickness of the frozen layers.
            zrange (list): a list of two float values [zmin, zmax]. Only atoms in zmin<z<zmax
                are used in calculating total mass.
            nwater (int): number of water molecules liberated (mass decreased) since last adjust_z call.
            filename (string): name prefix of the new files generated during this process.
                For example, the new gro file will be filename_adjZ_1.gro
            dzstep (float): At most this value is reduced from Z length each time.
                Reduction Z length by a larger amount will be done over multiple steps.
            mdp_min (str): path to the mdp file for energy minimization.
            mdp_nvt (str): path to the mdp file for NVT run.

        Returns:
            (float): Z length reduced.
        """

        logging.info(
            f"Running Adjust_Z, adjusting Z length after libreating {nwater} water molecules"
        )
        if not zl:
            zl = self.box[2]
        # Calculate total mass
        freelayer = list(self.coordinate.loc[
            (self.coordinate["z"] > zrange[0])
            & (self.coordinate["z"] < zrange[1])]["nr"].values)
        total_mass = self.atoms.loc[self.atoms["nr"].isin(
            freelayer)]["mass"].sum()
        logging.info(
            f"Total mass in layer {zrange[0]} - {zrange[1]} is {total_mass}")
        # Calculate how much Z needs to be reduced.
        z_reduce_ratio = 18.015 * nwater / (total_mass + 18.015 * nwater)
        deltaz = zl * z_reduce_ratio
        # Begin reducing Z length.
        to_reduce = deltaz
        nstep = 1
        while (True):
            if to_reduce > dzstep:
                to_reduce = to_reduce - dzstep
                self.box[2] = self.box[2] - dzstep
            else:
                self.box[2] = self.box[2] - to_reduce
            stepfn = filename + str(nstep)
            self.output_gro(stepfn + ".gro")
            self.gmx_run('opt', mdp_min, stepfn + "_min")
            self.gmx_run('md', mdp_nvt, stepfn + "_nvt")
            if to_reduce < dzstep:
                break
            nstep = nstep + 1
        logging.info(
            f"Reduce Z to ratio {1 - z_reduce_ratio}, by {deltaz} nm, over {nstep} steps."
        )
        return deltaz
