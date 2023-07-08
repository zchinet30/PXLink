# PXLink

PXLink (Polymer-Cross-Link) is a Python script that uses GROMACS (version 4.6 or higher) to automatically perfom crosslinking simulation. Currently it can use trimesic acid and m-Phenylenediamine as monomers to generate a crosslinked aromatic polyamide membrane model.

A detailed description of this script and demonstrations of its usage can be found in the following manuscript:
Chi Zhang, Guangle Bu, Md Symon Jahan Sajib, Lida Meng, Shiying Xu, Size Zheng, Lin Zhang, Tao Wei, PXLink: A simulation program of polymer crosslinking to study of polyamide membrane,  Computer Physics Communications, 2023, doi.org/10.1016/j.cpc.2023.108840

# System requirements:
`
The script is programmed and tested in Linux environment. 
1.	Install GROMACS 4.6 or newer. For installation, please refer to GROMACS website:  https://www.gromacs.org/
2.	Install Python 3.9 or newer. 
3.  Install the required Python packages: Numpy, Pandas, Matplotlib and NetworkX using a Python package manager such as pip or anaconda. For example, it is possible to install these packages with pip using the following command:
```C
$python3 -m pip install Numpy Pandas NetworkX Matplotlib
```


# How to use:

To use the script, the user needs to:
1. Place the following files in the same folder:
    - The PXLink script (PXLink.py)
    - The execution script (e,g, example_run_script.py)
    - The force field directory charmm36-mar2019.ff (We use a modified CHARMM36 force field which can found in the Example_crosslinking directory. If you want to use other force fields you might need to modify the mark_atoms, add_bond, count_atoms and output_contents methods in the PXLink script. Make sure this force field is included in the topology file.)
    - The initial system (files describing the initial system: init.top, init.gro, optionally the index file init.ndx).
    -The mdp files needed for GROMACS.
2. Prepare the running script, designate the input files and assign desired values to the controlling variables. The default values in example_run_script.py would run crosslinking until its DPC reaches 30%. These values can be changed according to the user’s requirements. Also, if the user wants to change simulation conditions such as time and temperature, please change the corresponding parameters in the mdp files. Refer to the comments in example_run_script.py for more details.
3. Run the execution script: 
```C
$Python3 example_run_script.py
```
Note: it is also possible to restart a run from halfway by using "top", "gro" and optionally "ndx" files generated during previous simulations. In that case, the user needs to modify the progress record variables to the values they left off (check the log file for the values). For the demo execution script example_run_script.py, these variables can be found in the "Progress record variables" part.


# File list:

PXLink.py: PXLink script used to perform simulation crosslinking
<br><br>

Example_crosslinking: folder of a small system (240 TMA and 360 MPD residues in a 4.48×4.48×6.00 nm3 box) as an example of the crosslinking process.

|- charmm36-mar2019.ff
<br>
folder of a customized version of CHARMM36 force field specifically for the crosslinking simulation of TMA (trimesic acid) and MPD (m-Phenylenediamine). Extra entries are added to the files "merged.rtp" (which describe the two monomer residues) and "ffbonded.itp" (extra angle and dihedral parameters needed for simulating the polymer), such that GROMACS can recognize the two monomer and the corresponding angles and dihedrals.

|- init.top
<br>
GROMACS topology file of the initial (equilibrated monomer mix) system.

|- init.gro
<br>
GROMACS coordinate file of the initial (equilibrated monomer mix) system.

|- init.ndx
<br>
GROMACS index file that is used to define the [frozen] group, so that the particles in the group will be frozen in place during MD runs, forming walls that separates different periodic images.

|- minim_frozen.mdp
<br>
GROMACS mdp file for energy minimization.

|- nvt_frozen.mdp
<br>
GROMACS mdp file for NVT runs. The NVT run at 560K lasts for 1 ns. Atoms in the [frozen] group are frozen in place.

|- nvt_frozen_continue.mdp
<br>
GROMACS mdp file for NVT runs. This one is used when the last step found no suitable carboxyl-amine pair that could form amide bond. Unlike nvt_frozen.mdp, this run keeps the atom velocities in the last run.

|- example_run_script.py
<br>
Demo execution script. The default settings would run simulation crosslinking until degree of crosslinking reachs 30%.

|- DPC30.top
<br>
GROMACS topology file of the crosslinked system with a DPC of 30%.

|- DPC30.gro
<br>
GROMACS coordinate file of the crosslinked system with a DPC of 30%.

Please note that the files, "DPC30.top" and "DPC30.gro", are for your reference. You may get different results from the demo running script due to the randomness of each MD run.
<br><br>

Example_analysis: demo folder of a large system that was analyzed in the paper (In Press)

|- dry_mem.top
<br>
GROMACS topology file of the large crosslinked system with a DPC of 70%.

|- dry_mem.gro
<br>
GROMACS coordinate file of the large crosslinked system with a DPC of 70%.

|- solv_mem.top
<br>
GROMACS topology file of the large crosslinked system with a DPC of 70%, after solvation.

|- solv_mem.top
<br>
GROMACS coordinate file of the large crosslinked system with a DPC of 70%, after solvation.
