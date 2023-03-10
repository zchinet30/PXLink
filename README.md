# PXLink

PXLink (Polymer-Cross-Link) is a Python script that uses GROMACS (version 4.6 or higher) to automatically perfom simulation crosslinking and generate atomistic model of aromatic polyamide membrane separation layer. Currently it can use trimesic acid and m-Phenylenediamine as monomers.

# System requirements:

This script has been tested in Linux. To use it, you need to:
1. Install GROMACS 4.6 or newer. For details on installation, please refer to GROMACS website: https://www.gromacs.org/
2. Install Python 3.9 or newer.
3. Install required Python libraries: Numpy, Pandas, Networkx and Matplotlib

# How to use:

1. Place the following files in the same folder:
    - [] The PXLink script
    - [] The running script e.g. example_run_script.py
    - [] The force field you are using (We use a modified CHARMM36 force field. If you want to use other force fields you might need to modify the mark_atoms, add_bond and count_atoms methods in the PXLink script. Make sure this force field is included in the topology file.)
    - [] The initial system (both top and gro files, and optionally ndx file if you use it)
    - [] The mdp files needed for Gromacs runs
2. Modify the running script, designate the files used and assign desired values to the control variables.
3. Run the running script using commands like Python3 example_run_script.py
Note: it is also possible to restart a run from halfway by using top, gro and optionally ndx files generated in a halfway step instead of the initial system files. In that case, you likely need to modify the progress record variables to the values you left off (check the log file for them). 


# File list:

PXLink.py: PXLink script used to perform simulation crosslinking

Example_crosslinking: folder of a small system as an example of the crosslinking process

|- charmm36-mar2019.ff: folder containing a modified version of CHARMM36 force field used in simulating the crosslinking of aromatic polyamide polymer from monomers TMA (trimesic acid) and MPD (m-Phenylenediamine). Extra entries are added to merged.rtp (which describe the two monomer residues) and ffbonded.itp (extra angle and dihedral parameters needed for simulating the polymer).

|- init.top: Gromacs topology file of the initial (equilibrated monomer mix) system.

|- init.gro: Gromacs coordinate file of the initial (equilibrated monomer mix) system.

|- init.ndx: Gromacs index file that puts the monomer residues in the top and bottom of the periodic box in a [frozen] group, so that these residues can be frozen in place during MD runs, forming a wall that separates different periodic images in Z direction.

|- minim_frozen.mdp: Gromacs mdp file used in energy minimization

|- nvt_frozen.mdp: Gromacs mdp file used in MD runs during simulation crossllinking. The NVT run at 560K lasts for 1 ns. Atoms in the [frozen] group are frozen in place.

|- nvt_frozen_continue.mdp: Gromacs mdp file used in MD runs when the last step find no carboxyl-amine pair suitable for amide bond formation. Unlike nvt_frozen.mdp, this run keeps the atom velocities in the last run.

|- example_run_script.py: example running script that will take the initial system and run simulation crosslinking until degree of crosslinking reachs 30%.

|- DPC30.top: Gromacs topology file of the crosslinked system with a DPC of 30%. Note that due to the randomess of MD runs, you would not get a completely identical system from the running script.

|- DPC30.gro: Gromacs coordinate file of the crosslinked system with a DPC of 30%.

Example_analysis: folder of a large system used in the analysis paragraph of the manuscript

|- dry_mem.top:  Gromacs topology file of the large crosslinked system with a DPC of 70%.

|- dry_mem.gro:  Gromacs coordinate file of the large crosslinked system with a DPC of 70%.

|- solv_mem.top:  Gromacs topology file of the large crosslinked system with a DPC of 70%, after solvation.

|- solv_mem.top:  Gromacs coordinate file of the large crosslinked system with a DPC of 70%, after solvation.
