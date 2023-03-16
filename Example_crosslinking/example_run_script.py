import PXLink as XL
import logging
import os

# ========the following section contains parameters to be assigned by user=========
# old_gmx: "True" if you are using GROMACS 4.x. "False" if you are using GROMACS 5 or newer.
old_gmx = True

# =====================run files==============================
# Name the input files required for the simulation.
# User can use files of systems output halfway during a crosslinking run to continue from there.
# In that case, do change the "Progress record variables" to where the crosslinking run stopped.
# gro: the coordinate file of the initial system
gro = 'init.gro'
# top: the topology file of the initial system. Should include the modified force field.
top = 'init.top'
# ndx: the index file containing the frozen groups. (Optional)
ndx = 'init.ndx'
# mdp_em: the mdp file for energy minimization
mdp_em = 'minim_frozen.mdp'
# mdp_NVT: the mdp file for NVT MD run
mdp_NVT = 'nvt_frozen.mdp'
# mdp_cont: the mdp file for NVT MD run, used when script fails to find C-N pair (Optional)
mdp_cont = 'nvt_frozen_continue.mdp'
# NOTE: 1. Parameters like simulation time and temperature are decided by these mdp files.
#       2. In mdp_NVT, the command "gen-vel" should be "yes" to generate velocities, while
#          in "mdp_cont" file, "gen-vel" should be "no".
#       3. If the user want to use the "constriant wall" method to create explicit membrane surfaces, mdp files should
#          use the "freezegrps" command to freeze the residues near Z-periodic boundaries.

# Name the dumping files:
# log: the detailed runtime log file
log = 'run_log.log'
# bondlog:the log file that records the length of each formed bond. (Optional)
bondlog = 'bonds.log'
# NOTE: The script will also dump a less detailed log file named (perfix) + "_script.log"

# ====================controlling parameters==============================
# Parameters that control the crosslinking run are defined in this part.
# The data structure of "dist", "zmin" and "zmax" is List.
# run_label: file name perfix of all files created during the run. Can include path.
run_label = 'Testrun'
# max_links: the simulation will end after this number of crosslinking bonds have been formed.
max_links = 1000
# max_dpc: the simulation will end after system reaches this degree of crosslinking.
max_dpc = 0.3
# dist: cutoff distance (nm) of carbon - nitrogen pairs allowed to form amide bonds.
dist = [0.3, 0.35, 0.4, 0.45, 0.5]
# set_shift: if the script fails to find any proper C-N pairs to form bond for this number of times,
#            move to the next dist[] value.
set_shift = 3
# NOTE: 1. At the beginning of the run, dist[0] is used to define the criteria of bond formation.
#       2. If the script failed to find any proper C-N pairs that satisfy the criteria for "set_shift"
#          times in a row, the next dist[], zmin[] and zmax[] values will be used until exhausted.

# =====================parameters about Z limits============================
# To prevent amide bonds from forming across periodic boundaries and connecting
# different periodic images, the script only allows the atoms in the "center region",
# i.e. zmin <= z <= zmax, to form amide bonds. Any atoms outside of this range will be ignored.
# This function will only be activated when "use_zlim" = "True".
# use_zlim: whether to activate the function that defines a restricted region for crosslinking.
use_zlim = False
# zmin: Atoms with Z coordinate below this value (nm) are ignored during crosslinking.
zmin = [0.3] * len(dist)
# zmax: Atoms with Z coordinate above this value (nm) are ignored during crosslinking.
zmax = [5.7] * len(dist)
# NOTE: 1. zmin and zmax are controlled by set_shift, just like dist.
#          That is to say, when the script fails to find C-N pairs set_shift times in a row,
#          it will use the next values in these lists.
#          For example, if dist = [0.3, 0.4], zmin = [1.0, 0.5], zmax = [3.0, 3.5] and set_shift = 3,
#          PXLink starts with a cutoff distance of 0.3 nm and tries to form amide bonds in 1.0 < z < 3.0 region;
#          when it fails to form bonds 3 times in a row, it will move on with a cutoff distance
#          of 0.4 nm, trying to form amide bonds in 0.5 < z < 3.5 region.
#       2. Another method to prevent amide bonds from forming across periodic boundaries is by
#          creating "constriant walls" of frozen residues near the Z-periodic boundaries,
#          which blocks movement and interaction of residues. This should be set up in the mdp files.
#          (Use freezegrps option in mdp files to freeze residues.)

# =====================parameters about NPT run=============================
# For a large system, the user may want to use either NPT or adjusting Z to keep system density stable.
# do_NPT: Whether to perform NPT runs to stablize system density.
do_NPT = False
# NPT_interval: an NPT run will be performed every this number of amide bonds have formed.
NPT_interval = 0
# NPT_after: NPT MD will not be performed until this number of amide bonds have formed.
NPT_after = 200
# mdp_NPT: the mdp file for NPT MD run.
mdp_NPT = ''
# NOTE: 1. The script will perform an NPT run to adjust the system density when every (NPT_after + NPT_interval * n)
#          amide bonds have formed, where n is any positive integer.
#       2. NPT_after is used with the "constriant wall" method, to prevent the frozen wall from moving
#          into center of the box in the NPT run. Recommend to set this value high enough
#          (e.g. 50%-60% of total carboxyl numbers). Also if the "constriant wall" method is used, the "constriant wall"
#          should be temporarily removed during the NPT run i.e., mdp_NPT should not use the "freezegrps" command.

# =====================parameters about adjusting Z=========================
# These parameters are used if you want to control the system density by adjusting the box length on z-axis (Z-length),
# using a "reduce box Z-length - energy minimization - NVT run" process automated by subroutine adjust_Z in PXLink.
# do_adjust_Z: whether to adjust the Z-length to stablize the system density.
do_adjust_Z = False
# adjust_Zmax: If this is True, when the script reduces Z-length, it will also reduce zmax.
adjust_Zmax = False
# adjust_interval: adjust the Z-length every this number of amide bonds have formed.
adjust_interval = 0
# adjust_after: the Z-length will not be adjusted until this number of amide bonds have formed.
adjust_after = 200
# mdp_adjZ_em: the mdp file for energy minimization after adjusting Z.
mdp_adjZ_em = ''
# mdp_adjZ_NVT: the mdp file for NVT run after adjusting Z.
mdp_adjZ_NVT = ''
# Thickness of the constrained wall, corresponding to the frozen groups.
z_freezelayer = 0.4
# NOTE: 1. adjust_interval and adjust_after works just like NPT_interval and NPT_after.
#       2.  Also if the "constriant wall" method is used, neither of the mdp files should
#           use the "freezegrps" command, just like mdp_NPT.

# =====================breakpoint variables==============================
# Breakpoint variables record the crosslinking process. They will be automatically adjusted
# during the simulation. But if the user wants to continue the simulation from a breakpoint,
# they need to change these variables manually to where the previous simulation terminated
# (which can be found in the scriptlog file).
# n_clink: the numbers of formed amide bonds; start from 0.
n_clink = 0
# loop: the number of polymerization loops; start from 1.
loop = 1
# p_set: the pointer position in dist[], zmin[] and zmax[]; start from 0.
p_set = 0
# no_pair_found: the number of times that the script failed to find any proper C-N pairs to form amide bond;
#                when "p_set" increases, this number will be reset to 0.
no_pair_found = 0
# continue_run: whether this run is a continuing run i.e., the script failed to find any proper C-N pairs
#               to form amide bond in the last polymerization loop. The initial value is "True".
continue_run = False
# dz: the total adjusted distance of Z-length. This is only needed when both "do_adjust_Z" and "adjust_Zmax"
#     are True.
dz = 0

# =======================Parameter assigning ends here=======================

# If a mdp_cont file is not provided, it will use the same mdp file as mdp_NVT.
if not mdp_cont:
    mdp_cont = mdp_NVT

# Stop GROMACS autobackup so that it won't break because of backing up mdout.mdp too many times.
os.environ['GMX_MAXBACKUP'] = '-1'
# Loads system information
Sys = XL.GromacsSys(topfile=top,
                    grofile=gro,
                    logfile=log,
                    ndxfile=ndx,
                    old_gmx=old_gmx)

if not zmin:
    zmin = [0] * len(dist)
if not zmax:
    zmax = [0] * len(dist)
if not len(dist) == len(zmin) == len(zmax):
    raise ValueError("Parameter sets should have same length!")

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
logger_script = logging.getLogger('L1')
handler1 = logging.FileHandler(run_label + '_scriptlog' + '.log')
handler1.setFormatter(formatter)
logger_script.setLevel('INFO')
logger_script.addHandler(handler1)

# Pre-crosslink equilibrium:
# logger_script.info(
#     "Start test run. A 2 ns NVT 300K equilibrium is run first."
# )
# Sys.gmx_run('md', 'nvt_preeq.mdp', run_label + '_initeq')

# Crosslinking begins.

# Sys.output_residue_network('test0.xml')

Sys.mark_atoms(nofreeze=True)
logger_script.info('Marked linking atoms with nofreeze=True\n')

if use_zlim:
    logger_script.info(
        f"Starting arguments:\ndist={dist[p_set]}\nzmin={zmin[p_set]}\nzmax={zmax[p_set]}"
    )
else:
    logger_script.info(f"Starting arguments:\ndist={dist[p_set]}\n")
if bondlog:
    with open(bondlog, 'a') as f:
        f.write('Index of atom A\tIndex of atom B\tA-B distance\t')
        f.write(
            'Whether bond goes across X=0\tWhether bond goes across Y=0\tDPC\n'
        )
        f.write('a  b   dist    across_X    across_Y   DPC\n')

z_range = [zmin[p_set], zmax[p_set]]
if do_adjust_Z:
    z_freespace = [z_freezelayer, Sys.box[2] - z_freezelayer]
    zl = z_freespace[1] - z_freespace[0]
    # n_water: the number of water molecules deleted since the last call of Sys.adjust_Z
    if n_clink > adjust_after + adjust_interval:
        n_water = (n_clink - adjust_after) % adjust_interval
    else:
        n_water = n_clink

while n_clink <= max_links and p_set < len(dist):
    logger_script.info(f'\nRunning loop {loop}, {n_clink} crosslink(s) done.')
    if not continue_run:
        # If a new bond formed during the last polymerization loop, an energy minimization will be performed,
        # and the NVT run restarts.
        # Run geometry optimization
        logger_script.info('Running geometry optimization.')
        Sys.gmx_run('opt', mdp_em, run_label + '_loop_' + str(loop) + '_opt')
        # Run NVT MD
        logger_script.info('Running NVT MD.')
        Sys.gmx_run('md', mdp_NVT, run_label + '_loop_' + str(loop) + '_nvt')
    else:
        # If no new bond formed during the last polymerization loop, then no energy minimization,
        # and NVT run should be a continuation of NVT run in the last polymerization loop.
        # Run NPT MD (continuation)
        logger_script.info('Continuing NVT MD.')
        Sys.gmx_run('md', mdp_cont, run_label + '_loop_' + str(loop) + '_nvt')
    # Search for carboxyl groups and amine groups that are close enough for crosslinking
    pair, pair_dist, across = Sys.CN_dist(r=dist[p_set],
                                          use_zlim=use_zlim,
                                          zlim=z_range)
    if pair:
        logger_script.info(f"Found C-N pair {pair}.")
        logger_script.info("Adding bond between the C-N pair.")
        continue_run = False
        # Add bond between the C-N pair found above
        Sys.add_bond(pair)
        n_clink += 1
        no_pair_found = 0
        Sys.output_contents(
            new_top=run_label + '_loop_' + str(loop) + '_topol.top',
            new_gro=run_label + '_loop_' + str(loop) + '_addbond.gro',
            change_files=True)
        if ndx:
            Sys.output_ndx(run_label + '_loop_' + str(loop) + '.ndx')
        logger_script.info("Bond added.")
        dpc = round(Sys.atom_count()["DPC_tri-linked-TMA_portion"], 3)
        if bondlog:
            with open(bondlog, 'a') as f:
                f.write(
                    f'{pair[0]}  {pair[1]}  {pair_dist}    {across[0]}  {across[1]}   {dpc}\n'
                )
        logger_script.info(f"There are now {Sys.count_clusters()} clusters.")
        logger_script.info(f"DPC is now {dpc}. Moving to next step.")
        if dpc > max_dpc:
            logger_script.info(
                f"Reached target DPC of {max_dpc}. Stop crosslinking.\n")
            break
        if n_clink == max_links:
            print("Max crosslink number reached. Break loop.\n")
            break
        # Perform NPT run if required
        if do_NPT and mdp_NPT and NPT_interval:
            if n_clink > NPT_after and n_clink % NPT_interval == 0:
                n_NPT = n_clink / NPT_interval
                logger_script.info(f"Running NPT run number {n_NPT}.")
                Sys.gmx_run('md', mdp_NPT, run_label + '_NPT_' + str(n_NPT))
        # Adjust Z-length if required
        if do_adjust_Z:
            n_water += 1
            if n_clink % adjust_interval == 0 and n_clink > adjust_after:
                z_freespace = [z_freezelayer, Sys.box[2] - z_freezelayer]
                zl = z_freespace[1] - z_freespace[0]
                adj_fn = run_label + '_loop_' + str(loop) + '_adjZ_'
                delta = Sys.adjust_z(zl=zl,
                                     zrange=z_freespace,
                                     nwater=n_water,
                                     filename=adj_fn,
                                     mdp_min=mdp_adjZ_em,
                                     mdp_nvt=mdp_adjZ_NVT)
                n_water = 0
                if adjust_Zmax:
                    dz = dz + delta
                    z_range = [zmin[p_set], zmax[p_set] - dz]
                    logger_script.info(
                        f'Reduced zmax by {dz} to {z_range[1]} now.')
    else:
        continue_run = True
        if bondlog:
            with open(bondlog, 'a') as f:
                f.write('\n')
        logger_script.info("No pair found! NVT MD will continue in next step.")
        no_pair_found += 1
    loop += 1

    # If the script failed to find any proper C-N pairs multiple times in a row, change the parameters.
    if no_pair_found >= set_shift:
        logger_script.info(
            f"Failed to find crosslink-able pairs in {set_shift} steps in a row. "
        )
        no_pair_found = 0
        p_set = p_set + 1
        dpc = round(Sys.atom_count()["DPC_tri-linked-TMA_portion"], 3)
        logger_script.info(f"Current degree of polymer crosslinking: {dpc}")
        if p_set >= len(dist):
            logger_script.info("No more crosslinks. Break loop.\n")
            break
        else:
            if use_zlim:
                z_range = [zmin[p_set], zmax[p_set]]
                if adjust_Zmax:
                    z_range = [zmin[p_set], zmax[p_set] - dz]
            logger_script.info(
                "Continue with adjusted crosslinking distance and Z range.")
            if use_zlim:
                logger_script.info(
                    f"dist={dist[p_set]} zmin={z_range[0]} zmax={z_range[1]}\n"
                )
            else:
                logger_script.info(f"dist={dist[p_set]}\n")

# Dump the modified coordinates and topology
logger_script.info("Outputting system after crosslink run.")
Sys.output_contents(new_top=run_label + '_topol_crosslink.top',
                    new_gro=run_label + '_conf_crosslink.gro')

# Remove the unreacted monomers. This includes both the monomers left unreacted and the frozen monomers.
# Comment out if you don't want to do it yet.
Unreacted_TMA = Sys.remove_residue('TMA')
Unreacted_MPD = Sys.remove_residue('MPD')
logger_script.info(f"Found {Unreacted_TMA} unreacted TMA monomers.")
logger_script.info(f"Found {Unreacted_MPD} unreacted MPD monomers.")
logger_script.info("Outputting system after removing unreacted monomers.")
Sys.output_contents(new_top=run_label + '_topol_crosslink_removed.top',
                    new_gro=run_label + '_conf_crosslink_removed.gro')

dpc = round(Sys.atom_count()["DPC_tri-linked-TMA_portion"], 3)
logger_script.info(f"Current degree of polymer crosslinking: {dpc}")
logger_script.info("Script run finished.")
Sys.output_residue_network(run_label + 'residue_network.xml')
