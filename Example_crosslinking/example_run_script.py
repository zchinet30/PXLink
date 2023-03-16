import XLink as XL
import logging
import os

# ========The following section contains parameters to be assigned by user=========
# old_gmx: Whether Gromacs 4.x or higher versions are used.
old_gmx = True

# =====================Run files==============================
# Designate the input files used for the run.
# User can use files of systems output halfway during a crosslinking run to continue from there.
# In that case, do change the "Progress record variables" to where the crosslinking run stopped.
# gro: gro file of the initial system
gro = 'init.gro'
# top: top file of the initial system. Should include the modified force field.
top = 'init.top'
# ndx: ndx file of the initial system, used to designate frozen wall residues (Optional)
ndx = 'init.ndx'
# mdp_em: mdp file for energy minimization
mdp_em = 'minim_frozen.mdp'
# mdp_NVT: mdp file for NVT MD run
mdp_NVT = 'nvt_frozen.mdp'
# mdp_cont: mdp file for NVT MD run, continuing from failing to find C-N pair (Optional)
mdp_cont = 'nvt_frozen_continue.mdp'
# NOTE: Parameters like simulation time and temperature are decided by these mdp files.
# Also in mdp_NVT, gen-vel should be yes, since these runs are started after PXLink modifying
# the system to add an amide bond; in mdp_cont, it can be no. If the user want to use
# the "constriant wall" method to create explicit membrane surfaces, these mdp files should
# use the freezegrps option to freeze the residues near Z-periodic boundaries.

# Output files:
# log: detailed log file to be output
log = 'run_log.log'
# bondlog: Output file that records the length of each bond when created. (Optional)
bondlog = 'bonds.log'
# NOTE: The script will also output a less detailed log file named (perfix) + "_script.log"

# =====================Run control parameters==============================
# Parameters that control the crosslinking run are defined in this part.
# dist, zmin and zmax are given in lists. At the beginning of the run, dist[0], zmin[0]
# and zmax[0] are used to set the conditions of amide bond forming. After the script fails
# to find carbon - nitrogen pairs that match the conditions for set_shift times in a row,
# the script will go on with dist[1], zmin[1] and zmax[1], and so on. For example, with the default
# settings, the script will form amide bonds between a pair of carbon and nitrogen if their distance is
# less then 0.3 nm and both atoms' Z coordinates are between 0.3 - 5.0 nm. After failing to find such
# pairs 3 times in a row, the script will move on to find pairs with distance less then 0.35 nm and
# both atoms' Z coordinates are between 0.3 - 5.0 nm.
# run_label: file name perfix of all files created during the run. Can include path.
run_label = 'Testrun'
# max_links: Crosslinking ends after this number of bonds are made.
max_links = 1000
# max_dpc: Crosslinking ends after reaching this degree of crosslinking.
max_dpc = 0.3
# dist: Max distance (nm) of carbon - nitrogen pair allowed for crosslinking.
dist = [0.3, 0.35, 0.4, 0.45, 0.5]
# set_shift: when failing to find C-N pairs in a row for this many times, use next dist value.
set_shift = 3

# =====================Parameters about Z limits============================
# In order to prevent amide bonds from forming across periodic boundary and produce membrane models
# with explicit membrane surfaces, the script can limit crosslinking in a "center region" zmin <= z <= zmax.
# Any atoms outside of this range are ignored when trying to find C-N pairs for crosslinking.
# This will only be used if use_zlim = True
# use_zlim: Whether to use zmin and zmax, to limit crosslinking in a given Z range.
use_zlim = False
# zmin: Atoms with Z coordinate below this value (nm) are ignored during crosslinking.
zmin = [0.3] * len(dist)
# zmax: Atoms with Z coordinate above this value (nm) are ignored during crosslinking.
zmax = [5.7] * len(dist)
# NOTE: zmin and zmax are controlled by set_shift, just like dist.
# That is to say, dist, zmin and zmax are lists, and when the script fails to find
# C-N pairs set_shift times in a row, it will use the next values in these lists.
# For example, if dist = [0.3, 0.4], zmin = [1.0, 0.5], zmax = [3.0, 3.5] and set_shift = 3,
# PXLink starts with a cutoff distance of 0.3 nm and tries to form amide bonds in 1.0 < z < 3.0 region;
# when it fails to form bonds 3 times in a row, it will move on with a cutoff distance
# of 0.4 nm, trying to form amide bonds in 0.5 < z < 3.5 region.
# NOTE: Another method to keep membrane surfaces is by creating "constriant walls" of frozen residues
# near the Z-periodic boundaries, which blocks movement and interaction of residues. This should be
# set up in the mdp files. (Use freezegrps option in mdp files to freeze residues.)

# =====================Parameters about NPT run=============================
# For a large system, you may want to use either NPT or adjusting Z to keep system density stable.
# This parameter is used if you want to control system density with NPT runs.
# If used, the script will do NPT runs regularly.
# This will only be used if do_NPT = True
# do_NPT: Whether to use NPT run to stablize system density.
do_NPT = False
# NPT_interval: an NPT run will be performed every NPT_interval amide bonds are created.
NPT_interval = 0
# NPT_after: NPT MD will not be run before NPT_after amide bonds are created.
NPT_after = 200
# mdp_NPT: mdp file for NPT MD run.
mdp_NPT = ''
# NOTE: It means that when the script forms (NPT_after + NPT_interval * n) amide bonds, where
# n is any positive integer, it will perform an NPT run to adjust system density.
# NPT_after is used with the "constriant wall" method, to prevent the frozen wall from moving
# into center of the box in the NPT run. Recommend to set this value high enough
# (e.g. 50%-60% of total carboxyl numbers). Also if the "constriant wall" method is used altogether,
# the "constriant wall" should be temporarily removed during the NPT run i.e., mdp_NPT should not use freezegrps.

# =====================Parameters about adjusting Z=========================
# Thses parameters are used if you want to control system density by adjusting Z box length,
# using a "reduce simulation cell Z-length - energy minimization - NVT run" process automated by subroutine
# adjust_Z in PXLink. They will only be used if do_adjust_Z = True
# do_adjust_Z: Whether to use NPT run to stablize system density.
do_adjust_Z = False
# adjust_Zmax: Whether to reduce zmax altogether with Z-length of simulation cell (if you also use "center region").
adjust_Zmax = False
# adjust_interval: Each time adjust_interval amide bonds are created, Z length will be adjusted.
adjust_interval = 0
# adjust_after: Z length will not be adjusted before adjust_after bonds are created.
adjust_after = 200
# mdp_adjZ_em: mdp file for energy minimization after adjusting Z.
mdp_adjZ_em = ''
# mdp_adjZ_NVT: mdp file for NVT run after adjusting Z.
mdp_adjZ_NVT = ''
# Thickness of "constriant wall". This is excluded when calculating density.
z_freezelayer = 0.4
# For the purpose of calculating density, only atoms within this range will be used
# NOTE: adjust_interval and adjust_after works just like NPT_interval and NPT_after
# in "Parameters about NPT run" section. Also if the "constriant wall" method is used,
# neither of the mdp files should use freezegrps option, just like mdp_NPT.

# =====================Progress record variables==============================
# Variables that records crosslinking progression. These variables automatically changes
# along crosslinking process, but if the user wants to continue from a half-way crosslinking run,
# they need to change these variables to where they left off (which can be found in the scriptlog file).
# n_clink: numbers of amide bonds formed. Starts from 0.
n_clink = 0
# loop: number of polymerization loops run, no matter whether a bond is formed. Starts from 1.
loop = 1
# p_set: which value in dist, zmin and zmax the script is using now. i.e., how many times have the
# script failed to find C-N pairs to form amide bond set_shift times in a row.
p_set = 0
# no_pair_found: tracks how many times the script has failed to find C-N pairs to form amide bond in a row.
# When no_pair_found equals set_shift, the script moves on with next cutoff distance value, and no_pair_found
# is reset to 0. Starts from 0.
no_pair_found = 0
# continue_run: whether the script fail to form an amide bond in the last loop, so that the next MD run should
# use mdp_cont instead of mdp_NVT. Starts from True.
continue_run = False
# dz: how much Z-periodic box length has been reduced due to adjust_Z. This is only needed when both
# do_adjust_Z and adjust_Zmax are used.
dz = 0

# =======================Parameter assigning ends here=======================

# mdp_cont is used when failing to creat a new bond, and MD run continues from the last one.
if not mdp_cont:
    mdp_cont = mdp_NVT

# No autobackup so that script doesn't end because of too many auto backups.
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
    # n_water: number of water liberated since last calling of Sys.adjust_Z
    if n_clink > adjust_after + adjust_interval:
        n_water = (n_clink - adjust_after) % adjust_interval
    else:
        n_water = n_clink

while n_clink <= max_links and p_set < len(dist):
    logger_script.info(f'\nRunning loop {loop}, {n_clink} crosslink(s) done.')
    # If a new crosslinking is made last loop, geometry optimization is needed, and NVT run restarts.
    if not continue_run:
        # Run geometry optimization (may not be good in actual run)
        logger_script.info('Running geometry optimization.')
        Sys.gmx_run('opt', mdp_em, run_label + '_loop_' + str(loop) + '_opt')
        # Run NPT MD
        logger_script.info('Running NVT MD.')
        Sys.gmx_run('md', mdp_NVT, run_label + '_loop_' + str(loop) + '_nvt')
    # If not, then no optimization, and NVT run should be a continuation of last run.
    else:
        # Run NPT MD (continuation)
        logger_script.info('Continuing NVT MD.')
        Sys.gmx_run('md', mdp_cont, run_label + '_loop_' + str(loop) + '_nvt')
    # Find carboxyl group and amine group close enough for crosslinking
    pair, pair_dist, across = Sys.CN_dist(r=dist[p_set],
                                          use_zlim=use_zlim,
                                          zlim=z_range)
    if pair:
        logger_script.info(f"Found C-N pair {pair}.")
        logger_script.info("Adding bond between the C-N pair.")
        continue_run = False
        Sys.add_bond(pair)  # Add bond between the C-N pair found above
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
        # Adjust Z length if required
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
        no_pair_found += 1  # Adjust arguments if fail to find pairs multiple times in a row
    loop += 1

    # If script fails to find new C-N pairs multiple times in a row, change parameters.
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

# Output modified coordinates and topology
logger_script.info("Outputting system after crosslink run.")
Sys.output_contents(new_top=run_label + '_topol_crosslink.top',
                    new_gro=run_label + '_conf_crosslink.gro')

# Remove unreacted monomers. This includes both monomers left unreacted and frozen monomers.
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
