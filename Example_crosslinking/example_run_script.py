import XLink as XL
import logging
import os

old_gmx = True  # True if you are using Gromacs 4.x; False if higher versions.
# =====================Run files==============================
# Designate the files used for the run.
# You can use files of systems halfway during a crosslinking run, to continue from there.
# In that case, do change the "Progress record variables" to where you left off.
gro = 'init.gro'  # gro file of the initial system.
top = 'init.top'  # top file of the initial system
log = 'run_log.log'  # log file to be output
ndx = 'Testrun/run28_init.ndx'  # ndx file of the initial system (Optional)
bondlog = 'init.ndx'  # Output file that records the length of each bond when created. (Optional)
# ndx file is used when designating the frozen wall residues using an index file.
mdp_em = 'minim_frozen.mdp'  # mdp file for energy minimization
mdp_NVT = 'nvt_frozen.mdp'  # mdp file for NVT MD run
mdp_cont = 'nvt_frozen_continue.mdp'  # mdp file for continuation MD run (Optional)
# mdp_cont is used when fail to creat a new bond, and MD run continues from the last one.
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
run_label = 'Testrun'  # File name perfix of all files created during the run. Can include path.
max_links = 1000  # Crosslinking ends after this number of bonds are made.
max_dpc = 0.3  # Crosslinking ends after reaching this degree of crosslinking.
# dist: Max distance (nm) of carbon - nitrogen pair allowed for crosslinking.
dist = [0.3, 0.35, 0.4, 0.45, 0.5]
set_shift = 3  # when failing to find C-N pairs in a row for this many times, use next set of dist, zmin and zmax.

# =====================Parameters about Z limits============================
# In order to prevent amide bonds from forming across periodic boundary,
# the script can limit crosslinking in a given Z range zmin <= z <= zmax.
# Any atoms outside of this range are ignored when trying to find C-N pairs for crosslinking.
# These are given in lists with the same length as dist, and will be controlled by set_shift just like dist.
# This will only be used if use_zlim = True
use_zlim = True  # Whether to use zmin and zmax, to limit crosslinking in a given Z range.
# zmin: Atoms with Z coordinate below this value (nm) are ignored during crosslinking.
zmin = [0.3] * len(dist)
# Zmax: Atoms with Z coordinate above this value (nm) are ignored during crosslinking.
zmax = [5.7] * len(dist)

# =====================Parameters about NPT run=============================
# For a large system, you may want to use either NPT or adjusting Z to keep system density stable.
# This parameter is used if you want to control system density with NPT runs.
# If used, the script will do an NPT run each time a number of bonds are created.
# This will only be used if do_NPT = True
do_NPT = False
NPT_interval = 0  # If not 0, an NPT run will be performed every this many bonds are created.
NPT_after = 200  # NPT MD will not be done before this many bonds are created.
# Recommend to set this value high enough (e.g. 50%-60% of total carboxyl numbers), so as
# to prevent the frozen wall from moving into center of the box in the NPT run.
mdp_NPT = ''  # mdp file for NPT MD run

# =====================Parameters about adjusting Z=========================
# Thses parameters are used if you want to control system density by adjusting Z box length.
# If used, the script will reduce Z length of the period box each time a number of bonds are created,
# then do a short NVT equilibrium.
# They will only be used if do_adjust_Z = True
do_adjust_Z = False  # If True, will call the adjust_Z method after each bond creation to reduce box Z to keep density.
adjust_Zmax = False  # If this is True, will reduce zmax altogether with Z length of peridic box.
adjust_interval = 0  # Each time this many bonds are created, Z length will be adjusted and a short NVT run is done.
adjust_after = 200  # Z length will not be adjusted before this many bonds are created.
# Recommend to set this value high enough (e.g. 50%-60% of total carboxyl numbers), so as
# to prevent the frozen wall from moving into center of the box in the following NVT run.
mdp_adjZ_em = ''  # mdp file for energy minimization after adjusting Z.
mdp_adjZ_NVT = ''  # mdp file for NVT after adjusting Z.
# Note that this shall not freeze any of the residues.
z_freezelayer = 0.4  # Thickness of frozen walls. (Top and bottom walls are assumed to have same thickness)
# For the purpose of calculating density, only atoms within this range will be used

# =====================Progress record variables==============================
# Variables that records how many crosslinking steps have been performed.
# Need to change to where you left off when continuing a run.
# Also note that if you continue from a step where NPT or adjusting Z should be run, it might be skipped.
n_clink = 0  # numbers of amide bonds formed.
loop = 1  # number of polymerization steps run. (1 + steps run)
p_set = 0  # script is using dist[p_set], zmin[p_set], zmax[p_set] now.
no_pair_found = 0  # steps of failing to find linkable C-N pairs in a row.
continue_run = False  # whether the MD simulation next step will be a continuation of the last.
# True when no crosslinking bond is formed in last step.
dz = 0  # Accumulative delta Z. Only used if do_adjust_Z and adjust_Zmax are both True.
# This variable is used to count how much Z length has reduced due to adjust_Z, thus reducing zmax.

if not zmin:
    zmin = [0] * len(dist)
if not zmax:
    zmax = [0] * len(dist)
if not len(dist) == len(zmin) == len(zmax):
    raise ValueError("Parameter sets should have same length!")

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
logger_script = logging.getLogger('L1')
handler1 = logging.FileHandler(run_label + '_scriptlog_' + '.log')
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
    logger_script.info(
        f"Starting arguments:\ndist={dist[p_set]}\n"
    )
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
        Sys.gmx_run('opt', mdp_em, run_label + '_step_' + str(loop) + '_opt')
        # Run NPT MD
        logger_script.info('Running NVT MD.')
        Sys.gmx_run('md', mdp_NVT, run_label + '_step_' + str(loop) + '_nvt')
    # If not, then no optimization, and NVT run should be a continuation of last run.
    else:
        # Run NPT MD (continuation)
        logger_script.info('Continuing NVT MD.')
        Sys.gmx_run('md', mdp_cont, run_label + '_step_' + str(loop) + '_nvt')
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
            new_top=run_label + '_step_' + str(loop) + '_topol.top',
            new_gro=run_label + '_step_' + str(loop) + '_addbond.gro',
            change_files=True)
        if ndx:
            Sys.output_ndx(run_label + '_step_' + str(loop) + '.ndx')
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
                adj_fn = run_label + '_step_' + str(loop) + '_adjZ_'
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
                    logger_script.info(f'Reduced zmax by {dz} to {z_range[1]} now.')
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
                    f"dist={dist[p_set]} zmin={z_range[0]} zmax={z_range[1]}\n")
            else:
                logger_script.info(
                    f"dist={dist[p_set]}\n")

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
