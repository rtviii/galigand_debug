How i generated the ligand. In pymol:
```
load 6XKI.pdb
sele resn V6D
save V6D.pdb, sele
```

Then parameter file generaton. Used OpenBabel to convert to MOL2:

`obabel -ipdb V6D.pdb -omol2 -O V6D.mol2`

Generated `.params` file:

`$ROSETTA/source/scripts/python/public/molfile_to_params.py  -n V6D -p V6D --clobber V6D.mol2`

To find out exactly which residue V6D ligand is in holo pdb file (`6XKI.pdb`):
```xml
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="ref2015" weights="ref2015"/>
  </SCOREFXNS>
  <RESIDUE_SELECTORS>
    <ResidueName name="v6d_ligand" residue_names="V6D"/>
  </RESIDUE_SELECTORS>
  <SIMPLE_METRICS>
    <SelectedResiduesMetric 
      name="v6d_residues" 
      residue_selector="v6d_ligand"
      rosetta_numbering="false"
    />
  </SIMPLE_METRICS>
  <MOVERS>
    <RunSimpleMetrics name="collect_metrics" metrics="v6d_residues"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="collect_metrics"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
```
Which, when ran on 6XKI.pdb (unaltered) produced:

```
SEQUENCE: 
SCORE: total_score dslf_fa13    fa_atr    fa_dun   fa_elec fa_intra_rep fa_intra_sol_xover4              fa_rep              fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb lk_ball_wtd       omega     p_aa_pp pro_close rama_prepro         ref yhh_planarity selection description 
SCORE:    5873.306     0.000 -2524.353   839.894  -429.463       23.205             116.262            6216.271            1734.467     -36.827     -56.083     -24.276     -96.107     -64.341      32.797     -41.734    11.212      34.459     137.912         0.011      503A 6XKI_0001
```



# Minimal galigand workflow:

```xml
<ROSETTASCRIPTS>

    <SCOREFXNS>
        <ScoreFunction name="ga_dock_score" weights="beta_genpot">
            <Reweight scoretype="fa_sol" weight="1.1"/>
            <Reweight scoretype="fa_rep" weight="0.5"/>
            <Reweight scoretype="gen_bonded" weight="1.0"/>
        </ScoreFunction>

        <ScoreFunction name="output_score" weights="beta_genpot">
             <Reweight scoretype="fa_elec" weight="1.2"/>
             <Reweight scoretype="fa_sol" weight="1.1"/>
             <Reweight scoretype="fa_rep" weight="1.0"/>
             <Reweight scoretype="gen_bonded" weight="1.0"/>
         </ScoreFunction>
    </SCOREFXNS>

    <!-- <RESIDUE_SELECTORS>
        <Chain name="rna_chain" chains="B"/>
        <ResidueName name="rna_selector" residue_name3="RGU,RCY,RAD,RU,RC,RA,RG" />
        <Chain name="ligand_chain" chains="L" />
    </RESIDUE_SELECTORS>

    <MOVE_MAP_FACTORIES>
        <MoveMapFactory name="freeze_rna_factory" bb="false" chi="false">
            <Backbone residue_selector="rna_chain" enable="false"/>
            <Chi residue_selector="rna_chain" enable="false"/>
        </MoveMapFactory>
    </MOVE_MAP_FACTORIES> -->

    <SCORINGGRIDS ligand_chain="L" width="45">
        <ClassicGrid grid_name="classic" weight="1.0"/>
    </SCORINGGRIDS>

    <MOVERS>
        <GALigandDock name="ga_dock_minimal"
                      runmode="dockrigid"         
                      scorefxn="ga_dock_score"    
                      grid_step="1"
                      padding="5.0">

        </GALigandDock>
    </MOVERS>

    <PROTOCOLS>
        <Add mover_name="ga_dock_minimal"/>
    </PROTOCOLS>

    <OUTPUT scorefxn="output_score"/>

</ROSETTASCRIPTS>
```
I run it via:
```sh

rosetta_scripts.default.macosclangrelease \
  -database $ROSETTA/database/ \
  -s 6XKI.pdb \
  -parser:protocol galigand_minimal.xml \
  -nstruct 10 \
  -gen_potential \
  -ex1 -ex2 \
  -use_input_sc \
  -extra_res_fa V6D.params \
  -out:suffix _galigand \
  -out:file:scorefile score_galig.sc \
  -ignore_unrecognized_res
```


This yields Rama torsion errors:

```
protocols.jd2.JobDistributor: [ WARNING ] 6XKI_galigand_0001 reported failure and will NOT retry
protocols.jd2.PDBJobInputter: PDBJobInputter::pose_from_job
protocols.jd2.PDBJobInputter: filling pose from saved copy 6XKI.pdb
protocols.jd2.PDBJobInputter: PDBJobInputter::pose_from_job
protocols.jd2.PDBJobInputter: filling pose from saved copy 6XKI.pdb
protocols.rosetta_scripts.ParsedProtocol: =======================BEGIN MOVER GALigandDock - ga_dock_minimal=======================
core.scoring.electron_density.ElectronDensity: Loading Density Map
core.scoring.electron_density.ElectronDensity: [ WARNING ] No density map specified
protocols.ligand_docking.GALigandDock.GALigandDock: Preparing grid using input pose
protocols.ligand_docking.GALigandDock.GridScorer:   ... built grid with 20.419 + 5 A radius; 1 A grid spacing
protocols.ligand_docking.GALigandDock.GridScorer: Initial dimension: 51 x 51 x 51
protocols.ligand_docking.GALigandDock.GridScorer: Initial origin = ( -52.1365,-62.3555,-61.4605 )
protocols.ligand_docking.GALigandDock.GridScorer: Received 10 residues for grid calculation:  RAD:LowerRNA:Virtual_Phosphate:5prime_end_OH RGU RAD RGU RAD RGU RAD RGU RAD RGU:UpperRNA:3prime_end_OH
protocols.ligand_docking.GALigandDock.GALigandDock: Build grid for 10 ligand residues and for 0 sidechain residues
protocols.ligand_docking.GALigandDock.GALigandDock: Mobile sidechains: 382+383+384+385+386+387+388+389+390+391
protocols.ligand_docking.GALigandDock.GridScorer: Building 51 x 51 x 51 grid; origin = ( -52.1365,-62.3555,-61.4605 )
protocols.ligand_docking.GALigandDock.GridScorer:    maxdis = 7.51668
protocols.ligand_docking.GALigandDock.GridScorer:    subhash_buffer = 1.73205
protocols.ligand_docking.GALigandDock.GridScorer:    hash_grid = 8
protocols.ligand_docking.GALigandDock.GridScorer: Calculating gridded energies for 17 unique atom types
protocols.ligand_docking.GALigandDock.GridScorer:   and lkb virtual parameters for 8 unique polar atom types
protocols.ligand_docking.GALigandDock.GridScorer: Grid calculation took 2.79136 seconds.
protocols.ligand_docking.GALigandDock.GALigandDock: Setting up Ligand aligner.
protocols.ligand_docking.GALigandDock.GridScorer: Recalculating grids with smoothing of 0.75
protocols.ligand_docking.GALigandDock.GridScorer:   ... built grid with 20.419 + 5 A radius; 2 A grid spacing
protocols.ligand_docking.GALigandDock.GridScorer: Initial dimension: 26 x 26 x 26
protocols.ligand_docking.GALigandDock.GridScorer: Initial origin = ( -52.1365,-62.3555,-61.4605 )
protocols.ligand_docking.GALigandDock.GridScorer: Received 10 residues for grid calculation:  RAD:LowerRNA:Virtual_Phosphate:5prime_end_OH RGU RAD RGU RAD RGU RAD RGU RAD RGU:UpperRNA:3prime_end_OH
protocols.ligand_docking.GALigandDock.GridScorer: Building 26 x 26 x 26 grid; origin = ( -52.1365,-62.3555,-61.4605 )
protocols.ligand_docking.GALigandDock.GridScorer:    maxdis = 7.51668
protocols.ligand_docking.GALigandDock.GridScorer:    subhash_buffer = 3.4641
protocols.ligand_docking.GALigandDock.GridScorer:    hash_grid = 8
protocols.ligand_docking.GALigandDock.GridScorer: Calculating gridded energies for 17 unique atom types
protocols.ligand_docking.GALigandDock.GridScorer:   and lkb virtual parameters for 8 unique polar atom types
protocols.ligand_docking.GALigandDock.GridScorer: Grid calculation took 0.538043 seconds.
core.energy_methods.GenericBondedEnergy: Creating new peptide-bonded energy container (381)
ligand_align: Defined 110 receptor Vsites.
ligand_align: Made total 37 phores from 79 clusters.
ligand_align: Too many sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC 18 -> 19
ligand_align: Defined 108 receptor Vsites.
ligand_align: Made total 37 phores from 79 clusters.
ligand_align: Too many sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC 18 -> 20
ligand_align: Defined 99 receptor Vsites.
ligand_align: Made total 35 phores from 73 clusters.
ligand_align: Too many sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC 18 -> 21
ligand_align: Defined 95 receptor Vsites.
ligand_align: Made total 31 phores from 69 clusters.
ligand_align: Too many sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC 18 -> 22
ligand_align: Defined 86 receptor Vsites.
ligand_align: Made total 30 phores from 65 clusters.
ligand_align: Too many sites, redefining with receptor sites with NNEIGH_CUT_EXPOSED_BB/SC 18 -> 23
ligand_align: Defined 78 receptor Vsites.
ligand_align: Made total 26 phores from 56 clusters.
protocols.ligand_docking.GALigandDock.GridScorer: Recalculating grids with smoothing of 0.375
protocols.ligand_docking.GALigandDock.GALigandDock: Construct a separate grid for LigandAligner, with 0.5 grid step and no movable scs.
ligand_align: Total 80 ligand phores defined:
ligand_align: Ligand phore 1: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 2: L.O5'/L.HO5'/L.OP2
ligand_align: Ligand phore 3: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 4: L.O3'/L.O2'/L.O4'
ligand_align: Ligand phore 5: L.O3'/L.HO2'/L.O5'
ligand_align: Ligand phore 6: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 7: L.N1/L.N3/L.N7
ligand_align: Ligand phore 8: L.N1/L.H61/L.H62
ligand_align: Ligand phore 9: L.H61/L.O6/L.H1
ligand_align: Ligand phore 10: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 11: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 12: L.O3'/L.O2'/L.O5'
ligand_align: Ligand phore 13: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 14: L.O2'/L.N3/L.O4'
ligand_align: Ligand phore 15: L.N3/L.H1/L.H22
ligand_align: Ligand phore 16: L.O6/L.N1/L.H61
ligand_align: Ligand phore 17: L.HO2'/L.O5'/L.O4'
ligand_align: Ligand phore 18: L.H1/L.H21/L.N3
ligand_align: Ligand phore 19: L.O5'/L.O3'/L.OP2
ligand_align: Ligand phore 20: L.O4'/L.O3'/L.O2'
ligand_align: Ligand phore 21: L.O3'/L.HO2'/L.O5'
ligand_align: Ligand phore 22: L.O2'/L.HO2'/L.O4'
ligand_align: Ligand phore 23: L.N1/L.N3/L.N7
ligand_align: Ligand phore 24: L.N7/L.H61/L.H62
ligand_align: Ligand phore 25: L.H61/L.O6/L.H1
ligand_align: Ligand phore 26: L.OP2/L.OP1/L.O5'
ligand_align: Ligand phore 27: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 28: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 29: L.O3'/L.HO2'/L.O5'
ligand_align: Ligand phore 30: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 31: L.N3/L.H1/L.H22
ligand_align: Ligand phore 32: L.O6/L.N7/L.H62
ligand_align: Ligand phore 33: L.O6/L.N1/L.N3
ligand_align: Ligand phore 34: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 35: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 36: L.O3'/L.O2'/L.O5'
ligand_align: Ligand phore 37: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 38: L.N1/L.N7/L.H61
ligand_align: Ligand phore 39: L.N7/L.H62/L.N7
ligand_align: Ligand phore 40: L.HO2'/L.O5'/L.O4'
ligand_align: Ligand phore 41: L.H61/L.H62/L.O6
ligand_align: Ligand phore 42: L.O4'/L.O3'/L.O2'
ligand_align: Ligand phore 43: L.O3'/L.HO2'/L.O5'
ligand_align: Ligand phore 44: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 45: L.O2'/L.HO2'/L.O4'
ligand_align: Ligand phore 46: L.N3/L.H1/L.H22
ligand_align: Ligand phore 47: L.O6/L.N7/L.H62
ligand_align: Ligand phore 48: L.O6/L.H1/L.N1
ligand_align: Ligand phore 49: L.H1/L.H21/L.N3
ligand_align: Ligand phore 50: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 51: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 52: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 53: L.N1/L.N3/L.N7
ligand_align: Ligand phore 54: L.N1/L.H61/L.H62
ligand_align: Ligand phore 55: L.O4'/L.O3'/L.O2'
ligand_align: Ligand phore 56: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 57: L.O2'/L.HO2'/L.O4'
ligand_align: Ligand phore 58: L.N3/L.H1/L.H22
ligand_align: Ligand phore 59: L.O6/L.H1/L.N1
ligand_align: Ligand phore 60: L.O6/L.H61/L.H62
ligand_align: Ligand phore 61: L.H1/L.H21/L.N3
ligand_align: Ligand phore 62: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 63: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 64: L.O3'/L.HO2'/L.O5'
ligand_align: Ligand phore 65: L.O3'/L.OP2/L.OP1
ligand_align: Ligand phore 66: L.N1/L.N3/L.N7
ligand_align: Ligand phore 67: L.N1/L.H61/L.N3
ligand_align: Ligand phore 68: L.N7/L.H62/L.O6
ligand_align: Ligand phore 69: L.H61/L.H1/L.H22
ligand_align: Ligand phore 70: L.O5'/L.O4'/L.O3'
ligand_align: Ligand phore 71: L.O4'/L.O2'/L.HO2'
ligand_align: Ligand phore 72: L.O3'/L.O2'/L.HO3'
ligand_align: Ligand phore 73: L.N3/L.H1/L.H21
ligand_align: Ligand phore 74: L.N3/L.N7
ligand_align: Ligand phore 75: L.N3/L.O4'
ligand_align: Ligand phore 76: L.N3/L.N7
ligand_align: Ligand phore 77: L.N3/L.N7
ligand_align: Ligand phore 78: L.HO3'/L.HO2'
ligand_align: Ligand phore 79: L.H21
ligand_align: Ligand phore 80: L.O5'
ligand_align: Found 101 ligand-receptor phore matches.
ligand_align: Large num. of n_phore_matches (101); Nautogen set by 70
protocols.ligand_docking.GALigandDock.GALigandDock: Automatically set nstruct-from-reference to 35 (from 70 trials) of total 100 left to sample.
protocols.ligand_docking.GALigandDock.GALigandDock: Est. time for matching: 70~140 seconds...

ERROR: Error in core::scoring::RamaPrePro::random_mainchain_torsions(): No mainchain score table for residue type RGU exists.
ERROR:: Exit from: src/core/scoring/RamaPrePro.cc line: 280
protocols.rosetta_scripts.ParsedProtocol: [ ERROR ] Exception while processing protocol:

File: src/core/scoring/RamaPrePro.cc:280
[ ERROR ] UtilityExitException
ERROR: Error in core::scoring::RamaPrePro::random_mainchain_torsions(): No mainchain score table for residue type RGU exists.


protocols.jd2.JobDistributor: [ ERROR ]

[ERROR] Exception caught by JobDistributor for job 6XKI_galigand_0002

[ ERROR ]: Caught exception:


File: src/core/scoring/RamaPrePro.cc:280
[ ERROR ] UtilityExitException
ERROR: Error in core::scoring::RamaPrePro::random_mainchain_torsions(): No mainchain score table for residue type RGU exists.


AN INTERNAL ERROR HAS OCCURED. PLEASE SEE THE CONTENTS OF ROSETTA_CRASH.log FOR DETAILS.

```
