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

