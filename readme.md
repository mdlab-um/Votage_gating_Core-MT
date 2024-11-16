# Voltage activation of BK channels
Data and Structure files related to voltage-gating simulations of Core-MT BK channels

# link to SI movies
https://drive.google.com/drive/folders/1TMUQh6gJTL85V_L7vUOJbX8oB733mGXv

*************************
Data formart

1  agr file
   can directly read by xmgrace , also human-readable and can be exported by groups (sear for S0/S1/S2 of target data set)

2  .dat/.xvg file
    human-readable file, can be plot by gnuplot or xmgrace

3.  .dx file
    grid file readable by vmd

*************************
*************************
1  inital (eqlibrated from cryo-EM close state ) and final structure (voltage-driven open structure) :

   in Inital_final_pdb folder, this are used to make  Fig 1 and 2E

*************************
2  Fig_2

Fig 2A : raw data : Fig_2A_gating_charge_move.agr  
         analysis script:  Fig_2A_gating_charge_move.py

Fig 2B :
         raw data :   Fig_2B_TM_move.agr 
         analysis script: Fig_2A_gating_charge_move.py  (same as 2A)  

Fig 2C:
         raw data :  Fig_2C_RMSD.agr  
         analysis script:  RMSD: use gromacs rms  fnction
         reference  structure:   Fig_2C_reference_open_structure.pdb

Fig 2D:
         raw data :  Fig_2D_helix_tilt.agr
         analysis script: Fig_2D_helix_tilt.py

Fig 2E :
         raw data, in upper folder "Inital_final_pdb"    

Fig 2F:
         rawdata :  Fig_2F_pore_profile.agr
         analysis script: Fig_2F_pore_profile_hole_input.inp   * this is input for progrome HOLE   

*************************
Fig 3:
         raw data: Fig_3_water_number.agr
                   Fig_3_conductance.agr
         analysis script:   Fig_3_water_number.py

*************************

Fig 4:
         raw data:  Fig_4A__zshift.agr
                    gnu type density data:  inactivated : Fig_4B_plot_inactivated_state_dashline.dat    
                                            activated   : Fig_4B_plot_activated_state_contour.dat 
 
         analysis script: to calculate the 2D movement:  Fig_4B_TM_1-5_2D.py
                          example gnu input to draw contour:  Fig_4B_plot_activated_state.gnu

*************************
Fig 5:
         raw data:   inactivated : Fig_5_potential_shifted_by_e-fild_inactive.dx
                     activated   : Fig_5_potential_shifted_by_e-fild_active.dx  

                                    *** these two file is too large , in case needed , please send email to coresponding author

         analysis script: origin vmd pmepot output does not have electic potential part,  need add the extra electric potential : Fig_5_add_external_potential_to_VMD_pmepot_output.py


*************************
*** special input files:
*************************

Pulling_simulation

   There is a steer MD stage                      : Protein_0-1000ns_1.0_steer.mdp
   and a constant reference position pulling stage: Protein_0-1000ns_1.0.mdp
   
   They are call  by                                0-1000ns.sh       

**** git

git init
git add  readme.txt
git commit -m "voltage gating of core MT"

 git add *
git commit -m "V5"


git remote add origin git@github.com:zhiguangjia/Votage_gating_Core-MT.git
git branch -M main
git push -u origin main

***  for later change:
#git add BC5_contact_probability.py
#git commit -m "Bc5 contact probability from a trajectory file"


#git add Fit_EC50.py
#git commit -m "Fit EC50 with fit_curve and double checked with numarical simulation"


username zhiguangjia
 token see laptoop code dataase

## Key contributors
Zhiguang Jia and Jianhan Chen

**Contact: Jianhan Chen** (jianhanc@umass.edu)


















