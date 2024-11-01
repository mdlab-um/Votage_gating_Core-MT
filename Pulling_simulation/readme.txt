#!/bin/bash


##########   model  

in  /home/zgjia/Project/BK_channel/Simulation_hBK/Anton/Core_MT_modeled_750mv_3  we have 750 mV driven open

   cp -r /home/zgjia/Project/BK_channel/Simulation_hBK/Anton/Core_MT_modeled_750mv_3/analysis .

   it start from   ../Core_MT_modeled/workdir.2_Eqlibrate.dms

                         <==  in attempt 5  it use     cp /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/0V_cter_restrain_symmetric/100ns.gro initial.gro

                                 <==  from   /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/Eqlibration_form_dewetted/step6.7_lipid_relaxtion.gro

                                           <==  from  /home/zgjia/Project/BK_channel/Simulation_test_with_EM_result/Seperate_S0_amber/close_PC_minitail_T/step6.6_equilibration.pdb

                                 ==> also used in /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/from_syn_restrain_0mv

 
1)

same  here

 #  cp /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/0V_cter_restrain_symmetric/100ns.gro initial.gro
 
 #  ref  /home/zgjia/Project/BK_channel/Simulation_hBK/Anton/Core_MT_modeled_750mv_3/gating_charge_ref_filter_S1-5.agr

    cp /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/0V_cter_restrain_symmetric/100ns.gro initial.gro



################################################################################


################################################################################

  res train in /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/0V_cter_restrain_symmetric/./Gen_restrain_less.py

   ./Gen_restrain_less.py   toppar/PROAB.itp  toppar/PROAB_rest.itp

    272-298 200

    >341 800

   in   use 0.1 kcal/mol  A   = 48 ?   

#   OK  use  50 here

  use 500 1 kcal 
  

################   pull

  ref  /home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_from_syn/restrain_pull_250mV_tripull_0.6kcal/commond


   ***************************

    idea is a test pullimng under electic environment

         **8 test in  /home/zgjia/test/test_gmx_pull_two_group



   ***************************  inital distance ??

   ./Gating_charge_move.py Eqilibrated_system_with_chain.pdb initial.gro

gate center  [76.40583483 75.4583327  47.29083443]
210
CZ  1.4791622161865234 A
CZ  1.1891651153564453 B
CZ  3.929166793823242 C
CZ  -0.5108356475830078 D
213
CZ  9.25916862487793 A
CZ  8.919164657592773 B
CZ  10.65916633605957 C
CZ  11.199163436889648 D

   we will target open state Anton

   ./Gating_charge_move.py Eqilibrated_system_with_chain.pdb ../from_syn_restrain_N750mv_anton/from_anton_open.gro


210                             increase compare above
CZ  -3.1558361053466797 A        4.6                              
CZ  -9.74583625793457 B         10.8
CZ  -7.875833511352539 C        11.8
CZ  -6.745832443237305 D        6.0       ~7.5    in paper 23 have ~7   average distance  -7.0
213
CZ  3.044168472290039 A         12
CZ  -3.3458385467529297 B       12 
CZ  -3.015836715698242 C        13  
CZ  -1.2158336639404297 D       10         11.7  in paper ~ 11          average distance -1.13  

check:
gmx_mpi grompp -f Protein_0-1000ns_1.0.mdp -o  Protein_0-1000ns_2020.tpr -c ../from_syn_restrain_N750mv_anton/from_anton_open.gro -r ../from_syn_restrain_N750mv_anton/from_anton_open.gro -n md.ndx -p Protein.top

       2         1         0      -0.316 nm         -0.700 nm
       3         1         0      -0.975 nm         -0.700 nm 
       4         1         0      -0.788 nm         -0.700 nm 
       5         1         0      -0.675 nm         -0.700 nm ==> -0.688  avcverage
       6         1         0       0.304 nm         -0.120 nm
       7         1         0      -0.335 nm         -0.120 nm
       8         1         0      -0.302 nm         -0.120 nm
       9         1         0      -0.122 nm         -0.120 nm ==>  -0.11375 average

*****************
check  inital:

  gmx_mpi grompp -f Protein_0-1000ns_1.0.mdp -o  Protein_0-1000ns_2020.tpr -c initial.gro  -r initial.gro -n md.ndx -p Protein.top

Pull group  natoms  pbc atom  distance at start  reference at t=0
       2         1         0       0.148 nm         -0.700 nm
       3         1         0       0.119 nm         -0.700 nm
       4         1         0       0.393 nm         -0.700 nm
       5         1         0      -0.051 nm         -0.700 nm ==> 0.15225
       6         1         0       0.926 nm         -0.120 nm
       7         1         0       0.892 nm         -0.120 nm
       8         1         0       1.066 nm         -0.120 nm
       9         1         0       1.120 nm         -0.120 nm ==> 1.001

   aim pull 210 down 0.8  to -0.7 

            213      1.0  to  0  

 for 213 all put rate   target  0.00  1.000/50000   0.00002

 for 210  

      2         1         0       0.148 nm       
      3         1         0       0.119 nm        use (0.15+0.700)/50000   0.000017

      4         1         0       0.393 nm            (0.4+0.700)/50000    0.000022
      5         1         0      -0.051 nm                 0.700/50000     0.000014


   ******  how lartge ??

   *  K:  200         ; kJ mol^-1 nm^-2   ==>  0.47801 kcal mol^-1 A^-2 
   *  K:  100         ; kJ mol^-1 nm^-2   ==>  0.239

         1005.6                                2.4
         2092.5                                5  





  #########################################

      condition summary 

  #########################################

 test


  gmx_mpi make_ndx -f Protein_pr_by_chain.tpr  -o md.ndx
 
   3 & r 286-288   <== ref
[ CZ_r_210 ]


  
   rename to group in mdp
   
 export PATH=/home/zgjia/Software/gromacs/gromacs2019/bin/:$PATH
 export LD_LIBRARY_PATH=/home/zgjia/Software/gromacs/gromacs2019//lib/:$LD_LIBRARY_PATH

  gmx_mpi grompp -f Protein_0-1000ns.mdp -o  Protein_0-1000ns_2020.tpr -c from_anton.gro -r from_anton.gro -n md.ndx  -p Protein.top

1 )  condition 1 

  #  cp initial.gro /home/zgjia/Umass_home/BK_channel/Simulation_hBK/Core_MT_from_syn/restrain_pull_200mV_tripull_0.6kcal/

    ref  /home/zgjia/Project/BK_channel/Simulation_hBK/6v3g_from_syn/200mV_pull/

############################################################



################

################   restrain  and  pull 


                                  RMSD?                                  align
                            _________________             _______________________________
                             S4             S5           P-helix   filter               S6
                       205       225   230      259     273     287     292    298     309       324
                        G               S        S       T       T       D      T       L

        Jianming also request TM 1 2 3  6  (split lower pore, upper pore ? )


   _________________             _______________________________
                                                                                                                                            S6
                      S1               S2               S3               S4                  S5              filter                    low       |         up
                  109     135     148     170      181      199    205       225         230      259       287     292          298       309        310       324
                   T       S       F       A        V        L      G         N           S        S         T       D            T         L          G         E
 real orien      down     up       up    down     down      up     up       down        down      up       down      up           up                           down
 sim  orien       up     down     down    up       up      down   down       up          up      down       up      down        down                            up


   ***************************

    idea is a test pullimng under electic environment 

 

   ***************************  

   ./Gen_restrain_less.py toppar/PROAB.itp  toppar/PROAB_rest.itp

########################33

    






############
###############################################################


  #  rsync -r ../Core_MT_close_PC_from_syn --progress  zhiguangjia_umass_edu@unity.rc.umass.edu:/home/zhiguangjia_umass_edu/Project/BK_channel/Simulation_hBK/Core_MT

  #  rsync --progress eq_from_syn.gro  zhiguangjia_umass_edu@unity.rc.umass.edu:/home/zhiguangjia_umass_edu/Project/BK_channel/Simulation_hBK/Core_MT/Core_MT_close_PC_from_syn

   rsync -r ../restrain_pull_0mV_210_213_sep_5kcal_steered_V3_1 --progress  pikes.ials.umass.edu:/home/zgjia/Project/BK_channel/Simulation_hBK/Core_MT_modeled/



###############################################################

############################################################################################################################################################

; 1st is V/nm  1000 mV/nm

# -120mv / (4 nm *1000) (or use system size?)  = 0.03

40/( 10.4 * 1000)

  150                      0.014423076923076924

  250   /( 10.4   * 1000)  0.02403846153846154

  500  /( 10.4   * 1000)   0.04807692307692308

  750 /( 10.4   * 1000)  0.07211538461538461


################################################################################

 00)                    preppare   # now trajectory treated in pikes

################################################################################


gmx_mpi trjconv -f ~/Pikes_home/BK_channel/Simulation_hBK/Core_MT_modeled/restrain_pull_0mV_210_213_sep_5kcal_steered_V3_1/Protein_0-1000ns.xtc -s Protein_pr_by_chain.tpr  -o sim_1ns.xtc -pbc mol  -dt 1000


 gmx_mpi trjconv -f  sim_1ns.xtc  -s Protein_pr_by_chain.tpr  -o 50ns.gro -pbc mol -tu ns -b 50 -e 50


 

  ** now run in 200 ns step to avoid unity harddisk error

    gmx_mpi trjconv -f  sim_1ns.xtc -s Protein_pr_by_chain.tpr -o 200ns.gro  -pbc mol  -tu ns -b 200 -e 200



###########################################   analysis

# vmd observation









########################################

    

    




~                                                                                                                                                                                             
~                                                                                                                                                                                             
~                                                                                                                                                                                             
~                                                                                                                                                                                             
~                                                                                                                                                                                             
~                                                                                                                                                                                             
~                                          
