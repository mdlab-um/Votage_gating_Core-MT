#!//usr/bin/python

import MDAnalysis
import numpy
import numpy as np
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
from collections import OrderedDict


structurefile = sys.argv[1]
targettrajectory = sys.argv[2]

outfilename   = sys.argv[3]

try:
    timestep     = float(sys.argv[4])
    first_frame  = int(sys.argv[5])
    skip         = float(sys.argv[6])
    final_frame  = int(sys.argv[7])

except:
    timestep     = 0.24   # suppose anton input  unit  ns
    first_frame  = 0
    skip         = 1
    final_frame = 999999999999


print  " "
print  " our CTD is in wrong direction (to up) , so in real , S4 moves up "
print  "                                            in here , it moves down"
print  "  so e fild shoud go dopwn (S4 is all-positive charged) "
print  " "
print  "  so the result need *-1 before use "
print  " "
print  " "


print  "  NOTE!  insertion / deletion need make a local version as resid changed !!  "


#               we need
#
#                       |       move                                rj already alinged on this, 1st frame should be ref?
#                       |    _________________             _______________________________
#                             S4             S5           P-helix   filter               S6
#                       205       225   230      259     273     287     292    298     309       324
#                        G               S        S       T       T       D      T       L
#
#        Jianming also request TM 1 2 3  6  (split lower pore, upper pore ? )
# 
#    
#   _________________             _______________________________                                                
#                                                                                                                                            S6
#                      S1               S2               S3               S4                  S5              filter                    inner       |      lower
#                  109     135     148     170      181      199    205       225         230      259       287     292          298       309        310       324
#                   T       S       F       A        V        L      G         N           S        S         T       D            T         L          G         E
# real orien      down     up       up    down     down      up     up       down        down      up       down      up           up                           down
# sim  orien       up     down     down    up       up      down   down       up          up      down       up      down        down                            up 
#
#  
#
#
#
#
#
#

#      we need set a reference value ....
#
#          horizantal movement should be define as the change in distance to filter ?

######################################     reference  part

reference     = MDAnalysis.Universe(structurefile, structurefile)

for ts in reference.trajectory:
    # should be single frame 

    filter_single =reference.select_atoms(" ( (resid 287:292 ) and name CA) ") 
    filter_single_ref_x = filter_single.center_of_mass()[0]
    filter_single_ref_y = filter_single.center_of_mass()[1] 
    
 
    S1 = reference.select_atoms(" ( (resid 109:135 ) and name CA) ")
    S2 = reference.select_atoms(" ( (resid 148:170 ) and name CA) ")
    S3 = reference.select_atoms(" ( (resid 181:199 ) and name CA) ")
    S4 = reference.select_atoms(" ( (resid 205:225 ) and name CA) ")
    S4_lower  = reference.select_atoms(" ( (resid 214:225 ) and name CA) ")
    S5 = reference.select_atoms(" ( (resid 230:259 ) and name CA) ") 
    S6 = reference.select_atoms(" ( (resid 298:324 ) and name CA) ")

    S6_inner = reference.select_atoms(" ( (resid 298:309 ) and name CA) ")
    S6_lower = reference.select_atoms(" ( (resid 310:324 ) and name CA) ")


    S1_ref_x = S1.center_of_mass()[0]
    S1_ref_y = S1.center_of_mass()[1]
    S1_ref_z = S1.center_of_mass()[2]
    S2_ref_x = S2.center_of_mass()[0]
    S2_ref_y = S2.center_of_mass()[1]
    S2_ref_z = S2.center_of_mass()[2]
    S3_ref_x = S3.center_of_mass()[0]
    S3_ref_y = S3.center_of_mass()[1]
    S3_ref_z = S3.center_of_mass()[2]
    S4_ref_x = S4.center_of_mass()[0]
    S4_ref_y = S4.center_of_mass()[1]
    S4_ref_z = S4.center_of_mass()[2]
    S4_lower_ref_x = S4_lower.center_of_mass()[0]
    S4_lower_ref_y = S4_lower.center_of_mass()[1]
    S4_lower_ref_z = S4_lower.center_of_mass()[2]

    S5_ref_x = S5.center_of_mass()[0]
    S5_ref_y = S5.center_of_mass()[1]
    S5_ref_z = S5.center_of_mass()[2]
    S6_ref_x = S6.center_of_mass()[0]
    S6_ref_y = S6.center_of_mass()[1]
    S6_ref_z = S6.center_of_mass()[2]

    S6_inner_x = S6_inner.center_of_mass()[0]
    S6_inner_y = S6_inner.center_of_mass()[1]
    S6_inner_z = S6_inner.center_of_mass()[2]
    S6_lower_x = S6_lower.center_of_mass()[0]
    S6_lower_y = S6_lower.center_of_mass()[1]
    S6_lower_z = S6_lower.center_of_mass()[2]


  #  S1_ref_filter_distance = math.sqrt( ((S1_ref_x-filter_single_ref_x)**2)+((S1_ref_y-filter_single_ref_y)**2) )

    ###  now define the tilt angle (vector)     *******************************   rememer  our  system is upside down, so reverse PDB up/down order   !!!!!


#################################   reference  end 

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)

totalframe = len(targetprotein.trajectory)

filter_horizon_move  = open(str(outfilename + '_filter_2D.xvg'),'w') # open file for writing
S1_horizon_move  = open(str(outfilename + '_S1_2D.xvg'),'w') # open file for writing
S2_horizon_move  = open(str(outfilename + '_S2_2D.xvg'),'w') # open file for writing
S3_horizon_move  = open(str(outfilename + '_S3_2D.xvg'),'w') # open file for writing
S4_horizon_move  = open(str(outfilename + '_S4_2D.xvg'),'w') # open file for writing
S4_lower_horizon_move  = open(str(outfilename + '_S4_lower_214_225_2D.xvg'),'w') # open file for writing
S5_horizon_move  = open(str(outfilename + '_S5_2D.xvg'),'w') # open file for writing
S6_horizon_move  = open(str(outfilename + '_S6_2D.xvg'),'w') # open file for writing
S6_inner_horizon_move  = open(str(outfilename + '_S6_inner_2D.xvg'),'w') # open file for writing
S6_lower_horizon_move  = open(str(outfilename + '_S6_lower_2D.xvg'),'w') # open file for writing


RCKchain = ['B','C','D','A']   ## may be wrong ....
VSDchain = ['A','B','C','D']  

# initalize value
framenumber = 0

for ts in targetprotein.trajectory:
   framenumber += 1

 #  if intitial_frame != 0  and ( framenumber <= intitial_frame or framenumber >= final_frame ):
 #     pass
   if (ts.frame == 0 or  float(ts.frame)%skip == 0 ) and (  ts.frame >= first_frame and ts.frame <= final_frame ) :

      currenttime = framenumber*timestep     

      filter_single =targetprotein.select_atoms(" ( (resid 287:292 ) and name CA) ")
      filter_single_current_x = filter_single.center_of_mass()[0]
      filter_single_current_y = filter_single.center_of_mass()[1]
      filter_single_current_z = filter_single.center_of_mass()[2]

      S1 = targetprotein.select_atoms(" ( (resid 109:135 ) and name CA) ")
      S2 = targetprotein.select_atoms(" ( (resid 148:170 ) and name CA) ")
      S3 = targetprotein.select_atoms(" ( (resid 181:199 ) and name CA) ")
      S4 = targetprotein.select_atoms(" ( (resid 205:225 ) and name CA) ")
      S4_lower = targetprotein.select_atoms(" ( (resid 214:225 ) and name CA) ")

      S5 = targetprotein.select_atoms(" ( (resid 230:259 ) and name CA) ")
      S6 = targetprotein.select_atoms(" ( (resid 298:324 ) and name CA) ")
      S6_inner = targetprotein.select_atoms(" ( (resid 298:309 ) and name CA) ")
      S6_lower = targetprotein.select_atoms(" ( (resid 310:324 ) and name CA) ")


      S1_current_x = S1.center_of_mass()[0]
      S1_current_y = S1.center_of_mass()[1]
      S1_current_z = S1.center_of_mass()[2]
      S2_current_x = S2.center_of_mass()[0]
      S2_current_y = S2.center_of_mass()[1]
      S2_current_z = S2.center_of_mass()[2]
      S3_current_x = S3.center_of_mass()[0]
      S3_current_y = S3.center_of_mass()[1]
      S3_current_z = S3.center_of_mass()[2]
      S4_current_x = S4.center_of_mass()[0]
      S4_current_y = S4.center_of_mass()[1]
      S4_current_z = S4.center_of_mass()[2]
      S4_lower_current_x = S4_lower.center_of_mass()[0]
      S4_lower_current_y = S4_lower.center_of_mass()[1]
      S4_lower_current_z = S4_lower.center_of_mass()[2]

      S5_current_x = S5.center_of_mass()[0]
      S5_current_y = S5.center_of_mass()[1]
      S5_current_z = S5.center_of_mass()[2]
      S6_current_x = S6.center_of_mass()[0]
      S6_current_y = S6.center_of_mass()[1]
      S6_current_z = S6.center_of_mass()[2]
      S6_inner_current_x = S6_inner.center_of_mass()[0]
      S6_inner_current_y = S6_inner.center_of_mass()[1]
      S6_inner_current_z = S6_inner.center_of_mass()[2]
      S6_lower_current_x = S6_lower.center_of_mass()[0]
      S6_lower_current_y = S6_lower.center_of_mass()[1]
      S6_lower_current_z = S6_lower.center_of_mass()[2]

    #  S1_current_filter_distance = math.sqrt( ((S1_current_x-filter_single_current_x)**2)+((S1_current_y-filter_single_current_y)**2) )
      filter_horizon_move.write( str(currenttime) + "  " + str(filter_single_current_x) + "  "  + str(filter_single_current_y) + "  "  + str(filter_single_current_z) + "  "  + "\n")
      S1_horizon_move.write( str(currenttime) + "  " + str(S1_current_x) + "  "  + str(S1_current_y) + "  "  + str(S1_current_z) + "  "  + "\n")
      S2_horizon_move.write( str(currenttime) + "  " + str(S2_current_x) + "  "  + str(S2_current_y) + "  "  + str(S2_current_z) + "  "  + "\n")
      S3_horizon_move.write( str(currenttime) + "  " + str(S3_current_x) + "  "  + str(S3_current_y) + "  "  + str(S3_current_z) + "  "  + "\n")
      S4_horizon_move.write( str(currenttime) + "  " + str(S4_current_x) + "  "  + str(S4_current_y) + "  "  + str(S4_current_z) + "  "  + "\n")
      S4_lower_horizon_move.write( str(currenttime) + "  " + str(S4_lower_current_x) + "  "  + str(S4_lower_current_y) + "  "  + str(S4_lower_current_z) + "  "  + "\n")

      S5_horizon_move.write( str(currenttime) + "  " + str(S5_current_x) + "  "  + str(S5_current_y) + "  "  + str(S5_current_z) + "  "  + "\n")
      S6_horizon_move.write( str(currenttime) + "  " + str(S6_current_x) + "  "  + str(S6_current_y) + "  "  + str(S6_current_z) + "  "  + "\n")

      S6_lower_horizon_move.write( str(currenttime) + "  " + str(S6_lower_current_x) + "  "  + str(S6_lower_current_y) + "  "  + str(S6_lower_current_z) + "  "  + "\n")
      S6_inner_horizon_move.write( str(currenttime) + "  " + str(S6_inner_current_x) + "  "  + str(S6_inner_current_y) + "  "  + str(S6_inner_current_z) + "  "  + "\n")

      #  now  loop four  pair RCK-VSD    

      if float(ts.frame)%100 == 0 :
          print  ' done frame ', framenumber





