#!//usr/bin/python

import MDAnalysis
import numpy
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

timestep      = float( sys.argv[4])
starting_time = float( sys.argv[5])  # simple add this namer to the first frame time, as our dcd is by 480 ns segment 

try:
    intitial_frame = int(sys.argv[6])
    final_frame    = int(sys.argv[7])   

except:
    intitial_frame = 0
    final_frame    = 99999999


print  " "
print  " our CTD is in wrong direction (to up) , so in real , S4 moves up "
print  "                                            in here , it moves down"
print  "  so e fild shoud go dopwn (S4 is all-positive charged) "
print  " "
print  "  so the result need *-1 before use "
print  " "
print  " "



print  "  NOTE!  insertion / deletion need make a local version as resid changed !!  "

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)



############################   for  inertion/deletion mutations, modify above range is enough   all  folowin part call this range


##############  check 207 210 213  CA

ofile207_CA  = open(str(outfilename + '_207_CA.xvg'),'w') # open file for writing
ofile210_CA  = open(str(outfilename + '_210_CA.xvg'),'w') # open file for writing
ofile213_CA  = open(str(outfilename + '_213_CA.xvg'),'w') # open file for writing

##############  check 207 210 213  CZ

ofile207_CZ  = open(str(outfilename + '_207_CZ.xvg'),'w') # open file for writing
ofile210_CZ  = open(str(outfilename + '_210_CZ.xvg'),'w') # open file for writing
ofile213_CZ  = open(str(outfilename + '_213_CZ.xvg'),'w') # open file for writing


#########  others

ofile153_CA  = open(str(outfilename + '_153_CA.xvg'),'w') # open file for writing
ofile167_CA  = open(str(outfilename + '_167_CA.xvg'),'w') # open file for writing
ofile186_CA  = open(str(outfilename + '_186_CA.xvg'),'w') # open file for writing
ofile219_CA  = open(str(outfilename + '_219_CA.xvg'),'w') # open file for writing

ofile153_CZ  = open(str(outfilename + '_153_CZ.xvg'),'w') # open file for writing
ofile167_CZ  = open(str(outfilename + '_167_CZ.xvg'),'w') # open file for writing
ofile186_CZ  = open(str(outfilename + '_186_CZ.xvg'),'w') # open file for writing
ofile219_CZ  = open(str(outfilename + '_219_CZ.xvg'),'w') # open file for writing


##############  check  S2-5
ofileS5  = open(str(outfilename + '_S5.xvg'),'w') # open file for writing
ofileS4  = open(str(outfilename + '_S4.xvg'),'w') # open file for writing
ofileS3  = open(str(outfilename + '_S3.xvg'),'w') # open file for writing
ofileS2  = open(str(outfilename + '_S2.xvg'),'w') # open file for writing
ofileS1  = open(str(outfilename + '_S1.xvg'),'w') # open file for writing






RCKchain = ['B','C','D','A']   ## may be wrong ....
VSDchain = ['A','B','C','D']  

# initalize value
framenumber = 0

print(' anton , when convert to dcd, all time lost , need restart')

for ts in targetprotein.trajectory:
   framenumber += 1



   if intitial_frame != 0  and ( framenumber <= intitial_frame or framenumber >= final_frame ):
      pass

   elif intitial_frame == 0 or ( framenumber > intitial_frame and framenumber < final_frame ):

    #  gateresidue = targetprotein.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")

      #x0 = gateresidue.center_of_mass()[0]
      #y0 = gateresidue.center_of_mass()[1]
    #  z0  = gateresidue.center_of_mass()[2]

    #  membrane    = targetprotein.select_atoms(" resname POP OPC OPE POPC POPE ")

      membrane    =  targetprotein.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")


      z0_membrane = membrane.center_of_mass()[2]

      ############################       CA part
 
      #  now  loop four  pair RCK-VSD    
      current_residue = 207
      ofile207_CA.write( str(ts.frame*timestep + starting_time ) + "    ")

      ## first, overall 

      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofile207_CA.write( str(-diff) + "    ")             

      for currect_VSDchain in VSDchain  :       

          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )             
       #  Clink_331to334   = targetprotein.select_atoms(" (resid %d:%d ) and segid %s "%(Clink_331to334_id[0], Clink_331to334_id[-1], currect_VSDchain) )

          diff          = Current_sel.center_of_mass()[2] -z0_membrane

          

          ofile207_CA.write( str(-diff) + "    ")

          # r                    = Clink_321.center_of_mass()-gateresidue.center_of_mass()
          # dist_321_filter     += numpy.linalg.norm(r)
          #  now there is question   how define the 
   
      ofile207_CA.write( "\n")

      ### next
      current_residue = 210
      ofile210_CA.write( str(ts.frame*timestep + starting_time ) + "    ")

      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofile210_CA.write( str(-diff) + "    ")

      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane

          ofile210_CA.write( str(-diff) + "    ")
      ofile210_CA.write( "\n")

      ### next
      current_residue = 213
      ofile213_CA.write( str(ts.frame*timestep  + starting_time) + "    ")

      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofile213_CA.write( str(-diff) + "    ")


      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane

          ofile213_CA.write( str(-diff) + "    ")
      ofile213_CA.write( "\n")



      ### next
      current_residue = 153
      ofile153_CA.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile153_CA.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile153_CA.write( str(-diff) + "    ")
      ofile153_CA.write( "\n")
      ### next
      current_residue = 167
      ofile167_CA.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile167_CA.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile167_CA.write( str(-diff) + "    ")
      ofile167_CA.write( "\n")
      ### next
      current_residue = 186
      ofile186_CA.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile186_CA.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile186_CA.write( str(-diff) + "    ")
      ofile186_CA.write( "\n")
      ### next
      current_residue = 219
      ofile219_CA.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile219_CA.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile219_CA.write( str(-diff) + "    ")
      ofile219_CA.write( "\n")


      ############################       CZ part

      current_residue = 207
      ofile207_CZ.write( str(ts.frame*timestep  + starting_time ) + "    ")

      Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofile207_CZ.write( str(-diff) + "    ")

      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) and segid %s "%(current_residue, currect_VSDchain) )
       #  Clink_331to334   = targetprotein.select_atoms(" (resid %d:%d ) and segid %s "%(Clink_331to334_id[0], Clink_331to334_id[-1], currect_VSDchain) )

          diff          = Current_sel.center_of_mass()[2] -z0_membrane

          ofile207_CZ.write( str(-diff) + "    ")

          # r                    = Clink_321.center_of_mass()-gateresidue.center_of_mass()
          # dist_321_filter     += numpy.linalg.norm(r)
          #  now there is question   how define the 

      ofile207_CZ.write( "\n")

      ### next
      current_residue = 210
      ofile210_CZ.write( str(ts.frame*timestep  + starting_time) + "    ")

      Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofile210_CZ.write( str(-diff) + "    ")

      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane

          ofile210_CZ.write( str(-diff) + "    ")
      ofile210_CZ.write( "\n")

      ### next
      current_residue = 213
      ofile213_CZ.write( str(ts.frame*timestep  + starting_time) + "    ")

      Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofile213_CZ.write( str(-diff) + "    ")

      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane

          ofile213_CZ.write( str(-diff) + "    ")
      ofile213_CZ.write( "\n")

      ### next
      current_residue = 153
      ofile153_CZ.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CG ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile153_CZ.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CG  ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile153_CZ.write( str(-diff) + "    ")
      ofile153_CZ.write( "\n")
      ### next
      current_residue = 167
      ofile167_CZ.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile167_CZ.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CZ Cz ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile167_CZ.write( str(-diff) + "    ")
      ofile167_CZ.write( "\n")
      ### next
      current_residue = 186
      ofile186_CZ.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CG ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile186_CZ.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CG ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile186_CZ.write( str(-diff) + "    ")
      ofile186_CZ.write( "\n")
      ### next
      current_residue = 219
      ofile219_CZ.write( str(ts.frame*timestep + starting_time ) + "    ")
      Current_sel   = targetprotein.select_atoms(" (resid %d and name CD ) "%(current_residue) )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofile219_CZ.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid %d and name CD ) and segid %s "%(current_residue, currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofile219_CZ.write( str(-diff) + "    ")
      ofile219_CZ.write( "\n")


      ############################    S4 part


      ofileS4.write( str(ts.frame*timestep  + starting_time) + "    ")

      Current_sel   = targetprotein.select_atoms(" ( resid 205:225 and name CA Ca  ) " )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane

      ofileS4.write( str(-diff) + "    ")

      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid 205:225 and name CA Ca  ) and segid %s "%( currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofileS4.write( str(-diff) + "    ")
      ofileS4.write( "\n")

      # next 
      ofileS1.write( str(ts.frame*timestep  + starting_time) + "    ")
      Current_sel   = targetprotein.select_atoms(" ( resid 111:136 and name CA Ca  ) " )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofileS1.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid 111:136 and name CA Ca  ) and segid %s "%( currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofileS1.write( str(-diff) + "    ")
      ofileS1.write( "\n")
      # next
      ofileS2.write( str(ts.frame*timestep  + starting_time) + "    ")
      Current_sel   = targetprotein.select_atoms(" ( resid 148:169 and name CA Ca  ) " )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofileS2.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid 148:169 and name CA Ca  ) and segid %s "%( currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofileS2.write( str(-diff) + "    ")
      ofileS2.write( "\n")
      # next
      ofileS3.write( str(ts.frame*timestep  + starting_time) + "    ")
      Current_sel   = targetprotein.select_atoms(" ( resid 174:199 and name CA Ca  ) " )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofileS3.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid 174:199 and name CA Ca  ) and segid %s "%( currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofileS3.write( str(-diff) + "    ")
      ofileS3.write( "\n")
      # next
      ofileS5.write( str(ts.frame*timestep  + starting_time) + "    ")
      Current_sel   = targetprotein.select_atoms(" ( resid 230:256 and name CA Ca  ) " )
      diff          = Current_sel.center_of_mass()[2] -z0_membrane
      ofileS5.write( str(-diff) + "    ")
      for currect_VSDchain in VSDchain  :
          Current_sel   = targetprotein.select_atoms(" (resid 230:256 and name CA Ca  ) and segid %s "%( currect_VSDchain) )
          diff          = Current_sel.center_of_mass()[2] -z0_membrane
          ofileS5.write( str(-diff) + "    ")
      ofileS5.write( "\n")


      if float(ts.frame)%100 == 0 :
          print  ' done frame ', framenumber





