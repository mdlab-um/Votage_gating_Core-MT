set term postscript eps color enhanced
set size square
set output "Core_MT_modeled_750mv_2_9500-10000ns_sep_25_2024.eps"

#set pm3d interpolate 0,0 
    #  interpolate is smooth

set pm3d at s
set pm3d map scansautomatic flush begin noftriangles nohidden3d implicit corners2color mean

set pm3d interpolate 0,0
#         interplote addition 1 point  for smooth
#

# jianhan /home/jianhanc/work/tutorial/charmm-vacuum/plot_cov.gnu
# set pm3d at s
# set pm3d scansautomatic flush begin noftriangles nohidden3d implicit corners2color mean  
set view map


set xtics 120
set ytics 120
set mxtics 5
set mytics 5
#set mxtics 2  divided major tics by 2?
#set grid back 64 56
set xrange [20:80]
set yrange [60:120]

#set label "F" at 78.7,80.4 tc rgb "black" font ",30" front 
#set label "6" at 70.7,85.0 tc rgb "black" font ",30" front     
#set label "5" at 62.0,79.6 tc rgb "black" font ",30" front     
#set label "4" at 52.6,77.9 tc rgb "black" font ",30" front     
#set label "3" at 44.7,86.2 tc rgb "black" font ",30" front     
#set label "2" at 40.2,79.8 tc rgb "black" font ",30" front     
#set label "1" at 46.9,72.7 tc rgb "black" font ",30" front     

#set grid ytics lc rgb "#bbbbbb" lw 0.5 lt 0
#set grid xtics lc rgb "#bbbbbb" lw 0.5 lt 0

#set contour 
#set style increment user
#set cntrparam levels incr -0.9,0.15,0.9

set cbrange [0:5.5] #
#set palette defined (0 "violet", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "orange", 6 "red", 7 "white")

set palette defined (0 "dark-blue", 0.5 "blue",  1 "light-blue", 1.5 "cyan", 2 "green", 2.5 "greenyellow", 3  "yellow", 3.5 "orange", 4 "red", 4.5 "light-red", 5 "pink", 5.5 'white')


set cbtics 1
set format cb "%3.1f"

splot 'Core_MT_modeled_750mv_2_9500-10000ns_sep_25_2024.dat' t ''

set out


