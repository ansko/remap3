#!/usr/bin/env gnuplot 


set terminal png


#plot 3795, "brief_results/11/11.lammps" u 1:2 w l t "modifier-modifier"
#plot -12777, "brief_results/11/11.lammps" u 1:3 w l t "modifier-polymer"
#plot -1351674, "brief_results/11/11.lammps" u 1:4 w l t "polymer-polymer"
#plot -1360656, "brief_results/11/11.lammps" u 1:5 w l t "soft-soft"


#plot 3795, "log.lammps" u 1:2 w l t "modifier-modifier
#plot -12777, "log.lammps" u 1:3 w l t "modifier-polymer"
#plot -1351674, "log.lammps" u 1:4 w l t "polymer-polymer"
plot -1360656, "log.lammps" u 1:5 w l t "soft-soft"
