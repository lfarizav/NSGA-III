#!/bin/bash
j=1
for i in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.97
  do 
     echo "The seed for dtlz4_3D is $i (j is $j)"
     ./nsga2r $i <input_data/dtlz4_3d.in
     cp plot.out fronts/dtlz4/3obj/plot_3d_$j.out
      cp plot_IGD.out fronts/dtlz4/3obj/plot_3d_IGD_$j.out
      cp plot_convergence.out fronts/dtlz4/3obj/plot_3d_convergence_$j.out
     j=$((j+1))
 done
j=1
for i in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.97
  do 
     echo "The seed for dtlz4_5D is $i (j is $j)"
     ./nsga2r $i <input_data/dtlz4_5d.in
     cp plot.out fronts/dtlz4/5obj/plot_5d_$j.out
     cp plot_IGD.out fronts/dtlz4/5obj/plot_5d_IGD_$j.out
     cp plot_convergence.out fronts/dtlz4/5obj/plot_5d_convergence_$j.out
     cp plot_pc.out fronts/dtlz4/5obj/plot_5d_pc_$j.out
     j=$((j+1))
 done
j=1
for i in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.97
  do 
     echo "The seed for dtlz4_8D is $i (j is $j)"
     ./nsga2r $i <input_data/dtlz4_8d.in
     cp plot.out fronts/dtlz4/8obj/plot_8d_$j.out
     cp plot_IGD.out fronts/dtlz4/8obj/plot_8d_IGD_$j.out
     cp plot_convergence.out fronts/dtlz4/8obj/plot_8d_convergence_$j.out
     cp plot_pc.out fronts/dtlz4/8obj/plot_8d_pc_$j.out
     j=$((j+1))
 done
j=1
for i in 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.97
  do 
     echo "The seed for dtlz4_10D is $i (j is $j)"
     ./nsga2r $i <input_data/dtlz4_10d.in
     cp plot.out fronts/dtlz4/10obj/plot_10d_$j.out
     cp plot_IGD.out fronts/dtlz4/10obj/plot_10d_IGD_$j.out
     cp plot_convergence.out fronts/dtlz4/10obj/plot_10d_convergence_$j.out
     cp plot_pc.out fronts/dtlz4/10obj/plot_10d_pc_$j.out
     j=$((j+1))
 done
