Code Usage
#############
This section demonstrates how to use the code to solve USGS/SECE benchmark problems TPV3, TPV5, TPV6, and TPV7. 


Input files
***************
The input file that we have created performs an automatic feed of the data into code during the execution. One can create a new input file according to the fault plane dimensions, material properties, loads and output locations following the example input file given below.

An example input file is given as follows::

   **Problem_Info**
   TPV_No   : Enter the TPV problem being solved (enter TPV3, TPV5, TPV6, and TPV7). Set as TEST if testing any other problem.
   TPV3
   **Fault Geometry**
   L1rpt    : The length of slip-weakening interface/rupture zone  size in X1 direction in mts
   30000.d0
   L3rpt    : The length of slip-weakening interface/rupture zone size in X3 direction in mts
   15000.d0
   Xr_Ratio: Ratio of (rupture zone + barrier zone)  to rupture zone. i.e., (rupture zone + barrier zone)/rupture zone
   2.56d0
   Z_SYMM   : Specify 'T' if symmetric boundary conditions need to used else give 'F'. For TPV3 - F  , and for TPV5,TPV6,TPV7 - T.
   F
   **Mesh Info**
   nele1    :  (nele1) Total number of elements in X1-direction (Should be the powers of 2)
   768
   nele3    : Total number of elements in X3-direction (Should be the powers of 2)
   384
   **Material Info**
   Csm      : The shear wave speed of lower half-space  in m/s
   3464d0
   Cdm      : Dilatational wave speed of lower half-space  in m/s
   6000d0
   Rom      : The density of lower half-space  in kg/m³
   2670d0
   nu       : Poisson's ratio
   0.25d0
   csratio  : csratio of lower half-space to top half-space
   1.0d0
   cdratio  : cdratio of lower half-space to top half-space
   1.0d0
   Roratio  : Roratio of lower half-space to top half-space
   1.0d0
   **Simulation Info**
   Tend     : Total rupture duration time in seconds
   12d0
   betaa    : Courant parameter (CFL parameter)
   0.3464d0
   gamma    : Non-dimensionalised time step (Refinement in Kernel interpolation.)
   0.005d0
   **Kernel Info**
   input    : Give input as 0 if kernels to be calculated and stored in folder ./Kernels
   input    : Give input as 1 if kernels are calculated a priori and stored in folder ./Kernels
   1
   **Background Frictional Properties**
   mus0     :  Static coefficient of friction
   0.677d0
   mur0     :  Dynamic coefficient of friction
   0.525d0
   DeltaC   : Critical slip-weakening distance 
   0.4d0
   ** Background Stresses** T0bg(1),T0bg(2),T0bg(3) in N/m^2
   70d6  -120d6   0d0
   **Nucleation Zone List** (Only rectangular and square nucleation zones are modeled).  Note: x12>x11 and x32>x31
   No_neucles : No.of Nucleation Zones. If no Zones, give one line as 0,0,0,0,T0bg(1),T0bg(2),T0bg(3). If yes, add lines with x11,x31,x12,x32,tau1, tau2, tau3.
   1
   -1.5d3   -1.5d3  1.5d3   1.5d3  81.6d6 -120d6  0d6
   **Asperity Zone List** (Only circular asperities are modeled)
   No_Asps : No.of Asperity Zones. If no Zones, No_Asps = 0 add nothing. If yes, add No_Asps lines with Asps_x, Asps_y,Asps_radii,Asps_mus0,Asps_mur0.
   0
   ** Coordinates of output station points** x_i Y_i
   NStations: No. of Station points. Followed by coordinates of the station points in metres
   2
   0d3   3d3
   7.5d3 0d3


Each parameter can be changed and a new input file can be fed for solving problems.



Kernel Data
***************

We have provided the precomputed Kernels for USGS/SCEC benchmark problems TPV3, TPV5, TPV6, and TPV7 with the mesh size that we have choosen to solve the problem. 

However, if you want to compute the kernels, a separate directory with source code is provided. One can generate the Kernels using the following steps::

  cd TBSBIEM
  cd src_Kernels
  make
  cd ..
  ./Kernel_Check
  
or one can generate the Kernels during the simulation as well.  
  
Just provide '0' at `Pre_Kernels` prompt in the input file, the code will generate and store in a directory given like `TPV5_Kernels` for example.


**Note:** If you provide '1' at `Pre_Kernels` prompt in the input file, the code will search for kernels with no kernels generated. Thus be carefull if you need to generate the kernels, provide '0' at `Pre_Kernels` prompt in the input file.


Running a Simulation
*********************
The code is provided with several input files to solve USGS/SCEC benchmark problems. Namely, TPV3, TPV5, TPV6 and TPV7. One can solve these benchmark problems changing a problem name in **TPV_Problem.in** file.

For example, if you want to solve TPV6 problem then store the string::

   TPV6

in the file **Test_Problem.in**

Compile and execute the code::

      cd TBSBIEM
      cd src
      make
      cd ..
      ./TBSBIEM-v1.1.0


The problem TPV6 will be solved and the data will be stored in the directory './data'.

Post-Processing
********************* 
For a quick plotting one can use the gnuplot script given below to plot the contour plots on the fault plane as::

   for i in {0001..0015}; do    gnuplot -e "set terminal jpeg; set hidden3d; set xlabel 'x1 (km)'; set ylabel 'x3 (km)'; set zlabel 'Slip (m)'; set xrange [-15:15]; set yrange [-7.5:7.5];   set zrange [0.0:0.5]; set cbrange [0.0:0.5]; set view map; splot './data/TPV3_Out$i.dat' u 2:3:4 ps 0.1 palette" > Slip_Top$i.jpeg; done 
   
The example script generate a contour plot of rupture front with duration intervel of 0.5 Sec on  the fault plan for TPV3 benchmark problem.

Or we have given gnuplot script in './Post_Processing' directory using which one can generate the contour plot with time intervel of 1 Sec. One need to change the variable *Problem_no* to respective TPV problem required to plot out of 3,5,6 and 7.

One can create an interesting video using ffmpeg cammand as::

   ffmpeg -r 10 -i Slip_Top%04d.jpeg  -vf "fps=10" Slip_Top.mp4
   
**Note:** To plot these figure and the video generation one needs to install gnuplot and ffmpeg. The installation commands are as follows::

   sudo apt update
   sudo apt upgrade
   sudo apt install gnuplot  
   sudo apt install ffmpeg -y


Comparison of Results
----------------------------
We have provided GNUPLOT scripts and the benchmark data of MDBSI code for TPV3, TPV5, TPV6, and TPV7 to plot the various field variables at specific station points. The plots generated from the scripts provides the comparison of results from TBSBIEM and MDSBI.
