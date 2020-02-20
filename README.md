# IMPES_3D
PURPOSE: Provide simulation solutions to 1D, 2D and 3D 2-PHASE reservoir problem.    
FEATURES:
- Works for Any Shape &amp; Size 1D, 2D &amp; 3D reservoir. 
- Matrixes are created in REAL-TIME.  Works for Multiple Layer Reservoir. 
- Calculates Areas based on direction in real-time 
- Calculates Volumes for each grid-block in real-time 
- Works with any size grid-blocks where DX != DY != DZ 
- Uses Up-winding technique to calculate inter-block relative permeability 
- Uses Harmonic Average to calculate inter-block transmissibility 
- Calculates Gravity for 3D simulations 
- Calculates J-Wells production based on radius, skin &amp; PWF 
- Handles Massive Output Data by dividing output files into smaller ones
- Works with Variable Absolute Permeabilities:    X,Y and Z Absolute Permeabilities can be assigned as initial conditions in input file. 
- X,Y and Z Absolute Permeabilities can be read from a Geostatistical Grid 
- Works with any number of user defined injectors and respective initial conditions
- Works with any number of user defined producers and respective initial conditions including J-wells 
INPUT:  
- All Initial Conditions are set in excel file 
- Number of injectors and Producers are set in excel file 
- Initial conditions for all Injectors &amp; Producers are set in excel file 
- Initial conditions for “J-Wells” (i.e. production data) are set in excel file  
OUTPUT:  
- Output 1: List of Pressure by grid-block &amp; iteration (CSV file)
- Output 2: List of Saturation by grid-block &amp; iteration (CSV file)
- Output 3 : List of Production by J-well &amp; iteration f (CSV file) 

NOTE: Since the number of data can be massive, all output files are divided into smaller ones to allow softwares like excel and Google spread sheet to being able to open them. 

Limitations: Variable Density, variable Porosity and Capillary pressure are not implemented in this simulator. 
 
THE PROGRAM:  
 
Main Algorithm: Solutions are found by using IMPES to implicitly solve for Pressure (equation a) and explicitly solve for Saturation (equations b). 
 
(a) (J + T + B/dt) P(n+1) = (B/dt) * P(n) + Q + G 
(b) S = (Sw - Swi)/ 1 - Swi -Swr and krw = .2S**3   and kro = (1-S)**3 
 
 
Transmissibility Matrix Configuration: Matrix is built in one single pass by acknowledging the fact that: 
 
(1) 2D is y * 1D  
(2) 3D is z * 2D 
 
where x, y and z represent the number of grid-blocks in each direction in space. This allowed to create the algorithm below: 
 
for i in range (0, Total # of blocks in Z-direction) 
 
    if i == 0, isFront = TRUE 
    if i == Z, isBack = TRUE 
 
    for j in range (0, Total # of blocks in Y-direction) 
 
         if j == 0, isTop = TRUE 
         if j == Y, isBottom = TRUE 
 
         for k in range (0, Total # of blocks in X-direction ) 
 
              if k == 1 isRight = TRUE 
              if k == x isLEFT = TRUE 
 
             Compute_Matrix() 
 
This simple method allows to configure transmissibility matrixes for any size and shape reservoir with absolute simplicity. 
 
Functions: 
 
- Brooks and Corey 
Brook and Corey type relationship is used when calculating relative permeability for 2-phase water and oil simulations 
 
- SPE relative permeability Model 1 
SPE relative permeability Model 1 was used when calculating relative permeability for 2-phase gas and oil simulations (Blunt, 2001) 
 
- Harmonic Mean 
Harmonic Mean was used to calculate inter-block transmissibility 
 
- Up-Winding or Up-Streaming 
Up-Winding was used to calculate the inter-block relative permeability in x, y and z-direction respectively 
 
- Dynamic Area and Volume Calculations 
Areas are calculated dynamically for x, y and z-direction respectively. Volumes are calculated in dynamically for each grid-block 
 
- Absolute Permeabilities 
X, Y and Z permeabilites can be assigned as Initial Conditions in Excel file or they can be read from a Geostatistical Grid. The name of the file containing the Geostatiscal Grid can be specified in the excel file containing initial conditions 
 
INPUT FILE: 
File name: Input.xlsx 
Total Sheets: 3 
Sheet Names: Initial, Injectors, Producers, and J 
Screen Shots below 
 
OUTPUT FILES: 
File names: 
(1) Pressure-Output-# 
(2) Saturation-Output-# 
(3) Jflow (where # in the file name indicates starting iteration, and J-flow represents the flow rate at that iteration in bbd) 
