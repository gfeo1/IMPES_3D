#not /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 10:07:34 2018

@author: giuseppefeo
"""

#/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 22 18:09:28 2018

@author: giuseppefeo
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:53:41 2018

@author: giuseppefeo
"""

import numpy as np
import pandas as pd
import math

#### ===================================================================== ####
#                     IMPES RESERVOIR SIMULATOR
#                      (1D & 2D & 3D) 2-PHASE
#### ===================================================================== ####

    
#### ========================  IMPORT INITIAL CONDITIONS ================= ####
# File Input.xlsx contains all the initial variables and conditions
#

# Assign spreadsheet filename to `file`
input_file = 'Input.xlsx'

# Load spreadsheet
input_xl = pd.ExcelFile(input_file)

# Load sheet Initial into a DataFrame by name: Initial_df
Initial_df = input_xl.parse('Initial')
    
isEnd = False
i = 0
    
while not isEnd:
    
    s = Initial_df.at[i,'Variable']
    if s == "DX":      #1
       DX = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "DY":      #2
       DY = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "XBLOCKS": #4
       XBLOCKS = int(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "YBLOCKS": #5
       YBLOCKS = int(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "PHI":    #6
       PHI = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "UW":     #7
       UW = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "SWI":    #8
       SWI = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "SWR":    #9
       SWR = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "CFW":    #10
       CFW = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "CFO":    #11
       CFO = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "CT":     #12
       CT = float(Initial_df.at[i,'Value'])
       i = i+1
       continue 
    if s == "P1":     #13
       P1 = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "S1":     #14
       S1 = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "BW":     #15
       BW = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "BO":     #16
       BO = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "DT":     #19
       DT = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "Ntime":  #20
       Ntime = Initial_df.at[i,'Value']
       Ntime = int(Ntime)
       i = i+1
       continue
    if s == "UO":     #21
       UO = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "ZBLOCKS":  #22
       ZBLOCKS = int(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "KX":     #23
       KX = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "KY":     #24
       KY = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "KZ":     #25
       KZ = float(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "IsGravity":     #25
       isGravity = int(Initial_df.at[i,'Value'])
       i = i+1
       continue
    if s == "Kfilename":
       Kfilename = Initial_df.at[i,'Value']
       i = i + 1
       continue
    if s == "end":
       IsEnd = True
       i = i+1
       break
    else:
       print("Variable No Define", Initial_df.at[i,'Variable'])
       exit()
       i = i+1
       
input_xl.close()      
    
#Calculate other Variables based on Initial Conditions
       
BlockTotal = XBLOCKS * YBLOCKS * ZBLOCKS        #Total Number of GridBlocks
LBLOCKS = int(XBLOCKS * YBLOCKS)                # GridBlocks per layer
L = DX * XBLOCKS                                #Total Lengh
RE = 0.2 * DX

#Read Geostatistical Grid (KData Array) 

if Kfilename != "NONE":
   input_Kdata = 'Kdata.xlsx'

   Kdata_xl = pd.ExcelFile(input_Kdata)

   #if BlockTotal != 2000:
   #   print("It does not match data")
   #   exit()
   
   Kdata_df = Kdata_xl.parse('KData')
   Kdata = np.zeros((BlockTotal,3))

   for i in range(0,BlockTotal):
    
       Kdata[i,0] = float(Kdata_df.at[i,'Kx'])
       Kdata[i,1] = float(Kdata_df.at[i,'Ky'])
       Kdata[i,2] = float(Kdata_df.at[i,'Kz'])
    
   Kdata_xl.close()
   
else:   

   ## If not Geostatistical Grid Kdata provided then 
   ## create Karray with values provided in Input File
   
   Kdata = np.zeros((BlockTotal,3))
   for i in range(0,BlockTotal): 
       Kdata[i,0] = KX
       Kdata[i,1] = KY
       Kdata[i,2] = KX


## Read DZ Values from Excel file
#NOTE: A = DY * DZ   #Area of Flow
#NOTE: V = A * DX    #Volume of each block

DZ_df = input_xl.parse('DZ')
myarray = [ 0 for row in range(BlockTotal)]
#buckets = [[0 for col in range(5)] for row in range(10)]
DZ = np.asarray(myarray)

Zfmatrix = np.zeros((BlockTotal,1))   ## Distance from top to center of gridblock
Vmatrix = np.zeros((BlockTotal,1))    ## Gridblock Volume
i = 0

for i in range(0,ZBLOCKS):
    
    t = float(DZ_df.at[i,'DZ'])
    
    if t < 0 and i != ZBLOCKS-1:
       print("DZ is missing layers")
       exit()
        
    if i == 0:
       ## First layer distance
       depth = t/2   
    else: 
       depth = depth + t  
    
    start = LBLOCKS * i
    for j in range(0,LBLOCKS):
       Zfmatrix[j + start] = depth
       DZ[j+start] = t
       Vmatrix[j + start] = DX * DY * t

Zfmatrix = np.asmatrix(Zfmatrix)
Vmatrix = np.asmatrix(Vmatrix) 

#Permeability of Water
krw1 = 0 #Initial Value of Water Permeabiliy
Krwarray = np.zeros((Ntime,BlockTotal))
for i in range (0, BlockTotal):
    Krwarray[0,i] = krw1

#Permeability of Oil
kro1 = 1 #Initial Permeability of Oil
Kroarray = np.zeros((Ntime,BlockTotal))
for i in range (0, BlockTotal):
    Kroarray[0,i] = kro1

#Pressure Array
## Initialize Pressure Array
Pressure = np.zeros((Ntime,BlockTotal))
for i in range (0, BlockTotal):
   Pressure[0,i] = P1
   
#Initialize Saturation Array
SatWater = np.zeros((Ntime,BlockTotal))
for i in range (0, BlockTotal):
   SatWater[0,i] = S1
   
## Read J Producers Information
Jmatrix = np.zeros((BlockTotal,BlockTotal)) 
Jmatrix = np.asmatrix(Jmatrix)
    
## Read J Values from Excel file
J_df = input_xl.parse('J')
isEnd = False
i = 0
Jarray = np.zeros((BlockTotal,5))

Jcount = 0 
while not isEnd:
    pos = int(J_df.at[i,'Pos'])
        
    if pos < 0:
       break
    
    Jcount = Jcount + 1
    
    Jarray[pos,0] = float(J_df.at[i,'RW'])  ## Radius
    Jarray[pos,1] = float(J_df.at[i,'S'])   ## Skin
    Jarray[pos,2] = float(J_df.at[i,'PWF']) ## Constant Pressure Well
    Jarray[pos,3] = int(J_df.at[i,'Direction'])
    
    if Jarray[pos,3] == 1:
       c = 2 * math.pi * DZ[pos] * ((Kdata[pos,1] * Kdata[pos,0])**(1/2))
       b = UO * BO * (math.log(RE/Jarray[pos,0]) + Jarray[pos,1])
       Jmatrix[pos,pos] = c/b
       Jarray[pos,4] = Jmatrix[pos,pos]
       i = i + 1
       continue
   
    if Jarray[pos,3] == 2:
       c = 2 * math.pi * DX[pos] * ((Kdata[pos,2] * Kdata[pos,1])**(1/2))
       b = UO * BO * (math.log(RE/Jarray[pos,0]) + Jarray[pos,1])
       Jmatrix[pos,pos] = c/b
       Jarray[pos,4] = Jmatrix[pos,pos]
       i = i + 1
       continue
   
Jmatrix = Jmatrix * 6.33E-3

## Create the Jflow Matrix for JFlow Output file 

Jflow = np.zeros((Jcount,3))
i = 0 
while not isEnd:
    pos = int(J_df.at[i,'Pos'])
        
    if pos < 0:
       break
   
    ## pos    
    Jflow[i,0] = int(pos)
    Jflow[i,1] = Jmatrix[pos,pos]
    Jflow[i,2] = Jarray[pos,2]
    
    i = i + 1
   
## Read Injector Information

Qmatrix = np.zeros((BlockTotal,1)) 
Qmatrix = np.asmatrix(Qmatrix)

Qwmatrix = np.zeros((BlockTotal,1))
Qwmatrix = np.asmatrix(Qwmatrix)

## Read Injector Values from Excel file
Injectors_df = input_xl.parse('Injectors')
isEnd = False
i = 0
    
while not isEnd:
    pos = Injectors_df.at[i,'Pos']
        
    if pos < 0:
       break
       
    p1 = float(Injectors_df.at[i,'Q'])
    
    pos = int(pos)
    
    if pos < BlockTotal:
       if Qmatrix[pos] != 0:
          print("Position: ", pos, " has already been assigned to another injector. /n Check your input.")
          exit()
       else:   
          Qmatrix[pos] = p1
          Qwmatrix[pos] = p1   
        
    i = i + 1

## Read Producer Information
Producers_df = input_xl.parse('Producers')
isEnd = False
i = 0
    
while not isEnd:
    pos = Producers_df.at[i,'Pos']
        
    if pos < 0:
       break
       
    p1 = float(Producers_df.at[i,'Qu'])
    
    pos = int(pos)
    if pos < BlockTotal:
        
       if Qmatrix[pos] != 0:
          print("Position: ", pos, " has already been assigned . /n Check your input.")
          exit()
       else:   
          Qmatrix[pos] = p1
    
        
    i = i + 1

#Add the J to the QMatrix    
for i in range(0,BlockTotal):
    if Jmatrix[i,i] != 0:
       Qmatrix[i] = Qmatrix[i] + Jmatrix[i,i] * Jarray[i,2]
       
       
# Close Input file
input_xl.close()       

## =================== End Create Global Variable Conditions ============= 

    
##========================= FUNCTION SECTION ==========================## 
       
def HarmonicMean(k1,k2):

    k = 2*(1/(1/k1 + 1/k2))

    return k      

def upstream(krarray,parray,i,j,IsZdirection,isYdirection,u):
## Input Variables
    
## karray: relative permeability
## parray: gridblock pressures
## i: gridblock position 1 to compare
## j: gridblock position 2 to compare
## IsZdiretion: are we comparing Z-direction
## IsYdirection: are we comparing Ydirection
## X-direcion is default
## u: the u of the phase in consideration
    
## Output Variable
#t: transmisibility
    
    
    
   if IsZdirection:
      k = HarmonicMean(Kdata[i,2], Kdata[j,2])
      
      if parray[i]>=parray[j]:      
          t = (k*DY*DX/(u * BO * DZ[i])) * krarray[i]
      else: 
          t = (k*DY*DX/(u * BO * DZ[j])) * krarray[j]
          
      return t    
          
   if isYdirection:       
      
      k = HarmonicMean(Kdata[i,1], Kdata[j,1])   
       
      if parray[i]>=parray[j]:
          t = (k*DX*DZ[i]/(u * BO * DY)) * krarray[i]
      else: 
          t = (k*DX*DZ[j]/(u * BO * DY)) * krarray[j]
          
      return t
    
   ## Is X-Direction - Default
   
   k = HarmonicMean(Kdata[i,0], Kdata[j,0])
   
   if parray[i]>=parray[j]:
       t = (k*DY*DZ[i]/(u * BO * DX)) * krarray[i]
   else: 
       t = (k*DY*DZ[j]/(u * BO * DX)) * krarray[j]     
          
   return t

#=================== END upstream =================#
   
def Compute_Matrix_1D(tmatrix,karray,parray,start,isTop,isBottom,isFront,isBack,u):
## Input Variables
    
## karray: relative permeability
## parray: gridblock pressures 
## start: starting position in gridblock
## isTop, isBottom, isFront, isBack : boundary conditions
## u: the u of the phase in consideration
    
## Output Variables: 
## tmatrix : final matrix being built

 
    for j in range (0,XBLOCKS):
        
        i = j + start
        
        if (j%XBLOCKS == 0):
        ## First Block in the row - XDirection
        
            ## Left - xDirection
            t1 = 0
           
            ## Back - yDirection
            if isBack:
               t2 = 0
            else:   
               t2 = upstream(karray,parray,i,i+XBLOCKS,False,True,u)
               tmatrix[i,i+XBLOCKS] = -t2
           
            #Right - xDirection
            t3 = upstream(karray,parray,i,i+1,False,False,u)
            tmatrix[i,i+1] = -t3
           
            #Front - yDirection
            if isFront:
               t4 = 0
            else:
               t4 = upstream(karray,parray,i,i-XBLOCKS,False,True,u)
               tmatrix[i,i-XBLOCKS] = -t4
            
             # TOP - zDirection
            if isTop:
               t5 = 0
            else:   
               t5 = upstream(karray,parray,i,i-LBLOCKS,True,False,u)
               tmatrix[i,i-LBLOCKS] = -t5
           
            # BACK - Vertical
            if isBottom:
               t6 = 0
            else:
               t6 = upstream(karray,parray,i,i+LBLOCKS,True,False,u)
               tmatrix[i,i+LBLOCKS] = -t6
           
            ## Set Diagonal
            tmatrix[i,i] = t1 + t2 + t3 + t4 + t5 + t6
            
            continue
        
        if j < XBLOCKS and ((j+1)%XBLOCKS != 0):
           ## Side Block
           
            ## Left - xDirection
            t1 = upstream(karray,parray,i,i-1,False,False,u)
            tmatrix[i,i-1] = -t1
            
            ## Back - yDirection
            if isBack:
               t2 = 0
            else:     
               t2 = upstream(karray,parray,i,i+XBLOCKS,False,True,u)
               tmatrix[i,i+XBLOCKS] = -t2
            
            ## Right
            t3 = upstream(karray,parray,i,i+1,False,False,u)
            tmatrix[i,i+1] = -t3
            
            ## Front - yDirection
            if isFront:
               t4 = 0
            else:   
               t4 = upstream(karray,parray,i,i-XBLOCKS,False,True,u)
               tmatrix[i,i-XBLOCKS] = -t4
               
            # Top - zDirection
            if isTop:
               t5 = 0
            else:  
               t5 = upstream(karray,parray,i,i-LBLOCKS,True,False,u)
               tmatrix[i,i-LBLOCKS] = -t5
           
            # Bottom - zDirection
            if isBottom:
              t6 = 0
            else:
              t6 = upstream(karray,parray,i,i+LBLOCKS,True,False,u)
              tmatrix[i,i+LBLOCKS] = -t6
             
            #Set Diagonal
            tmatrix[i,i] = t1 + t2 + t3 + t4 + t5 + t6 
            
            continue
        
        if (j+1)%XBLOCKS == 0:
        ## Last Block in the row - xDirection
            
            ## Left - xDirection
            t1 = upstream(karray,parray,i,i-1,False,False,u)
            tmatrix[i,i-1] = -t1
            
            ## Back - yDirection
            if isBack:
               t2 = 0
            else:   
               t2 = upstream(karray,parray,i,i+XBLOCKS,False,True,u)
               tmatrix[i,i+XBLOCKS] = -t2
            
            ## Right - xDirection
            t3 = 0
            
            ## Front - yDirection
            if isFront:
               t4 = 0
            else:
               t4 = upstream(karray,parray,i,i-XBLOCKS,False,True,u)
               tmatrix[i,i-XBLOCKS] = -t4
            
            # Top - zDirection
            if isTop:
               t5 = 0
            else:   
               t5 = upstream(karray,parray,i,i-LBLOCKS,True,False,u)
               tmatrix[i,i-LBLOCKS] = -t5
           
            # Bottom - zDirection
            if isBottom:
              t6 = 0
            else:
              t6 = upstream(karray,parray,i,i+LBLOCKS,True,False,u)
              tmatrix[i,i+LBLOCKS] = -t6
             
            #Set Diagonal 
            tmatrix[i,i] = t1 + t2 + t3 + t4 + t5 + t6
            
            continue
           
    return

#=================== END Compute_Matrix_1D =================#  

def Compute_Matrix_2D(tmatrix,karray,parray,start,isTop,isBottom,u):
## Input Variables
    
## karray: relative permeability
## parray: gridblock pressures 
## start: starting position in gridblock
## isTop, isBottom : boundary conditions
## u: the u of the phase in consideration
    
## Output Variables: 
## tmatrix : final matrix being built
    
    for j in range(0,YBLOCKS):
        
        ## First layer is Front layer
        if j == 0:
          isFront = True
        else:
          isFront = False
          
        ## Last layer is bottom layer  
        if j == YBLOCKS - 1:
          isBack = True
        else:
          isBack = False
          
        ## Set starting point for next layer
        xstart = start + j * XBLOCKS
        
        Compute_Matrix_1D(tmatrix,karray,parray,xstart,isTop,isBottom,isFront,isBack,u)
        
    return
        
#=================== END Compute_Matrix_2D =================# 
        
def Compute_Matrix_3D(karray,parray,u):
## Input Variables
    
## karray: relative permeability
## parray: gridblock pressures 
## u: the u of the phase in consideration
    
## Output Variables: 
## tmatrix : final matrix being built
    
    ## Is this the bottom layer?
    isBottom = False
    
    ## Is this the top layer?
    isTop = False
    
    ## Initialize tmatrix
    tmatrix = np.zeros((BlockTotal,BlockTotal))
    
    for i in range(0,ZBLOCKS):
        
        ## First layer is top layer
        if i == 0:
          isTop = True
        else:
          isTop = False
          
        ## Last layer is bottom layer  
        if i == ZBLOCKS - 1:
          isBottom = True
        else:
          isBottom = False
          
        ## Set starting point for next layer
        ystart = LBLOCKS * i
        
        Compute_Matrix_2D(tmatrix,karray,parray,ystart,isTop,isBottom,u)
        
    return tmatrix    
       
#=================== END Compute_Matrix_3D =================#   
        
def CalBMatrix(vmatrix):
## Input Variables
## vmatrix: Volumes for each gridblock
    
## Output Variables: 
## tmatrix : final Bmatrix
   
   ## Create B Matrix
   tmatrix = np.zeros((BlockTotal,BlockTotal))
   for i in range (0, BlockTotal):
       
       b = CT * PHI * vmatrix[i]  #ft^3/psi
       tmatrix[i,i] = b

   tmatrix = np.asmatrix(tmatrix)
   
   return tmatrix 
## ========== end ==========================
   
def InitPressureMatrix(Pi):
## Input Variables
## Pi
    
## Output Variables: 
## tmatrix : final Pressure Matrix
    
   tmatrix = np.zeros((BlockTotal,1))
   
   for i in range (0, BlockTotal):
       tmatrix[i] = Pi[i]
   
   tmatrix = np.asmatrix(tmatrix)
   
   return tmatrix
## ========== end ======================
   
def Compute_Saturation(swmatrix,bw,dt,volmatrix,twmatrix,pmatrix,qwmatrix):
## Solve explicitly for saturation
    
## Input Variables
## swmatrix: water saturation matrix
## bw property
## dt: delta time
## volmatrix: each gridblock volume
## twmatrix : water transmisibility matrix
## pmatrix : pressure(n+1) of each gridblock
## qwmatrix: Injectors matrix
    
## Output Variables: 
## tmatrix : final Saturation Matrix    
 
    
    s = np.zeros((BlockTotal,1))
    s = np.asmatrix(s)
    for i in range(0,BlockTotal):
       s[i,0] = swmatrix[i]
 
    #tmatrix = s + ((bw * dt)/(vol * PHI))*((-twmatrix * 6.33E-3 * pmatrix) + qwmatrix)
    
    a = ((bw * dt)/(volmatrix * PHI))
    b = ((-twmatrix * 6.33E-3 * pmatrix) + qwmatrix)
    
    for i in range(0,BlockTotal):
        b[i] = a[i] * b[i]
    
    tmatrix = s + b
    
    tmatrix = np.asmatrix(tmatrix)
    return tmatrix
## =============== end ===============

## == Calculate K for Water    
def CalKWaterBrookCorey(swater):
## Input Variables
## swater: Water Saturation for each gridblock
    
## Output Variables: 
## karray : relative water permeability for each gridblock    
    
    ##This is an array
    krarray = np.zeros((BlockTotal))
    for i in range(0,BlockTotal):
      if swater[i] <= .8:
         s = (swater[i] - SWI)/(1-SWI-SWR)
         krarray[i] = .2 * (s ** 3)
      else:
         krarray[i] = .2      
    
    return krarray
## ==============  end ================
    
def CalKOilBrookCorey(swater):
## Input Variables
## swater: Water Saturation for each gridblock
    
## Output Variables: 
## karray : relative oil permeability for each gridblock     
    
    ##This is an array
    krarray = np.zeros((BlockTotal))
    
    for i in range(0,BlockTotal):
        
       sw = (swater[i] - SWI)/(1-SWI-SWR)
       krarray[i] = (1 - sw) ** 3  
    
    return krarray
## ==============  end ================ 
    
def CalKGas(sgas):
## SPE relative permeability Model 1
## Input Variables
## sgas: Gas Saturation for each gridblock
    
## Output Variables: 
## karray : relative permeability for each gridblock 
    
    ##This is an array
    sarray = np.zeros((BlockTotal))
    
    for i in range(0,BlockTotal):
     
      sg = -2E-7*sgas[i]**4 - 1E-5*sgas[i]**3 + .0025*sgas[i]**2 - .0945*sgas[i] + 1.0915  
      sarray[i] = sg  
       
    return sarray
## ==============  end ================

def CalKOil(soil):
## SPE relative permeability Model 1
## Input Variables
## sgas: Gas Saturation for each gridblock
    
## Output Variables: 
## karray : relative permeability for each gridblock 
    
    
    ##This is an array
    sarray = np.zeros((BlockTotal))
    
    for i in range(0,BlockTotal):
      so = 3E-8*soil[i]**5 - 4E-7*soil[i]**4 + 3E-6*soil[i]**3 - 8E-6*soil[i]**2 +1E-5*soil[i] - 6E-6    
      #so = 2E-6 * (soil[i])**4 - 8E-5 * (soil[i])**3 + .0012 * (soil[i])**2 - 0.0065 * soil[i] + .0093  
      sarray[i] = so  
    
    return sarray
## ==============  end ================
    
def SaveMatrix(fhandle,filename, ma, isStart):
## Save current data. Re-opens a new file handle when threshold is reached
    
    mdata = np.squeeze(np.asarray(ma))
    
    if isStart:
        
       fhandle.close()
       fhandle = open(filename,"w+")
       
       ## Write Header
       tstring = ""
       for i in range(0,BlockTotal):
           tstring = tstring + str(i) + ","
       tstring = tstring + "\r\n"
       fhandle.write(tstring)
    
    
    tstring = ""
    for i in range(0,BlockTotal):
        tstring = tstring + str(mdata[i]) + ","
    tstring = tstring + "\r\n"
    fhandle.write(tstring)
    
    return fhandle
        
## ====================== END FUNCTIONS =======================
  
##==================== START MAIN PROGRAM =====================
    
## Create Gravity Matrix
Gmatrix = np.zeros((BlockTotal,2)) 
Gmatrix = np.asmatrix(Gmatrix) 

# Create Bmatrix
Bmatrix = CalBMatrix(Vmatrix)  

## Open OutPut CSV Text Handles

## (1) PRESSURE
fpressure = open("Pressure-Output-0.csv","w+")
## Write Header
tstring = ""
for i in range(0,BlockTotal):
    tstring = tstring + str(i) + ","
tstring = tstring + "\r\n"
fpressure.write(tstring)

## (2) Saturation
fsaturation = open("Saturation-Output-0.csv","w+")
## Write Header
tstring = ""
for i in range(0,BlockTotal):
    tstring = tstring + str(i) + ","
tstring = tstring + "\r\n"
fsaturation.write(tstring)

## (3)  JFlow
## Open File and Create Headers
## Write Header

fJflow = open("JFlow-Output.csv","w+")
tstring = ""
for i in range(0,Jcount):
    tstring = tstring + "J - " + str(int(Jflow[i,0,])) + ","
tstring = tstring + "\r\n"
fJflow.write(tstring)
FileSizeCount = 0

for n in range(0,Ntime-1):
    
    if n%10 == 0:
       ## Let us knkow how far we are with calculations 
       print(n)  
    
    #Create Transmisibility of Water matrix
    Twmatrix = Compute_Matrix_3D(Krwarray[n], Pressure[n], UW)
    
    #Create Trasmisibility of Oil Matrix
    Tomatrix = Compute_Matrix_3D(Kroarray[n], Pressure[n], UO)
    
    # Add water and oil matrixes
    Tmatrix  = (Twmatrix + Tomatrix) * 6.33E-3
    
    # If initial conditions asked for, compute gravity
    if isGravity == 1: 
       Gmatrix = .433 * Tmatrix * Zfmatrix
       Gmatrix = np.asmatrix(Gmatrix)
    
    # Create matrix with Initial Pressures for each iteration
    Pmatrix = InitPressureMatrix(Pressure[n])
    
    # QMATRIX and Jmatrix are Created when setting Initial Conditions
    
    a = Tmatrix + (Bmatrix / DT) + Jmatrix
    b = ((Bmatrix / DT) * Pmatrix) + Qmatrix #+ Gmatrix
    
    Pn1pmatrix = np.linalg.solve(a, b)
    
    #(1) Save Pressures
    Pn1parray = np.squeeze(np.asarray(Pn1pmatrix))
    for i in range(0,BlockTotal):
        Pressure[n+1,i] = Pn1parray[i]
    
    #(2) Save Pressures    
    if n%500 == 0 and n!=0:  
       fpressure = SaveMatrix(fpressure,"Pressure-Output-" + str(n) + ".csv", Pn1pmatrix, True) 
    else:   
       SaveMatrix(fpressure,"Pressure-Output-" + str(n) + ".csv", Pn1pmatrix, False)   
       
    #Save Q JFlow
    
    tstring = ""
    Psum = 0
    Plen = len(Pn1pmatrix)
    for i in range(0,Plen):
        Psum = Psum + np.squeeze(np.asarray(Pn1pmatrix[i]))
    Pavg = Psum/Plen    
    for i in range(0,Jcount):
        QJ = float(Jflow[i,1]) * 6.33E-3 * (Pavg - float(Jflow[i,2]))
        tstring = tstring + str(QJ) + ","
    tstring = tstring + "\r\n"
    fJflow.write(tstring)    
        
    Satnextmatrix = Compute_Saturation(SatWater[n],BW,DT,Vmatrix,Twmatrix,Pn1pmatrix,Qwmatrix)
    
    #(1) Save Saturations
    Satnextarray = np.squeeze(np.asarray(Satnextmatrix))
    for i in range(0,BlockTotal):
        SatWater[n+1,i]= Satnextarray[i]
     
    #(2) Save Saturations    
    if n%500 == 0 and n!=0:  
       fsaturation = SaveMatrix(fsaturation,"Saturation-Output-" + str(n) + ".csv", Satnextmatrix, True) 
    else:   
       SaveMatrix(fsaturation,"Saturation-Output-" + str(n) + ".csv", Satnextmatrix, False)   
    
    #Calculate KWater Permeability for next iteration
    Krwarray[n+1] = CalKWaterBrookCorey(SatWater[n+1])
    #TKoutput_df = pd.DataFrame(Krwarray, index=range(0,Ntime), columns=range(0,BlockTotal))
    
    #Calcaculate next oil permeability for next iteration
    Kroarray[n+1] = CalKOilBrookCorey(SatWater[n+1])
    
## Close Output CSV files
fpressure.close()
fsaturation.close()
fJflow.close()