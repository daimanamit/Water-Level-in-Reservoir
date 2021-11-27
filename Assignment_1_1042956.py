# -*- coding: utf-8 -*-
"""

@author: ada003
"""

import numpy as np
import matplotlib.pyplot as plt
#----------------------------------------------------------------------------------------
# Parameters
#----------------------------------------------------------------------------------------
b = 90.0        # Spillway width
h0 = 40.0       # initiall water elevation
a = 40.5        # Spillway elevation
m = 0.5         # Spillway constan
Q = 7           # Flow of the bottom outlet m3/s
dt = 100      # time step in seconds
T = 100*3600      # Total time for computation in hours

#-------------------------------------------------------------------------------------------
#Intensity in function of time
#--------------------------------------------------------------------------------------------
I =     [11, 11, 40, 76, 100, 85, 60, 35, 25, 20, 10, 0] # m3 / s
tI =    [0, 18000, 36000, 54000, 72000, 90000, 108000, 126000, 144000, 162000, 180000, 360000] # Second

#-------------------------------------------------------------------------------------------
#Reservoir Area in function of water level
#-------------------------------------------------------------------------------------------
Ar =    [0, 1.3, 5.4]  # Area in Km2
Ar_h =  [20, 40, 50]  # Elevation of the reservoir in masl

#-------------------------------------------------------------------------------------------
# Number of steps
#-------------------------------------------------------------------------------------------
Nt = int(T/dt)+1

#Nt[-1] = T*60
tmat=list(np.zeros(Nt))
for i1 in range(Nt):
   tmat[i1]=(i1*dt)
#-------------------------------------------------------------------------------------------
#Water level arrays
#-------------------------------------------------------------------------------------------
H_Ex = []         # list values of h Explicit scheme
H_Im = []          # list values of h Implicit scheme

#------------------------------------------------------------------------------------------
# Over flow form mthe weir
#------------------------------------------------------------------------------------------
Qw_Ex = [] # Flow over the weir in m3/s explicit
Qw_Im = [] # Flow over the weir in m3/s implicit
Volume_W = [ ] # Total water volume over the weir

#-----------------------------------------------------------------------------------------
# initial values
#-----------------------------------------------------------------------------------------
H_Ex.append(h0)
H_Im.append(h0)
Qw_Ex.append(0)
Qw_Im.append(0)
Volume_W.append(((Q)*(dt)))

Error = 0.001

#-----------------------------------------------------------------------------------------
# Solution for Explicit scheme
#-----------------------------------------------------------------------------------------
for i in range(1,Nt):
    if H_Ex[i-1] > a:
        
        h_Temp = H_Ex[i-1] + ((dt*(np.interp((i-1)*dt,tI,I) - Q - m*b*(H_Ex[i-1] - a)**1.5)) /(np.interp(H_Ex[i-1] , Ar_h , Ar) * 1000000))
        Qw_Ex.append(m*b*(H_Ex[i-1] - a)**1.5)
        Volume_Temp = (abs(Qw_Ex[i])*(dt)) + ((Q)*(dt))
        Volume_W.append(Volume_Temp)
        
    else:
       h_Temp = H_Ex[i-1] + ((dt*(np.interp((i-1)*dt,tI,I)-Q))/(np.interp(H_Ex[i-1] , Ar_h , Ar) * 1000000))
       Qw_Ex.append(0)
       Volume_W.append(((Q)*(dt)))
       
    H_Ex.append(h_Temp)

#-------------------------------------------------------------------------------------------
# Solution for Implicit scheme
#-------------------------------------------------------------------------------------------    
for i in range(1,Nt):
    E = 1
    H_Iter =  H_Im [i-1]
    while E > Error: 
        if H_Iter > a:        
            h_Temp = H_Im [i-1] + ((dt*(np.interp(i*dt,tI,I) - Q - m*b*(H_Iter - a)**1.5)) /(np.interp(H_Iter , Ar_h , Ar) * 1000000))
            Qw_Temp = (m*b*(H_Iter - a)**1.5)
        else:
           h_Temp = H_Im [i-1] + ((dt*(np.interp(i*dt,tI,I) - Q))/(np.interp(H_Iter , Ar_h , Ar) * 1000000))
           Qw_Temp = 0
           
        E = abs(abs(h_Temp) - abs(H_Iter))
         
        H_Iter = h_Temp
     
    Qw_Im.append(Qw_Temp)       
    H_Im.append(h_Temp)

#----------------------------------------------------------------------------------------------
# plotting
#---------------------------------------------------------------------------------------------- 
# Plotting Outflowing Volume of Water
plt.figure(figsize = [12,6])
plt.title('Total Outflowing Volume of Water from Reservoir', size = 20, color = 'g')
plt.ylabel('Volume of Water  (m^3)')
plt.xlabel('Time in Seconds')
plt.plot(tmat,Volume_W, linewidth = 0.8, label = 'Outflowing Volume')
plt.grid(color='lightgrey', linestyle=':', linewidth=1)
plt.legend(loc='best')
plt.show()

#Plotting Discharge from Weir
plt.figure(figsize = [12,6])
plt.title('Discharge from Weir', size = 20, color = 'g')
plt.ylabel('Discharge (m^3/s)')
plt.xlabel('Time in Seconds')
plt.plot(tmat,Qw_Im, linewidth = 0.8, label = 'Discharge (m^3/s)')
plt.grid(color='lightgrey', linestyle=':', linewidth=1)
plt.legend()
plt.show()

# Plotting Explicit Solution
plt.figure(figsize = [12,6])
plt.title('Explicit Solution', size = 20, color = 'g')
plt.ylabel('Water Level in masl')
plt.xlabel('Time in Seconds')
plt.plot(tmat, H_Ex, linewidth = 0.8, label = 'Explicit Solution')
plt.grid(color='lightgrey', linestyle=':', linewidth=1)
plt.legend()
plt.show()       

# Plotting Implicit Solution
plt.figure(figsize = [12,6])
plt.title('Implicit Solution', size = 20, color = 'g')
plt.ylabel('Water Level in masl')
plt.xlabel('Time in Seconds')
plt.plot(tmat, H_Im, linewidth = 0.8, label = 'Implicit Solution', color = 'g')
plt.grid(color='lightgrey', linestyle=':', linewidth=1)
plt.legend()
plt.show()       

# Plotting Emplicit and Implicit Solution
plt.figure(figsize = [12,6])
plt.title('Explicit and Implicit Solution, Δt = 100 Seconds', size = 20, color = 'g')
plt.ylabel('Water Level in masl')
plt.xlabel('Time in Seconds')
plt.plot(tmat, H_Ex, linewidth = 0.8, label = 'Explicit Solution')
plt.plot(tmat, H_Im, linewidth = 0.8, label = 'Implicit Solution')
plt.grid(color='lightgrey', linestyle=':', linewidth=1)
plt.legend() 
plt.show()     
    