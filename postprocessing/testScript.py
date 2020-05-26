import matplotlib
#matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
import matplotlib.markers as mmarkers
import matplotlib.lines as mlines
import numpy as np
from numpy import inf
import sys
import math as m
import scipy.interpolate 
import os
from mpl_toolkits.axes_grid1 import AxesGrid
import codecs
import statistics
#import mpl_tookits.axisartist as AA
import subprocess as sp
import re


###############################################################################################
## Funktionen ################################################################################
###############################################################################################

def bubbleVel(U,V):
	return np.divide(np.sum(np.multiply(np.multiply(V,(1-fI)),U)),np.sum(np.multiply((1-fI),V)))

def bubbleCenter(x,y,z,v,fI):
    xb=np.divide(np.sum(np.multiply(np.multiply((1-fI),V),x)),np.sum(np.multiply((1-fI),V)))
    yb=np.divide(np.sum(np.multiply(np.multiply((1-fI),V),y)),np.sum(np.multiply((1-fI),V)))
    zb=np.divide(np.sum(np.multiply(np.multiply((1-fI),V),z)),np.sum(np.multiply((1-fI),V)))
    return xb,yb,zb

def bubbleVel(a,U):
    Ububble= np.divide(np.sum(np.multiply((1-a),U)),np.sum(1-a))
    return Ububble	

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
    
def include(filename):
    if os.path.exists(filename): 
        execfile(filename)



# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * iterate through case directory and find time folders * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
folder1='/home/Documents/Master/Masterarbeit/Pauls-4.0/run/'
folder2='water'
SimPath = os.join.path(folder1,folder2)
timei = []

for f in os.listdir(SimPath):
    if isFloat(f):
        timei.append(str(f))
time = sorted(timei)
tlen = len(time)


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * execute external application (writeCellCentres) -> only possible if openfoam is sourced in shell * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * -> can be done for all times in loop using time array  * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
if not os.path.exists(os.path.join(SimPath,'0','ccx')):
    
    cmdWriteCell = 'writeCellCentres -case ' + SimPath + ' -time 0'
    print('writing Cell Centers for time = 0')
    os.system(cmdWriteCell) 

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * read x,y files from case * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
meshSize0 = np.genfromtxt(os.path.join(SimPath,'0','ccx'),skip_header=20,max_rows=1)
meshSize=int(meshSize0)

x=np.genfromtxt(os.path.join(SimPath,'0','ccx'),skip_header=22,max_rows=meshSize)
y=np.genfromtxt(os.path.join(SimPath,'0','ccy'),skip_header=22,max_rows=meshSize)
#z=np.genfromtxt(os.path.join(SimPath,'0','ccz'),skip_header=22,max_rows=meshSize)
V=np.genfromtxt(os.path.join(SimPath,'0','V'),skip_header=22,max_rows=meshSize)
fluidIndicator=np.genfromtxt(os.path.join(SimPath,'0','fluidIndicator'),skip_header=22,max_rows=meshSize)

t=1

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * read fields (here alpha and U * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
FoamFile = os.path.join(SimPath,time[t])
FileName = os.path.join(FoamFile,'U')
with open(FileName, 'rb') as f:
    clean_lines = (line.replace(b'(',b' ').replace(b')',b' ') for line in f)
    U0 = np.genfromtxt(clean_lines, skip_header=22,max_rows=meshSize)

U=np.array(U0)



# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * calculate integral data * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *

Ububble = bubbleVel(U,V,fluidIndicator)
xb,yb,zb = bubblePos(x,y,z,V,fluidIndicator)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * prepare Plot data * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *

# if 3d
#idx = np.where(np.abs(z-zb)==np.abs(z-zb).min())
# Us= U[:,0]
# UX= Us[idx]
U_X = U[:,0]
U_Y = U[:,1]
U_y_rel = U[:,1]/Ububble
N = np.sqrt(U_X**2+U_Y**2)  # there may be a faster numpy "normalize" function
Uxn, Uyn = U_X/N, U_Y/N


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * contour plot example * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *

# create Mesh Grid
xiQ, yiQ = np.linspace(x.min(), x.max(), nMeshContourX), np.linspace(y.min(), y.max(), nMeshContourY)
xiQ, yiQ= np.meshgrid(xiQ, yiQ)
# interpolate velocity directions on mesh ( for vector plot )
qUx = scipy.interpolate.griddata((x, y), Uxn, (xiQ, yiQ), method='linear')
qUz = scipy.interpolate.griddata((x, y), Uzn, (xiQ, yiQ), method='linear')
# set ranges for color bar ( for contour plot)
Umin,Umax = -1,1
nMeshScalarPlot = 30
Uloc = np.linspace(Umin, Umax, nMeshScalarPlot)
UscalarLevels1 = FixedLocator(Uloc,nbins=nMeshScalarPlot).tick_values(Umin,Umax)

# create plot variable for contour plot
plotU = scipy.interpolate.griddata((x, y), U_y_rel, (xi, yi), method='linear')
# contour plot of relative velocity
Uplot1 = ax1.contourf(plotU,UscalarLevels1, cmap = 'seismic',extent=[x.min(), x.max(), y.min(), y.max()],extend='both')

# vector plot of normed velocity direction
QS = ax1.streamplot(xiQ,ziQ,qUx,qUz,color='k', linewidth=0.5, density=densityStream)

yMin =
yMax =
xMin =
xMax =
# domain size
ax1.set_ylim(yMin, yMax)
ax1.set_xlim(xMin, xMax)

ax0.axis('off')    
cbar = fig.colorbar(Uplot1, ax=ax0, extend='both',pad=-2,aspect=30,fraction=0.7)
cbar.set_ticks([-1,0,1])
cbar.set_ticklabels(['-1','0','1'])
cbar.set_label('U/U$_b$ / - ')
filename = SaveFolder  + filename + time[t] +'.png'
plt.savefig(filename, bbox_inches = 'tight',dpi=900)
plt.clf()

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *
# * * * * * parse through logFile * * * * *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * *

### on hlrn the logfile is named: logfile.oxxxxxxxx and xxxxxxx is a number that is increasing for each submitted job
### to append the files in the order they were created they need to be sorted

regNum = '([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)'
filenamesV = []
idxV = []
for f in os.listdir(SimPath):
    if 'logfile.o' in f:
        filenamesV.append(f)
        idxV.append(int(re.findall('logfile.o(\d+)',f)[0]))

reg = '^Time = ' + regNum + '\n'
progT=re.compile(reg)
reg = '^Bubble Velocity Vector = \(' + regNum +  ' ' + regNum + ' ' + regNum + '\)\n'
progUb=re.compile(reg)
