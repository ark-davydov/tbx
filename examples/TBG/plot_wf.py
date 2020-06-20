import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.interpolate import griddata
from matplotlib import rc

font_default=24
#rc('text', usetex=True)
#plt.rcParams.update({'font.size': font_default})
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')


nwan=2
nsym=12
rcut=300.0

def plot_frame(ax,xshift,yshift):
  frame_width=1.0
  col='k'
  if (yshift==0):
    col='r'
  x0,y0=0,0
  x1,y1=avec[0,0],avec[0,1]
  x0,x1=x0+xshift,x1+xshift
  y0,y1=y0+yshift,y1+yshift
  ax.plot([x0,x1],[y0,y1],c=col,lw=frame_width)
  x0,y0=0,0
  x1,y1=avec[1,0],avec[1,1]
  x0,x1=x0+xshift,x1+xshift
  y0,y1=y0+yshift,y1+yshift
  ax.plot([x0,x1],[y0,y1],c=col,lw=frame_width)
  avec12=avec[0,:]+avec[1,:]
  x0,y0=avec[0,0],avec[0,1]
  x1,y1=avec12[0],avec12[1]
  x0,x1=x0+xshift,x1+xshift
  y0,y1=y0+yshift,y1+yshift
  ax.plot([x0,x1],[y0,y1],c=col,lw=frame_width)
  x0,y0=avec[1,0],avec[1,1]
  x1,y1=avec12[0],avec12[1]
  x0,x1=x0+xshift,x1+xshift
  y0,y1=y0+yshift,y1+yshift
  ax.plot([x0,x1],[y0,y1],c=col,lw=frame_width)

def readwan(fname):
   avec=np.zeros((2,2))
   f=open(fname,'r')
   a,b,c=f.readline().strip().split()
   avec[0,:]=np.array([float(a),float(b)])
   a,b,c=f.readline().strip().split()
   avec[1,:]=np.array([float(a),float(b)])
   a,b,c=f.readline().strip().split()
   at,nats=f.readline().strip().split()
   nr,nrwan=f.readline().strip().split()
   natoms=int(nats)
   natoms_per_layer=int(float(nats)/2)
   nrwan=int(nrwan)
   atc1=[]
   atc2=[]
   wan1=[]
   wan2=[]
   for ir in range(nrwan):
     x,y,z=f.readline().strip().split()
     vcur=np.matmul(np.array([float(x),float(y)]),avec)
     dcurfix=np.matmul(np.array([0.5,0.5]),avec)
     for i in range(natoms):
       a,b,c,d,e=f.readline().strip().split()
       z=float(c)
       atcur=[float(a)+vcur[0],float(b)+vcur[1]]
       dxy=np.linalg.norm(np.array(atcur)-dcurfix)
       if dxy<rcut:
         if (z<0):
           atc1.append(atcur)
           wan1.append(np.complex(float(d),float(e)))
         else:
           atc2.append(atcur)
           wan2.append(np.complex(float(d),float(e)))
   return avec,np.array(atc1).transpose(),np.array(atc2).transpose(),\
             np.abs(np.array(wan1)),np.abs(np.array(wan2))

def plotcmap(x,y,z,ninterp,ax):
  xy=np.vstack((x,y)).transpose()
  xyz=np.vstack((x,y,z)).transpose()
  data=pd.DataFrame(xyz,columns = list("XYZ"))
  numcols, numrows = ninterp, ninterp
  xi = np.linspace(data.X.min(), data.X.max(), numcols)
  yi = np.linspace(data.Y.min(), data.Y.max(), numrows)
  xi, yi = np.meshgrid(xi, yi)
  zi = griddata(xy, z, (xi, yi),method='linear')
  cs=ax.contourf(xi, yi, zi, vmin=0,vmax=data.Z.max())
  return cs

def plot_figure(avec,atc1,atc2,wan1,wan2,prefix,wfname):
    sz=0.01
    fig, ax = plt.subplots(num=None, figsize=(20,20), dpi=100, facecolor='w', edgecolor='white')
    ax.set_aspect(aspect=1)
    ax.axis('off')
    wmin=min(np.hstack((wan1,wan2)))
    wmax=max(np.hstack((wan1,wan2)))
    cs=plotcmap(atc1[0,:]-xdist,atc1[1,:],wan1,200,ax)
    cs=plotcmap(atc2[0,:]+xdist,atc2[1,:],wan2,200,ax)
    cbaxes = fig.add_axes([0.20, 0.25, 0.6, 0.02])
    cb=plt.colorbar(cs,orientation='horizontal',shrink=0.6,cax=cbaxes)
    line='$p_z$ component weight of '+wfname
    print(prefix)
    print(line)
    
    cb.set_label(line,size=40,labelpad=20)
    cb.ax.tick_params(
      axis='x',          # changes apply to the x-axis
      which='both',      # both major and minor ticks are affected
      bottom=False,      # ticks along the bottom edge are off
      top=True,         # ticks along the top edge are off
      labelbottom=False, # labels along the bottom edge are off
      labeltop=True) # labels along the bottom edge are off

    plot_frame(ax,-xdist,0)
    plot_frame(ax,+xdist,0)
    mxy=np.amax(atc1[1,:])
    dxfix=0.2*xdist
    ax.text(-xdist-dxfix,mxy+0.1*mxy,r'layer 1',fontsize=40)
    ax.text( xdist-dxfix,mxy+0.1*mxy,r'layer 2',fontsize=40)
    plot_centers(ax,xdist)
    plt.savefig(prefix+'.png',bbox_inches='tight',dpi=200)
    plt.cla()
    plt.clf()
    plt.close()

def plot_centers(ax,xdist):
   centers=[]
   with open('wannier_centres.xyz') as f:
      ncenters=float(f.readline().strip())
      f.readline()
      for i in range(int(ncenters/2)):
         f.readline()
      for i in range(int(ncenters/2)):
         _,x,y,z=f.readline().strip().split()
         centers.append(np.array([float(x),float(y),float(z)]))
   for ce in centers:
      dcurfix=np.matmul(np.array([0.5,0.5]),avec)
      for i in range(-4,4):
         for j in range(-4,4):
            vc=np.matmul(np.array([i,j]),avec)
            x,y=ce[0]+vc[0],ce[1]+vc[1]
            if (np.linalg.norm(np.array([x,y])-dcurfix)<rcut):
              ax.scatter(x-xdist,y,s=200,c='gray',zorder=100,alpha=0.7)
            x,y=ce[0]+vc[0],ce[1]+vc[1]
            if (np.linalg.norm(np.array([x,y])-dcurfix)<rcut):
              ax.scatter(x+xdist,y,s=200,c='gray',zorder=100,alpha=0.7)
            
  
   

xdist=360
nwan=12
iwan=1
for i in range(iwan,nwan+1):
  prefix='wfmloc'+str(i)
  wfname=r'$w_{'+str(i)+'}$'
  fname=prefix+'.dat'
  if os.path.exists(fname):
    avec,atc1,atc2,wan1,wan2=readwan(fname)
    plot_figure(avec,atc1,atc2,wan1,wan2,prefix,wfname)

for i in range(iwan,nwan+1):
  prefix='wftrial'+str(i)
  wfname=r'$w_{'+str(i)+'}$'
  fname=prefix+'.dat'
  if os.path.exists(fname):
    avec,atc1,atc2,wan1,wan2=readwan(fname)
    plot_figure(avec,atc1,atc2,wan1,wan2,prefix,wfname)


for i in range(iwan,nwan+1):
  for isym in range(1,nsym+1):
    prefix='wfs'+str(i)+'_'+str(isym)
    wfname=r'$w_{'+str(i)+'}(g_{'+str(isym)+'})$'
    fname=prefix+'.dat'
    if os.path.exists(fname):
      avec,atc1,atc2,wan1,wan2=readwan(fname)
      plot_figure(avec,atc1,atc2,wan1,wan2,prefix,wfname)

for i in range(iwan,nwan+1):
  for isym in range(1,nsym+1):
    prefix='wfd'+str(i)+'_'+str(isym)
    fname=prefix+'.dat'
    wfname=r'$w_{'+str(i)+'}(g_{'+str(isym)+'})$'
    if os.path.exists(fname):
      avec,atc1,atc2,wan1,wan2=readwan(fname)
      plot_figure(avec,atc1,atc2,wan1,wan2,prefix,wfname)

