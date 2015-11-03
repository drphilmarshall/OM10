from numpy import *
'''
  NAME 
     host
   
  PURPOSE
     to get properties of host galaxies, given redshift, host galaxy type and i band magnitude of QSO

  INPUT:
     z (redshift), mQi (i band magnitude of QSO), type ("e"=early type, "l"= late type)
  
  OUTPUT:
     magnitude of host in different bands:
     mhi, mhR, mhg, mhz, mhr

  AUTHORS
     Kai & Adri
     basing on Adri's relatinship among different bands.
   
  HISTORY
     2014-07-23  Kai & aanello

'''
####Input:
z=1.
mQi=20.
type="e"
#########

if z<0.5 or z>3.965:
   print "warning! out of redshift range."

f=open('$OM10_DIR/data/ETGcols.txt','r')
ETG=loadtxt(f)
f.close()
g=open('$OM10_DIR/data/LTGcols.txt','r')
LTG=loadtxt(g)
g.close()
k=open('$OM10_DIR/data/QSOcols.txt','r')
QSO=loadtxt(k)
k.close()

MRhQ=0.+1.*random.randn()
if type=="e":
   G=ETG
else:
   G=LTG


diff=list(abs(G[:,0]-z))
select = diff.index(min(diff))
mhi=mQi+G[select,1]-QSO[select,1]+MRhQ+G[select,5]
mhR=mhi-G[select,1]
mhg=mhi+G[select,2]+G[select,3]    
mhz=mhi-G[select,4]
mhr=mhi+G[select,3]
print "mhi, mhR, mhg, mhz, mhr"
print mhi, mhR, mhg, mhz, mhr