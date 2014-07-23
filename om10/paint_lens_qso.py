from numpy import *
import math
"""
    NAME
        paint_lens_qso

    PURPOSE
        to get the complete properties of lens and qso, given redshifts, velocity dispersion, i band magnitude of qso.

    
    INITIALISATION
        change the inputs at the beginning: z_om10_g (redshift of lens),
        sigma_om10 (velocity dispersion),z_om10_q (redshift of qso),i_q (i band magnitude).
        

    OUTPUT:
        case1: if it prints "out of region!", your inputs are out of
            z_om10_g:0-1, sigma_om10: 80-400, z_om10_q:0-5, i_q: 15-22

        case2:if it prints "non-zero", the outputs are directly from SDSS

        case3:if it prints "zero", there is no SDSS data aroud your inputs,
             outputs are obtained based on weighted values from all SDSS data. This is the same as case1.

        Lens Properties          
        redshift  sigma  g Reff  r   i    z   g mag   r mag   i mag   z   w1    w2    w3   w4
        QSO Properties          
        redshift   g     r     i     z     w1      w2       w3     w4  
   
    AUTHORS
      Kai Liao & aagnello

    COMMENTS
        this script is based on Adri's mathematica version, Kai converted it to Python and added some modification.  

    HISTORY
      2014-07-23  Kai & aanello
    """


######### Input:
z_om10_g=0.9
sigma_om10=60.
z_om10_q=4.5
i_q=14.5


## read data from SDSS
f=open('$OM10_DIR/data/LRGo.txt','r')
lrg=loadtxt(f)
f.close()
#print lrg[0,0],lrg.shape
g=open('$OM10_DIR/data/QSOo.txt','r')
qso=loadtxt(g)
g.close()
#print qso[0,0],qso.shape

######
zmin=0.
zmaxq=5.
zmaxg=1.
imin=15.
imax=22.
sming=80.
smaxg=400.
binsz=50
binsobs=20
step1=float((zmaxg-zmin)/binsz)
step2=float((smaxg-sming)/binsobs)
step3=float((zmaxq-zmin)/binsz)
step4=float((imax-imin)/binsobs)

RSDSS=1.5 #aperture size for SDSS 1.5 arcseconds
ng=zeros((500,200))
nq=zeros((500,200))
mlrg=zeros((50,20,14))
slrg=zeros((50,20,14))
mqso=zeros((50,20,9))
sqso=zeros((50,20,9))

##########

for i in range(lrg.shape[0]):
    
    sigma=lrg[i,1]*(2.45/2.0)**0.5*(lrg[i,3]/RSDSS)**-0.066
    kz=math.ceil((lrg[i,0]-zmin)/step1)-1
    ko=math.ceil((sigma-sming)/step2)-1
    
    if -1<kz<50 and -1<ko<20:
       ng[kz,ko]+=1
       mlrg[kz,ko,:]+=lrg[i,:]
       slrg[kz,ko,:]+=lrg[i,:]**2
       #slrg[kz,ko]+=


#print ng[16,9]
#meang=mlrg[16,9,:]/ng[16,9]
#print meang
#print sqrt(slrg[16,9,:]/ng[16,9]-meang**2)

###########

for j in range(qso.shape[0]):
    
    kz=math.ceil((qso[j,0]-zmin)/step3)-1
    ko=math.ceil((qso[j,3]-imin)/step4)-1
    
    if -1<kz<50 and -1<ko<20:
       nq[kz,ko]+=1
       mqso[kz,ko,:]+=qso[j,:]
       sqso[kz,ko,:]+=qso[j,:]**2


#print sum(nq)

#meanq=mqso[21,14,:]/nq[21,14]
#print meanq
#print sqrt(sqso[21,14,:]/nq[21,14]-meanq**2)



kz_g=math.ceil((z_om10_g-zmin)/step1)-1
ko_g=math.ceil((sigma_om10-sming)/step2)-1
if kz_g<50 and -1<ko_g<20 and ng[kz_g,ko_g]>0:
   
   print 'non-zero'
   meang=mlrg[kz_g,ko_g,:]/ng[kz_g,ko_g]
   errorg=sqrt(slrg[kz_g,ko_g,:]/ng[kz_g,ko_g]-meang**2)
   om10g=meang+(-1.+2.*random.rand(14))*errorg
else:
   if kz_g>49 or ko_g>19 or ko_g<0:
      print 'out of the region!'
   else:
      print 'zero'
   meangg=zeros((50,20,14))
   errorgg=zeros((50,20,14))
   disg=zeros((50,20))
   wmeang=zeros(14)
   werrorg=zeros(14)
   wg=0
   
   for i in range(50):
       for j in range(20):
           if ng[i,j]>0:
              meangg[i,j,:]=mlrg[i,j,:]/ng[i,j]
              errorgg[i,j,:]=sqrt(slrg[i,j,:]/ng[i,j]-meangg[i,j,:]**2)
           else:
              meangg[i,j,:]=mlrg[i,j,:]
              errorgg[i,j,:]= meangg[i,j,:]
           disg[i,j]=abs(z_om10_g-(i*step1+step1/2))+abs(sigma_om10-(j*step2+step2/2))
           wmeang += meangg[i,j,:]*e**(-disg[i,j])
           wg +=e**(-disg[i,j])
           werrorg += errorgg[i,j,:]*e**(-disg[i,j])
           
   om10g=wmeang/wg+(-1.+2.*random.rand(14))*werrorg/wg
print "Lens Properties"           
print "redshift  sigma  g Reff  r   i    z   g mag   r mag   i mag   z   w1    w2    w3   w4"   
print om10g


######################################3


kz_q=math.ceil((z_om10_q-zmin)/step3)-1
ko_q=math.ceil((i_q-imin)/step4)-1
if kz_q<50 and -1<ko_q<20 and ng[kz_q,ko_q]>0:

   print 'non-zero'
   meanq=mqso[kz_q,ko_q,:]/nq[kz_q,ko_q]
   errorq=sqrt(sqso[kz_q,ko_q,:]/nq[kz_q,ko_q]-meanq**2)
   om10q=meanq+(-1.+2.*random.rand(9))*errorq
else:
   if kz_q>49 or ko_q>19 or ko_q<0:
      print 'out of the region!'
   else:
      print 'zero'
   meanqq=zeros((50,20,9))
   errorqq=zeros((50,20,9))
   disq=zeros((50,20))
   wmeanq=zeros(9)
   werrorq=zeros(9)
   wq=0
   
   for i in range(50):
       for j in range(20):
           if nq[i,j]>0:
              meanqq[i,j,:]=mqso[i,j,:]/nq[i,j]
              errorqq[i,j,:]=sqrt(sqso[i,j,:]/nq[i,j]-meanqq[i,j,:]**2)
           else:
              meanqq[i,j,:]=mqso[i,j,:]
              errorqq[i,j,:]= meanqq[i,j,:]
           disq[i,j]=abs(z_om10_q-(i*step3+step3/2))+abs(i_q-(j*step4+step4/2))
           wmeanq += meanqq[i,j,:]*e**(-disq[i,j])
           wq +=e**(-disq[i,j])
           werrorq += errorqq[i,j,:]*e**(-disq[i,j])
           
   om10q=wmeanq/wq+(-1.+2.*random.rand(9))*werrorq/wq
print "QSO Properties"           
print "redshift   g     r     i     z     w1      w2       w3     w4"   
print om10q










