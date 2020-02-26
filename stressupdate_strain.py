#


import numpy as np 
import matplotlib.pyplot as plt
#**********************
#C Stop program if neither plane strain nor axisymmetric state



def sy(ep):
    a = 200
    b = 100
    c = .9
    syield = a + b *ep#*(1-np.exp(-c*ep))
    dsyield = b#*c*np.exp(-c*ep)
    return syield, dsyield

R1, R2 , R3 =1.0, 2.0, 3.0
#C Initialise some algorithmic and internal variables
#C Set some material properties
YOUNG=300*1000
POISS=.3
#NHARD=I
#C Shear and bulk moduli and other necessary constants
GMODU=YOUNG/(R2*(R1+POISS))
BULK=YOUNG/(R3*(R1-R2*POISS))
R2G=R2*GMODU
R3G=R3*GMODU
STRAT = np.zeros(4) 
RSTAVA = np.zeros(4)
PSTAVA = np.zeros(4)
EET =   np.zeros(4)
STRES = np.zeros(4)
epn = 0
#C Elastic predictor: Compute elastic trial state
#C ----------------------------------------------
#C Volumetric strain and pressure stress
sarray = np.array([])
earray = np.array([])
E11 = 0
DGAMA = 0
for i in range(100):
    DELE = 0.00001
    E11 = E11 +DELE
    STRAT[0] = STRAT[0] + DELE #DELSTR[0]
    STRAT[1] = STRAT[1] #- POISS*DELE#DELSTR[1]
    STRAT[2] = STRAT[2] + 0#DELSTR[2]
    STRAT[3] = STRAT[3] #DELSTR[3]
    
    EEV=STRAT[0]+STRAT[1]+STRAT[3]
    P=BULK*EEV
    EEVD3=EEV/R3
    EET[0]=STRAT[0]-EEVD3
    EET[1]=STRAT[1]-EEVD3
    EET[3]=STRAT[3]-EEVD3
    
    EET[2]=STRAT[2]/R2
    
    VARJ2T=R2G*R2G*(EET[2]**2+0.5*(EET[0]**2+\
                        EET[1]**2+EET[3]**2))
    QTRIAL=(R3*VARJ2T)**(1/2)
    SIGMAY= 200#sy(epn)[0]
    TOL = 1*10**-6
    #C -------------------------------
    PHI=QTRIAL-SIGMAY
    if PHI/SIGMAY <= TOL:
        STRES[0]=R2G*EET[0]+P
        STRES[1]=R2G*EET[1]+P
        STRES[2]=R2G*EET[2]
        STRES[3]=R2G*EET[3]+P
        
        RSTAVA[0]=STRAT[0]
        RSTAVA[1]=STRAT[1]
        RSTAVA[2]=STRAT[2]
        RSTAVA[3]=STRAT[3]
        


    else:
        SIGMAY, DSIGMAY = 200, 100
        DGAMA = 0
        DENOM = -R3G-DSIGMAY
        epnn = epn
        #for i in range(100):
            #C Compute residual derivative
            
            #C Compute Newton-Raphson increment and update variable DGAMA
        DGAMA = -PHI/DENOM
        #DGAMA=DGAMA+DDGAMA
            #C Compute new residual
        epnn=epnn+DGAMA
            #SIGMAY, DSIGMAY = sy(epnn)
          #  PHI=QTRIAL-R3G*DGAMA-SIGMAY
           # DENOM = -R3G-DSIGMAY
            #C Check convergence
           # RESNOR=abs(PHI/SIGMAY)
           # if (RESNOR<=TOL):
           #     break
            #else:
             #   print('convergence failed')
              #  print(RESNOR)
            
            #C update accumulated plastic strain
        #epn = EPBAR 
#C update stress components
        epn = epnn
        FACTOR=R2G*(R1-R3G*DGAMA/QTRIAL)
        STRES[0]=FACTOR*EET[0]+P
        STRES[1]=FACTOR*EET[1]+P
        STRES[2]=FACTOR*EET[2]
        STRES[3]=FACTOR*EET[3]+P
        #C compute converged elastic (engineering) strain components
        FACTOR1=FACTOR/R2G
        RSTAVA[0]=FACTOR1*EET[0]+EEVD3
        RSTAVA[1]=FACTOR1*EET[1]+EEVD3
        RSTAVA[2]=FACTOR1*EET[2]*R2
        RSTAVA[3]=FACTOR1*EET[3]+EEVD3   
        ed1, ed2, ed3, ed4 = FACTOR*EET[0],FACTOR*EET[1],FACTOR*EET[2],FACTOR*EET[3]
        ednorm = (ed1**2 + ed2**2 + ed4**2 + 2*ed3**2)**(1/2)
        PSTAVA[0]= PSTAVA[0] + (3/2.0)**(1/2)*DGAMA*ed1/ednorm
        PSTAVA[1]= PSTAVA[1] + (3/2.0)**(1/2)*DGAMA*ed2/ednorm
        PSTAVA[2]= PSTAVA[2] + (3/2.0)**(1/2)*DGAMA*ed3/ednorm
        PSTAVA[3]=  PSTAVA[3] +(3/2.0)**(1/2)*DGAMA*ed4/ednorm
    START = RSTAVA
    sarray = np.append(sarray , STRES[0])
    earray = np.append(earray , E11)
    plt.plot(earray, sarray)
print('RSTAVA', RSTAVA)
print('PSTAVA', PSTAVA)
print('E11', E11)
print(' SUM', RSTAVA + PSTAVA)
    
        
        
        
        
        