import numpy as np

def parameter():
    E = 200*1000
    nu = .3
    return E, nu
def TangentMatrix():
    E = parameter()[0]
    nu = parameter()[1]
    #D = E/(1-nu**2)*np.array([[1, nu, 0],[nu, 1, 0],[0, 0, (1-nu)/2]]) for plane stress
    D = E/((1-2*nu)*(1+nu))*np.array([[1-nu, nu, 0],[nu, 1-nu, 0],[0, 0, (1-2*nu)/2]])# for plane strain
    return D
def sigmay(ep):
    sy = 620 + 3300*(1- np.exp(-2*ep))
    dsy = 3300*2*np.exp(-2*ep)
    return sy, dsy

def stressupdate(straininc, stressold, epold, eqepold):
    E = parameter()[0]
    nu = parameter()[1]
    G = E/(2*(1+nu))
    tol = 1*10**-8
    D = TangentMatrix()
    stresstr = np.zeros(4)
    stressnew = np.zeros(4)
    epnew = np.zeros(4)
    stresstr[0:3] = stressold[0:3] + np.matmul(D, straininc[0:3])
    #stresstr[3] = stressold[3] + nu*(stresstr[0] + stresstr[1])
    stresstr[3] = stressold[3] + E*nu/((1-2*nu)*(1+nu))*(straininc[0] + straininc[2] )
    pressure = (1/3.0)*(stresstr[0] + stresstr[1] + stresstr[3])
    
    sdev1 = stresstr[0] - pressure
    sdev2 = stresstr[1] - pressure
    sdev12 = stresstr[2]
    sdev3 = stresstr[3] - pressure
    snorm = (sdev1**2 + sdev2**2 +sdev3**2 + 2* sdev12**2)**(1/2)
    qtrial = (3/2)**(1/2)*snorm
    syield = sigmay(eqepold)[0]
    ftrial = qtrial - syield
    #print(ftrial)
    if ftrial/syield <=tol:
        stressnew = stresstr
        epnew = epold
        eqepnew = eqepold
    else:
        syield, H = sigmay(eqepold)
        dgama =0
        fntrial = qtrial -3*G*dgama - syield
        #denom = 3*G +H
        count =0
        while abs(fntrial/syield) >=tol:
            count+=1
           # H = sigmay(eqepnew)[1] 
            ddgama = -(qtrial-3*G*dgama-syield)/(-3*G - H)
            dgama = dgama + ddgama 
            eqepnew = eqepold + ddgama
            syield,H = sigmay(eqepnew)
            fntrial = qtrial -3*G*dgama - syield
            print(abs(fntrial))
        s1 = (1- dgama*3*G/qtrial)*sdev1
        s2 =(1- dgama*3*G/qtrial)*sdev2
        s12 = (1- dgama*3*G/qtrial)*sdev12
        s3 =(1- dgama*3*G/qtrial)*sdev3
        #print(s1)
       # print('plastic')
        stressnew[0] = s1 + pressure
        stressnew[1] = s2 + pressure
        stressnew[2] = s12 
        stressnew[3] = s3 + pressure
        #eqepnew = eqepold + dgama 
        snewnorm = (s1**2 + s2**2 +s3**2 + 2* s12**2)**(1/2)
        epnew[0] = epold[0] + dgama* (3/2)**(1/2)*s1/snewnorm
        epnew[1] = epold[1] + dgama* (3/2)**(1/2)*s2/snewnorm
        epnew[2] = epold[2] + dgama* (3/2)**(1/2)*s12/snewnorm
        epnew[3] = epold[3] + dgama* (3/2)**(1/2)*s3/snewnorm
        
        
    return stressnew, epnew, eqepnew


straininc, stressold, epold, eqepold = np.zeros(4), np.zeros(4), np.zeros(4), 0
stress = np.array([])
strain =np.array([])
straininc1 =0 
str1 = np.array([])
n =5000
for i in range(n):
    straininc[0] = .0006
    straininc[1] = -0.3*.0006
    straininc[3] = 0#-0.3*.0006
    snew, epnew, eqepnew = stressupdate(straininc, stressold, epold, eqepold)
    stressold, epold, eqepold = snew, epnew, eqepnew
    stress= np.append(stress, snew)
    strain = np.append(strain, epnew)
    straininc1 = straininc1 + straininc[0]
    str1 =  np.append(str1, straininc1)
stress= stress.reshape(n, 4)
strain = strain.reshape(n, 4)
plt.plot(str1, stress[:, 0])   
    
    
    
    