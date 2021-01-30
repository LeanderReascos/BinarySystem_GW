from system_dinamics import *
import sys
from spin_weighted_Yml import *

cos = np.cos
sin = np.sin
sqrt = np.sqrt
pi = np.pi

G = 6.674184e-11 #Constante da gravitacao universal
c = 299792458 #Velocidade da Luz
gammaE = 0.577216 #Euler Constant


    

def calc_WaveForm(r,theta,phi,psi,x,x0,m1,m2):
    M = m1+m2
    u = m1*m2/M
    n = u/M
    dm = m1-m2
    
    ls = [2,3,4]
    
    H = 0+0j
    
    h22 = -sqrt(pi/5)*G*u/(c**2*r)*np.exp(-2j*psi)*x*(\
            1-(107/42-55/42*n)*x + (2*pi+6j*np.log(x/x0))*x**1.5\
            - (2173/1512+1096/216*n-2047/1512*n**2)*x**2\
            -((107/21-34/21*n)*pi+24j*pi+1j*(107/7-34/7*n)*np.log(x/x0))*x**2.5\
            +(27027409/646800-856/105*gammaE+2/3*pi**2-1712/105*np.log(2)-428/105*np.log(x)\
             -18*(np.log(x/x0))**2-(278185/33264-41/96*pi**2)*n-20261/2772*n**2+114635/99792*n**3\
             +1j*428/105*pi+12j*pi*np.log(x/x0))*x**3)
    h21 = -1j*8/3*sqrt(pi/5)*G*u/(c**2*r)*dm/M*np.exp(-1j*psi*x**1.5)
    h33 = 3j*sqrt(pi/7)*G*u*dm/(c**2*r*M)*np.exp(-3j*psi)*x**1.5
    h32 = -8/3*sqrt(pi/7)*G*u/(c**2*r)*np.exp(-2j*psi)*(1-3*n)*x**2
    h31 = -1j/3*sqrt(2*pi/35)*G*u*dm/(c**2*r*M)*np.exp(-1j*psi)*x**1.5
    h44 = 64/9*sqrt(pi/7) *G*u/(c**2*r)*np.exp(-4j*psi)*(1-3*n)*x**2
    h42 = -8/63*sqrt(pi)*G*u/(c**2*r)*np.exp(-2j*psi)*(1-3*n)*x**2

    hs = {'h21':h21,'h22':h22,'h31':h31,'h32':h32,'h33':h33,'h44':h44,'h42':h42}
    
    for l in ls:
        if l != 4:
            for m in range(1,l+1):
                h = hs['h'+str(l)+str(m)]
                H += sYlm(-2,l,m,theta,phi)*h + sYlm(-2,l,-m,theta,phi)*(-1)**(-l)*np.conjugate(h)
        else:
            h = hs['h'+str(l)+'2']
            m = 2
            H += sYlm(-2,l,m,theta,phi)*h + sYlm(-2,l,-m,theta,phi)*(-1)**(-l)*np.conjugate(h)
            h = hs['h'+str(l)+'4']
            m = 4
            H += sYlm(-2,l,m,theta,phi)*h + sYlm(-2,l,-m,theta,phi)*(-1)**(-l)*np.conjugate(h)
    return H



S = System_2Bodies(file='Export_Data.txt')
args = S.get_freqsGW()
x = args[:,0]
psi = args[:,1]
m1,m2 = S.get_m1m2()
H = calc_WaveForm(4e5,0,0,psi,x,S.get_x0(),m1,m2)
print(len(H))

T = np.arange(0,S.get_tf(),S.get_h())
plt.plot(T,np.abs(H)*np.sin(x*T))
plt.show()