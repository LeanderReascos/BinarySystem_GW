import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd
import sys


G = 6.674184e-11 #Constante da gravitacao universal
c = 299792458 #Velocidade da Luz
gamma_E = 0.577216 #Euler Constant

def calc_mag(r):
    return np.sqrt(np.sum(r**2))

def calc_versor(r):
    return r/calc_mag(r)

def calc_distance(r1,r2):
    r = r2-r1
    return calc_mag(r)

def vectorial_product(r1,r2):
    vx,vy,vz = r1
    ux,uy,uz = r2
    x = vy*uz-vz*uy
    y = -(vx*uz-vz*ux)
    z = vx*uy-vy*ux
    return np.array([x,y,z])


class Particle():
    def __init__(self,mass,position=[0,0,0],velocity=[0,0,0],name=None,color=None,append_Array = False):
        self.__r = np.array(position)
        self.__v = np.array(velocity)
        self.__m = mass
        self.__name = name
        self.__color = color
        if append_Array:
            self.__positions = [np.array(position)]
            self.__velocities = [np.array(velocity)]
        else:
            self.__positions = []
            self.__velocities = []
    
    def get_Color(self):
        return self.__color
    
    def get_Name(self):
        return self.__name
    
    def get_NewtonianPotential(self,r,m):
        return -G*m*self.__m/calc_distance(self.__r,r)
    
    def get_KineticEnergy(self):
        m_A = self.__m
        v_A = self.__v
        return m_A*(c**2-1/2*v_A**2-1/8*v_A**4/c**2)
        
    def get_mass(self):
        return self.__m
    
    def get_velocity(self):
        return self.__v
    
    def get_position(self):
        return self.__r
    
    def get_positionArray(self):
        return np.array(self.__positions)
    
    def get_velocityArray(self):
        return np.array(self.__velocities)
    
    def set_position(self,r):
        self.__r = r
        self.__positions.append(r)
        
    def set_velocity(self,v):
        self.__v = v
        self.__velocities.append(v)


class System_2Bodies():
    
    def __init__(self,body1=None,body2=None,f0=None,w=np.array([0,0,1]),file=" "):
        if file == " ":
            m1 = body1.get_mass()
            m2 = body2.get_mass()
            self.__M = m1+m2
            self.__u = m1*m2/self.__M
            self.__n = self.__u/self.__M
            self.__R = (m1*body1.get_position()+m2*body2.get_position())/self.__M
            self.__a = body2.get_position()-body1.get_position()
            
            if f0 is None:
                self.__v = body2.get_velocity()-body1.get_velocity()#
                v = calc_mag(self.__v)
                self.__x = v**2/c**2*(1+2/3*(3-self.__n)*v**2/c**2)#(G*self.__M*v/(c**3*calc_mag(self.__a)))**(2/3)#
                x0 = self.__x
                f0 = x0/(2*np.pi)
                tc= 5/256*G*self.__M/(c**3*self.__n)*(np.pi*G*self.__M*f0/c**3)**(-8/3)
                print(x0,v,tc)
            else:
                self.__v = vectorial_product(w,calc_versor(self.__a))
                x0 = 2*np.pi*f0
                tc= 5/256*G*self.__M/(c**3*self.__n)*(np.pi*G*self.__M*f0/c**3)**(-8/3) #5/256*G*self.__M/(c**3*self.__n)*1/x0**4*(1+(743/252+11/3*self.__n)*x0)
                print(x0,tc)
            
            self.__x0 = x0
            self.__Bodies = [body1,body2]
            phi = 0
            self.__phi = 0
            self.__xs = []
            self.calc_vectors([self.__a,self.__v],[x0,phi])
        else:
            self.import_data(file)

    def import_data(self,file_name):
        with open(file_name,'r') as f:
            line = f.readline()
            count_bodies = 0
            body = False
            values = np.array(['BODY','DATA','SIMULATION'])
            aux = np.arange(1,4)
            Bodies = []
            self.__time = []
            self.__xs = []
            zone = ""
            while line:
                vals = line.split(',')
                if (values == vals[0]).any():
                    zone = values[np.sum((values == vals[0])*aux)-1]
                    line = f.readline()
                    continue 

                if zone == 'BODY':
                    count_bodies += 1
                    name,mass,color = [v.split(':')[1] for v in vals[:-1]]
                    Bodies.append(Particle(mass=float(mass),name=name,color=color))

                elif zone == 'DATA':
                    t,x,y,z,vx,vy,vz = [float(v) for v in vals[:-1]]
                    self.__time.append(t)
                    r = np.array([x,y,z])
                    v = np.array([vx,vy,vz])
                    Bodies[count_bodies-1].set_position(r)
                    Bodies[count_bodies-1].set_velocity(v)

                elif zone == 'SIMULATION':
                    if vals[0].split(':')[0] == 'Time_Simulation':
                        self.__tf, self.__h = [float(v.split(':')[1]) for v in vals[:-1]]
                    elif vals[0].split(':')[0] == 'x0':
                        self.__x0 =  float(vals[0].split(':')[1])
                    else:
                        x,phi = [float(v) for v in vals[:-1]]
                        self.__xs.append(np.array([x,phi]))
                line = f.readline()
            self.__Bodies = Bodies
            f.close()
        print('Read file')

    def get_tf(self):
        return self.__tf

    def get_h(self):
        return self.__h

    def get_x0(self):
        return self.__x0

    def get_m1m2(self):
        return [b.get_mass() for b in self.__Bodies]

    def export_data(self,tf,h):
        file = open('Export_Data.txt','w')
        for i,body in enumerate(self.__Bodies):
            print(f'BODY,',file=file)
            print(f'Name:{body.get_Name()},m{i+1}:{body.get_mass()},color:{body.get_Color()},',file=file)
            r = body.get_positionArray()
            r = np.array(r)
            v = body.get_velocityArray()
            df = pd.DataFrame(r,columns=['XX','YY','ZZ'])
            df = df.dropna(how='any')
            t = np.arange(0,tf+2*h,h)
            print(f'DATA,',file=file)
            for j in range(len(df['XX'])):
                x,y,z =[df['XX'][j],df['YY'][j],df['ZZ'][j]]
                print(f'{t[j]},{x},{y},{z},{v[j][0]},{v[j][1]},{v[j][2]},',file=file)
            tf = t[j]    
        print(f'SIMULATION,',file=file)  
        print(f'Time_Simulation:{tf},h:{h},',file=file)   
        print(f'x0:{self.__x0},',file=file)   
        freqs = self.get_freqsGW()
        for f,p in freqs[:j]:
            print(f'{f},{p},',file=file)
        
    def calc_MotionEquation(self,ks=None):
        if ks is None:
            ks = np.array([[0,0,0],[0,0,0]])
        kr = ks[0]
        kv = ks[1]
        
        a = self.__a+kr
        v = self.__v+kv
        dr = v
        M = self.__M
        n = self.__n
        mag_v = calc_mag(v)
        mag_a = calc_mag(a)
        dv = -G*M/mag_a**3*a*(1-(4+2*n)*G*M/(c**2*mag_a)+(1+3*n)*mag_v**2/c**2-3/2*n*(np.dot(v,a))**2/(c**2*mag_a**2))+(4-2*n)*G*M/(c**2*mag_a**3)*np.dot(v,a)*v

        return np.array([dr,dv],float)
        
    def calc_Args(self,ks=None):
        if ks is None:
            ks = np.array([0,0])
        kr = ks[0]
        kv = ks[1]
        
        x = self.__x + kr
        M = self.__M
        n = self.__n
        dx = 64/5*n*c**3/(G*M)*x**5*(1-(743/336+11/4*n)*x+4*np.pi*x**1.5\
                                    +(34103/18144+13661/2016*n+59/18*n**2)*x**2\
                                    -(4159/672+189/8*n)*np.pi*x**2.5\
                                    +(16447322263/139708800+16/3*np.pi**2-1712/105*gamma_E-856/105*np.log(16*x)\
                                     +(-56198689/217728+451/48*np.pi**2)*n+541/896*n**2-5605/2591*n**3)*x**3\
                                    -(4415/4032-358675/6048*n-91495/1512*n**2)*np.pi*x**3.5)
        #print(dx)
        dphi = c**3/(G*M)*x**1.5
        return np.array([dx,dphi],float)

    def get_bodies(self):
        return self.__Bodies

    def get_freqsGW(self):
        return np.array(self.__xs)
    
    def calc_vectors(self,r,args):
        a,v = r
        x,phi = args
        self.__xs.append(args)
        M = self.__M
        n = self.__n
        self.__x = x
        self.__phi = phi
        
        mag_v = np.sqrt(x*c**2*(1-2/3*(3-n)*x))
        w = x**1.5*c**3/(G*M)
        mag_a = mag_v/w
        
        self.__v = mag_v*calc_versor(v)
        self.__a = mag_a*calc_versor(a)
        
        b1,b2 = self.__Bodies
        c1 = -b2.get_mass()/M
        c2 = b1.get_mass()/M
        cs = [c1,c2]
        
        for i,b in enumerate(self.__Bodies):
            b.set_position(cs[i]*self.__a)
            b.set_velocity(cs[i]*self.__v)
    
    
    def evolution_dt(self,h):
        k0 = h*self.calc_MotionEquation()
        k1 = h*self.calc_MotionEquation(k0/2)
        k2 = h*self.calc_MotionEquation(k1/2)
        k3 = h*self.calc_MotionEquation(k2)
        
        a0 = h*self.calc_Args()
        a1 = h*self.calc_Args(a0/2)
        a2 = h*self.calc_Args(a1/2)
        a3 = h*self.calc_Args(a2)
        
        args = np.array([self.__x,self.__phi])
        r = np.array([self.__a,self.__v])
        return [r+1/6*(k0+2*k1+2*k2+k3),args+1/6*(a0+2*a1+2*a2+a3)]
    


