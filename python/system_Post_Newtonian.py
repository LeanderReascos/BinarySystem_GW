import numpy as np
import matplotlib.pyplot as plt

G = 6.674184e-11 #Constante da gravitacao universal
c = 299792458 #Velocidade da Luz

def calc_magnitude(r):
    return np.sqrt(np.sum(r**2))

def calc_distance(r1,r2):
    r = r2-r1
    return calc_magnitude(r)


class Particle():
    def __init__(self,mass,position=[0,0,0],velocity=[0,0,0],name=None,color=None,append_Array = True):
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
        

class System():
    def __init__(self,num_particles=0,M=None,R=None,V=None,x0=None):
        self.__Particles = []
        for i in range(num_particles):
            self.__Particles.append(Particle(M[i],R[i],V[i]))
        
    def calc_centerMass(self):
        M = 0
        u_ = 0
        R = np.zeros(3)
        for p in self.__Particles:
            m = p.get_mass()
            M += m
            u_ += 1/m
            R += m*p.get_position()
        
        self.__M = M
        self.__R = R
        self.__u = 1/u_
        
    def add_Particle(self,particle):
        self.__Particles.append(particle)
    
    def calc_MotionEquation(self,ks=None):
        if ks is None:
            ks = np.array([[[0,0,0]]*len(self.__Particles),[[0,0,0]]*len(self.__Particles)])
        dR = []
        dV = []
        kr = ks[0]
        kv = ks[1]
        for i,A in enumerate(self.__Particles):
            dvA = 0
            r_A = A.get_position()+kr[i]
            v_A = A.get_velocity()+kv[i]
            dR.append(v_A)

            for B in self.__Particles:
                sum1B = 0
                sum2B = 0
                sum3B = 0
                sum4B = 0
                r_B = B.get_position()
                v_B = B.get_velocity()
                r_AB = r_A-r_B
                if A != B:
                    potB = -B.get_NewtonianPotential(r_A,1)
                    sum2B = potB
                    sum1B = potB*r_AB/calc_magnitude(r_AB)**2
                    sum3B = potB/calc_magnitude(r_AB)**2*np.dot(r_AB,(4*v_A-3*v_B))*(v_A-v_B)/c**2
                    sum4B = 5*A.get_NewtonianPotential(r_B,1)/c**2 + calc_magnitude(v_A)**2/c**2 - 4*np.dot(v_A,v_B)/c**2 \
                            + 2*calc_magnitude(v_B)**2/c**2-3/2*(np.dot(v_B,r_AB)/(c*calc_magnitude(r_AB)))**2
                sum1C = 0
                sum2C = 0
                sum3C = 0
                for C in self.__Particles:
                    if A != C:
                        r_C = C.get_position()
                        v_C = C.get_velocity()
                        r_AC = r_A-r_C  
                        potC = -C.get_NewtonianPotential(r_A,1)
                        sum1C += potC/c**2
                        if C != B:
                            r_BC = r_B-r_C
                            potCB = -C.get_NewtonianPotential(r_B,1)
                            sum2C += -potCB/c**2 + 1/2*potCB/c**2*np.dot(r_AB,r_BC)/calc_magnitude(r_BC)**2
                            sum3C += potCB/c**2/calc_magnitude(r_BC)**2*r_BC
                            
                dvA += -sum1B*(1-4*sum1C+sum2C+sum4B)-7/2*sum2B*sum3C+sum3B

            dV.append(dvA)
        return np.array([dR,dV],float)
    
    
    def get_actualR(self):
        R = []
        V = []
        for p in self.__Particles:
            R.append(p.get_position())
            V.append(p.get_velocity())
        return np.array([R,V],float)
    
    def evolution_dt(self,h):
        k0 = h*self.calc_MotionEquation()
        k1 = h*self.calc_MotionEquation(k0/2)
        k2 = h*self.calc_MotionEquation(k1/2)
        k3 = h*self.calc_MotionEquation(k2)
        
        return self.get_actualR()+1/6*(k0+2*k1+2*k2+k3)
    
    def modify_R(self,r):
        R = r[0]
        V = r[1]
        for i,p in enumerate(self.__Particles):
            p.set_position(R[i])
            p.set_velocity(V[i])
