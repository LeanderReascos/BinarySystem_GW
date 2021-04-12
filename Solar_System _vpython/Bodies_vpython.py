import numpy as np
import vpython as vp
import sys as s

earth_mass = 5.9722e24
f = 10e5
rescale = 100
uA = 1.4710e11/f
uR = 6371e3

def  norm(v):
    return np.sqrt(v.x**2+v.y**2+v.z**2)

class Body:
    def __init__(self,filename):
        with open(filename, 'r') as f:
            line = f.readline()
            print('Reading file: '+filename)
            while line:
                tag,value = line.split(':')
                value = value[:-1]
                if tag == 'name':
                    self.__name = value
                elif tag == 'mass':
                    self.__nEarthMass = float(value)/uA
                elif tag == 'radius':
                    self.__radius = float(value)/uA
                elif tag == 'texture':
                    self.__texture = value
                elif tag == 'color':
                    r,g,b = [float(v) for v in value[1:-1].split(',')]
                elif tag == 'data':
                    self.__r = np.load(value)
                    self.__r /= uA
                    self.__r /= rescale
                line = f.readline()
            f.close()
        self.__vectorColor = vp.vector(r,g,b)
        print(f"\tRadius: {self.__radius}, Pos: {self.__r[0]}")
        self.__body = vp.sphere(pos=vp.vector(self.__r[0,0],self.__r[0,1],0), radius=self.__radius)
        self.__trail = vp.attach_trail(self.__body,radius= self.__radius*0.8,color=vp.vector(r,g,b),retain=2000)

        if self.__texture == 'no':
            self.__body.color = vp.vector(r,g,b)
        else:
            self.__body.texture = self.__texture

    def get_numberPoints(self):
        return len(self.__r)
    
    def set_pos(self,i):
        self.__body.pos = vp.vector(self.__r[i,0],self.__r[i,1],0)
    
    def set_retain(self,val):
        self.__trail.clear()
        self.__trail.stop()
        self.__trail = vp.attach_trail(self.__body,radius= self.__radius*0.8,color=self.__vectorColor,retain=val)

class System:
    def __init__(self,files_bodies, seconds,h,times):
        self.__h = h
        self.__Bodies = []
        self.__times = np.array(times)
        self.__scene = vp.scene
        self.__scene.width = 1880#width
        self.__scene.height = 1040#height
        self.__scene.range = 4
        self.__scene.camera.pos = vp.vector(f/rescale,0,0)
        self.__scene.camera.axis = vp.vector(-1,0,0)
        self.__scene.camera.up = vp.vector(0,0,1)

        for file_name in files_bodies:
            self.__Bodies.append(Body(file_name))

        self.__N =  self.__Bodies[0].get_numberPoints()
        self.__fps = np.array([ times[i]/(h*s) for i,s in enumerate(seconds[:-1])]+[self.__N/seconds[-1]])

    def run(self):
        frame = 0
        i = 0
        v = vp.vector(0,0,0)
        active = False
        dt = [0.1,0.1]
        vel = [2000,100]
        while frame < self.__N:
            s.stdout.write(f'\rframe: {frame} totalFrames: {self.__N} fps: {self.__fps[i]} v = {v}')
            vp.rate(self.__fps[i])
            k = vp.keysdown() # a list of keys that are down
            v = vp.vector(0,0,0)
            a = vp.vector(0,0,0)
            if 'p' in k:
                active = True
            if 'w' in k: 
                v = self.__scene.camera.axis/norm(self.__scene.camera.axis)*vel[i]
            if 's' in k: 
                v = -self.__scene.camera.axis/norm(self.__scene.camera.axis)*vel[i]
            if 'a' in k: 
                aux = self.__scene.camera.axis/norm(self.__scene.camera.axis)*vel[i]
                v = vp.vector(-aux.y,aux.x,aux.z)
            if 'd' in k: 
                aux = self.__scene.camera.axis/norm(self.__scene.camera.axis)*vel[i]
                v = vp.vector(aux.y,-aux.x,aux.z)
            if 'up' in k:
                v.z += vel[i]/10
            if 'down' in k:
                v.z -= vel[i]/10
            self.__scene.camera.pos += v*dt[i]
            if active:    
                for body in self.__Bodies:
                    body.set_pos(frame)
                if i<len(self.__times) and frame*self.__h <= self.__times[i]+self.__h and frame*self.__h >= self.__times[i]-self.__h:
                    for body in self.__Bodies[:5]:
                        body.set_retain(20)
                    i+=1
                frame += 1



names = ['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune']
files = ['data/'+v for v in names]
sys = System(files,[20,25],4,[365*5])
sys.run()




    