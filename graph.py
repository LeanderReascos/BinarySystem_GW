import numpy as np
import vpython as vp
from system_dinamics import *



class System_VPython(System_2Bodies):
    """docstring for System_VPython"""
    def __init__(self, filename,unities_distance=1e3,unities_radius=1.9891e30,width=1000,height=500):
        super(System_VPython, self).__init__(file=filename)

        
        self.__Bodies = self.get_bodies()
        self.__Array = [b.get_positionArray()/unities_distance for b in self.__Bodies]
        self.__N = len(self.__Array[0])
        print("Number of points: ", self.__N)
        max_distance = np.max(np.abs(self.__Array[0]))

        self.__scene = vp.scene
        self.__scene.width = width
        self.__scene.height = height
        self.__scene.range = max_distance*2
        self.__scene.camera.pos = vp.vector(-2*max_distance,0,max_distance)
        self.__scene.camera.axis = vp.vector(2*max_distance,0,-max_distance)
        self.__scene.camera.up = vp.vector(0,0,1)

        aux = [p.get_positionArray()[0]/unities_distance for p in self.__Bodies]
        print(f'Max: {max_distance}\nInitial:\n\tBody1: {aux[0]}\tBody2: {aux[1]}')
        pos_0 = [vp.vector(v[0],v[1],v[2]) for v in aux]
        rad = [p.get_mass()/unities_radius*2 for p in self.__Bodies]
        self.__Spheres = [vp.sphere(pos=pos_0[i], radius=rad[i], retain=20, make_trail= True) for i in range(len(aux))]

    def make_animation(self,fps):
    	import sys
    	if self.__N % fps != 0:
    		for i,a in enumerate(self.__Array):
    			a = list(a)
    			a += [a[-1]]*(self.__N%fps)
    			self.__Array[i] = np.array(a)
    	self.__N = len(self.__Array[0])
    	seconds = self.__N//fps
    	frames = seconds*fps
    	print(f'Second: {seconds}\nFrames: {frames}\nN: {self.__N}')
    	f = 1
    	while f <= frames:
    		vp.rate(fps)
    		sys.stdout.write(f'\r{f}/{frames}')
    		sys.stdout.flush()
    		for i,sphere in enumerate(self.__Spheres):
    			x,y,z = self.__Array[i][f]
    			sphere.pos = vp.vector(x,y,z)
    		f += 1
    		

file_name = 'Export_Data.txt'
S = System_VPython(file_name)
S.make_animation(60)