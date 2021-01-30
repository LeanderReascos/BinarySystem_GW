from system_dinamics import *
import sys

file_name = 'bodies_input.txt'

Bodies = []

M_Sun = 1.9891e30
UA = 1.4710e11
c = 299792458

unidades = {'Msun':M_Sun,'UA':UA,'c':c}

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def convert(string):
   num = ''
   uni = ''
   for s in string:
      if is_number(s) or s == '-' or s == '.' or s=='e':
         num += s
      else:
         uni += s
   if uni != '':
      return float(num)*unidades[uni]
   else:
      return float(num)

def makeSystem(file_name):

    Bodies = []
    with open(file_name,'r') as f:
       line = f.readline()
       while line:
          name,mass,r0,v0,color,_ = line.split(',')
          _,n = name.split(':')
          _,m = mass.split(':')
          _,r = r0.split(':')
          _,v = v0.split(':')
          _,color = color.split(':')

          x,y,z = r.split(' ')
          vx,vy,vz = v.split(' ')
          x = convert(x)
          y = convert(y)
          z = convert(z)
          vx = convert(vx)
          vy = convert(vy)
          vz = convert(vz)
          m = convert(m)

          body = Particle(m,[x,y,z],[vx,vy,vz],name=n,color=color)
          print(body.get_mass(),body.get_position(),body.get_velocity(),body.get_Name(),body.get_Color())
          Bodies.append(body)
          line = f.readline()

    return System_2Bodies(body1=Bodies[0],body2=Bodies[1])

def export_data(System,tf,h):
    t = 0
    while t <= tf:
        sys.stdout.write(f'\r{np.round(t/tf*100,2)}/{100}%')
        sys.stdout.flush()
        r,args = System.evolution_dt(h)
        if np.isnan(r).any():
          print('NAN ERROR',t/tf)
          break
        System.calc_vectors(r,args)
        t += h

    T = np.arange(0,t+h,h)
    Xs = System.get_freqsGW()
    fig1,ax1 = plt.subplots()
    i = len(Xs)
    ax1.plot(T[:i],Xs)

    Bodies = System.get_bodies()
    distance = []
    r1 = Bodies[0].get_positionArray()
    r2 = Bodies[1].get_positionArray()

    for i,r in enumerate(r1):
        distance.append(calc_distance(r,r2[i]))
    i = len(Xs)
    fig2,ax2 = plt.subplots()
    ax2.plot(T[:i],distance)

    fig3 = plt.figure()
    ax3 = plt.axes(projection='3d')
    for i,body in enumerate(Bodies):
            r = body.get_positionArray()
            ax3.plot3D(r[:,0],r[:,1],r[:,2],color=body.get_Color())
            ax3.scatter3D(r[-1,0],r[-1,1],r[-1,2],color=body.get_Color(),label=body.get_Name())

    ax3.legend()

    plt.show()

    System.export_data(tf,h)

export_data(makeSystem(file_name),1.6,8e-4)
