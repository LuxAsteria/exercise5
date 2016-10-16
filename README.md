# exercise5
ps:All pictures in this assignment are drawn by myself
***
###Abstract
Here in this time's homework, the Euler method is applied to this question, which requires calculating the cannon shell's trajectories. This exercise takes 5 steps. First, enumerating the trajectory under the circumstances that ignoring the friction of natrual wind, taking the default B<sub>2</sub> value, too. Then, taking the natural wind's influence into cosideration. Last but not least, discussing the situation that the value of B<sub>2</sub> changing with the temperature.
***
###Background
**Realistic Projectile Motion__Question 2.11 and 2.12**:
  >Calculate cannon shell trajectories applying Euler method, using the adiabatic model of air density and considering the effect of the Earth's revolution about its own axis( Coriolis force)
***
###Exercise
####Analysis:
For a start, establishing a regtangular coordinates system as shown in the picture. Given our planet is approximately a spheroid, so we chose the system as shown in the following picture.
![v coordinate](https://github.com/LuxAsteria/test3/blob/master/v%20coordinates.png)

and we can analyze the velocity of wind in this way,too.
![vwind coordinate](https://github.com/LuxAsteria/test3/blob/master/vwind%20coordinates.png)

Here are several factors that will affect the movement of a cannon shell. First, air friction.This friction is associated with the relative velocity between the cannon shell and the wind.
As is illustrated in the textbook<sup>1</sup>, we have:
![latex](https://github.com/LuxAsteria/test3/blob/master/f%20drag%20common.png)

Here, V<sup>2</sup> refers to the square of the composite velocity of wind and the cannon shell.
F<sub>drag</sub> could be affected by pressure and temperature.
Consider the effect of pressure, we have[1]:

![latex rho0](https://github.com/LuxAsteria/test3/blob/master/f%20drag%20pressure.png)

And we have<img src="http://latex.codecogs.com/gif.latex?\frac{\rho}{\rho_{0}}=(1-\frac{az}{T_{0}})^{\alpha}" alt="" title="" />

Then,consider the effect of temperature, we have:

![latex temperature](https://github.com/LuxAsteria/test3/blob/master/temperature.png)
Where we have <img src="http://latex.codecogs.com/gif.latex?B_{2}^{ref}=0.00004/m" alt="" title="" />

Next step, think about the Coriolis force.This force is caused by the autorotation of Earth and the movement of cannon. When the velocity of the cannon is not perpendicular to the angular velocity, the force won't be zero and would always be orthogonal to the composite velocity.
Here we have:
<img src="http://latex.codecogs.com/gif.latex?F_{cor}=-2m\vec{\omega}\times\vec{v}" alt="" title="" />
where m is the mass of the cannon shell.

####Getting the exressions:
From the *second Law of Newton's*,it's not difficult to get:

<img src="http://latex.codecogs.com/gif.latex?\vec{F}=m/vec{a}=m\frac{d\vec{v}}{dt}=\vec{F_{drag}}+\vec{F_{cor}}+\vec{G}" alt="" title="" />
Project the velocity to 3 coordinates when we ignore the effect of natural wind, we have[2]:

<img src="http://latex.codecogs.com/gif.latex?\frac{dV_{x}}{mdt}=-(1-\frac{az}{T_{0}})^{\alpha}B_{2}VV_{x}-2sin(\gamma)\omega\frac{\sqrt{2}}{2}(V_{y}-V_{z})" alt="" title="" />

<img src="http://latex.codecogs.com/gif.latex?\frac{dV_{y}}{mdt}=-(1-\frac{az}{T_{0}})^{\alpha}B_{2}VV_{y}+2\omega(sin(\gamma)V_{x}\frac{\sqrt{2}}{2}+cos(\gamma)V_{z})" alt="" title="" />

<img src="http://latex.codecogs.com/gif.latex?\frac{dV_{z}}{mdt}=-(1-\frac{az}{T_{0}})^{\alpha}B_{2}VV_{z}-2\omega(sin(\gamma)V_{x}\frac{\sqrt{2}}{2}+cos(\gamma)V_{y})-g" alt="" title="" />

Thinking about the natural wind, so some terms should be added to the former formula.It's easy to get those terms, just projecting the V<sub>wind</sub> to 3 coordinates, and taking them in the formula [1], appending them to the end of every expression in [2].

####The introduction of Euler method had been mentioned in [exercise 4](https://github.com/LuxAsteria/exercise-4), for more details please click the superlink.
***
####step 1: Considering about the temperature and the natural wind
![t_w](https://github.com/LuxAsteria/test3/blob/master/t_w%2019(wz)%20100(z).png)

Here's the code:
```
class cannons_c:
    """with no coriolis force"""
    def __init__(self,v=700,T=300,angz=100,angxy=30,timestep=0.1,time_duration=120,omeg=0.0000727,b=0.00004,k=10000,apha=2.5,a=0.0065,g=9.8,anglati=40,angwxy=40,angwz=139,vwind=6,Tr=340):
        self.angz=[angz]
        self.angxy=[angxy]
        self.anglati=anglati
        self.v=[v]
        self.T=T
        self.time=time_duration
        self.t=[0]
        self.dt=timestep
        self.vx=[v*math.cos(angxy)*math.sin(angz)]
        self.vy=[v*math.sin(angxy)*math.sin(angz)]
        self.vz=[v*math.cos(angz)]
        self.xx=[0]
        self.xy=[0]
        self.xz=[0]
        self.angxy=[angxy]
        self.angz=[angz]
        self.omeg=omeg
        self.nsteps=int(time_duration//timestep+1)
        self.b=b
        self.k=k
        self.apha=apha
        self.a=a
        self.g=g
        self.angwxy=angwxy
        self.angwz=angwz
        self.vwind=vwind
        self.Tr=Tr
        self.vcx=[v*math.cos(angxy)*math.sin(angz)]
        self.vcy=[v*math.sin(angxy)*math.sin(angz)]
        self.vcz=[v*math.cos(angxy)]
        self.xc=[0]
        self.yc=[0]
        self.zc=[0]
        self.vc=[]
    def cal_c(self):
        for i in range(self.nsteps):
            if self.zc[i]>=0:
                cv=math.sqrt((self.vcx[i])**2+(self.vcy[i])**2+(self.vcz[i])**2)
                self.vc.append(cv)
                cvx=self.vcx[i]-self.dt*self.b*((1-self.a*self.zc[i]/self.Tr)**(self.apha))*((self.Tr/self.T)**(self.apha))*((self.vc[i]*self.vcx[i])+(self.vwind)**2*math.sin(self.angwz)*math.cos(self.angwxy))
                self.vcx.append(cvx)
                cvy=self.vcy[i]-self.dt*self.b*(1-self.a*self.zc[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*((self.vc[i]*self.vcy[i])+(self.vwind)**2*math.sin(self.angwz)*math.sin(self.angwxy))
                self.vcy.append(cvy)
                cvz=self.vcz[i]-self.dt*self.b*(1-self.a*self.zc[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*((self.vc[i]*self.vcz[i])+(self.vwind)**2*math.cos(self.angwz))-self.dt*self.g
                self.vcz.append(cvz)
                xc=self.xc[i]+self.dt*self.vcx[i]
                yc=self.yc[i]+self.dt*self.vcy[i]
                zc=self.zc[i]+self.dt*self.vcz[i]
                self.xc.append(xc)
                self.yc.append(yc)
                self.zc.append(zc)
            else:
                xc=self.xc[i]+self.dt*self.vcx[i]
                yc=self.yc[i]+self.dt*self.vcy[i]
                zc=self.zc[i]+self.dt*self.vcz[i]
                self.xc.append(xc)
                self.yc.append(yc)
                self.zc.append(zc)
                break
    def show(self):
        fig=plt.figure(1)
        ax=fig.gca(projection='3d')
        ax.plot(self.xc,self.yc,self.zc)
        ax.set_zlim3d(0,10000)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
y=cannons_c()
y.cal_c()
y.show()
```
####step 2: Taking Coriolis into consideration but exclude the influence of wind :
Consider the force in the following coordinates:
![coriolis](https://github.com/LuxAsteria/test3/blob/master/coriolis%20force.png)

Here's the screenshot:
![t_c](https://github.com/LuxAsteria/test3/blob/master/t_c_39(wz)100(z)%203d.png)
![t_c2d](https://github.com/LuxAsteria/test3/blob/master/t_c_139(wz)100(z).png)

Here's the code:
```
    def calculate(self):
        for i in range(self.nsteps):
            if self.xz[i]>=0:
                simv=math.sqrt(self.vx[i]**2+self.vy[i]**2+self.vz[i]**2)
                self.v.append(simv)
                simvx=self.vx[i]-self.dt*(self.b*(1-self.a*self.xz[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*(self.vx[i]*self.v[i]))+self.dt*2*self.omeg*(math.sin(self.anglati)*math.sqrt(2)/2*self.vz[i]+self.vy[i]*math.cos(self.anglati))
                self.vx.append(simvx)
                simvy=self.vy[i]-self.dt*(self.b*(1-self.a*self.xz[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*(self.vy[i]*self.v[i]))-self.dt*2*self.omeg*(self.vx[i]*math.cos(self.anglati)+self.vz[i]*math.sin(self.anglati)*math.sqrt(2)/2)
                self.vy.append(simvy)
                simvz=self.vz[i]-self.dt*self.g-self.dt*(self.b*(1-self.a*self.xz[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*(self.vz[i]*self.v[i]))+2*self.dt*self.omeg*(math.sin(self.anglati)*self.vy[i]*math.sqrt(2)/2-self.vx[i]*math.sqrt(2)/2*math.sin(self.anglati))
                self.vz.append(simvz)
                simx=self.xx[i]+self.dt*self.vx[i]
                simy=self.xy[i]+self.dt*self.vy[i]
                simz=self.xz[i]+self.dt*self.vz[i]
                self.xx.append(simx)
                self.xy.append(simy)
                self.xz.append(simz)
            else:
                simx=self.xx[i]+self.dt*self.vx[i]
                simy=self.xy[i]+self.dt*self.vy[i]
                simz=self.xz[i]+self.dt*self.vz[i]
                self.xx.append(simx)
                self.xy.append(simy)
                self.xz.append(simz)
                break
   
```
####step 3:Considering Coriolis force, natural wind and the effect of temperature:

The screenshot of the simulation:

![pic](https://github.com/LuxAsteria/test3/blob/master/t_c_w_39(wz)100(z)_3d.png)
![pic2d](https://github.com/LuxAsteria/test3/blob/master/t_c_w%2039(wz)100(z).png)

Here's the code:
```
 def calculate(self):
        for i in range(self.nsteps):
            if self.xz[i]>=0:
                simv=math.sqrt(self.vx[i]**2+self.vy[i]**2+self.vz[i]**2)
                self.v.append(simv)
                simvx=self.vx[i]-self.dt*(self.b*(1-self.a*self.xz[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*(self.vx[i]*self.v[i]+self.vwind**2*math.sin(self.angwz)*math.cos(self.angwxy)))+self.dt*2*self.omeg*(math.sin(self.anglati)*math.sqrt(2)/2*self.vz[i]+self.vy[i]*math.cos(self.anglati))
                self.vx.append(simvx)
                simvy=self.vy[i]-self.dt*(self.b*(1-self.a*self.xz[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*(self.vy[i]*self.v[i]+self.vwind**2*math.sin(self.angwz)*math.sin(self.angwxy)))-self.dt*2*self.omeg*(self.vx[i]*math.cos(self.anglati)+self.vz[i]*math.sin(self.anglati)*math.sqrt(2)/2)
                self.vy.append(simvy)
                simvz=self.vz[i]-self.dt*self.g-self.dt*(self.b*(1-self.a*self.xz[i]/self.Tr)**self.apha*(self.Tr/self.T)**self.apha*(self.vwind**2*math.cos(self.angwz)+self.vz[i]*self.v[i]))+2*self.dt*self.omeg*(math.sin(self.anglati)*self.vy[i]*math.sqrt(2)/2-self.vx[i]*math.sqrt(2)/2*math.sin(self.anglati))
                self.vz.append(simvz)
                simx=self.xx[i]+self.dt*self.vx[i]
                simy=self.xy[i]+self.dt*self.vy[i]
                simz=self.xz[i]+self.dt*self.vz[i]
                self.xx.append(simx)
                self.xy.append(simy)
                self.xz.append(simz)
            else:
                simx=self.xx[i]+self.dt*self.vx[i]
                simy=self.xy[i]+self.dt*self.vy[i]
                simz=self.xz[i]+self.dt*self.vz[i]
                self.xx.append(simx)
                self.xy.append(simy)
                self.xz.append(simz)
                break
```

####step4:Comparing the data:
First, get the maximum of height and distance:
```
def maxi(self):
        print(max(self.xz))
        print(max(self.xx))
```
The data are put in the table:

First, show the differeces caused by changes in angel beteewn velocity(wind velocity) and the z coordinate

![table](https://github.com/LuxAsteria/test3/blob/master/table.png)

Second, show the differences caused by the different angle between angular velocity and z coordinate:
![table](https://github.com/LuxAsteria/test3/blob/master/lati.png)

####setp 5:Found in pictures:
![pic](https://github.com/LuxAsteria/test3/blob/master/t_c_39(wz)120(z)3d.png)

It's apparent that the trajectory would tilt when we set the simulation in 3d projection and taking the following factors:Coriolis force, air firction,etc. into consideration.

***
###Conclusion
In this assignment,I've simulated the trajectory of a cannon shell, considering the factors following: air friction,including the fluent caused by cannon shell's movement and natural wind, Coriolis force, the effect of temperature and pressure. Conspicuously, if we exclude the influence of natural wind, the trajectory projected to the xy plane is approximately a line. However, if we take it into consideration, that's definitely not a line but a curve. If the Coriolis force is eliminated, as is shown above,the distance and height would be shorter, or in other words, *the trajectory would be flatter*, so does the wind's influence is excluded.However, the latter one's influence is smaller relatively.
Besides, the effect of the change in angels is conspicuous,too.
Although I've considered as much factors as possible, I still neglected one factor:the curve of the ground, which means the enumeration of maximum distance may not be such specific. In comparation with the radius of Earth, however, if the accuracy of statistic is not demanded, this factor could still be put aside.
***
###Acknowledge
[1]Computational Physics, Edition 2, Nicholas J.Giordano & Hisao Nakanishi
[2]Professor Cai Hao's basic codes

###PS:All the pictures are drawn by myself
#PPS:ALL RIGHTS RESERVED!
