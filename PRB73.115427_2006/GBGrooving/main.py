#!/usr/bin/python
import numpy as np
import random as rd
import os,math,sys
#need to change get_gm function for anisotropic, don't forget
#gm1p gm2p as derivitation of surface energy to contact angle theta


class grooving_iso(object):
    '''
    def mesh_adaptive(self):
        for i in range(0,self.mp/4):
            self.msh[i][0]=i/16.0
        for i in range(self.mp/4,self.mp/2):
            self.msh[i][0]=1.0/16.0+(i-self.mp/4)/8.0
        for i in range(self.mp/2,3*self.mp/4):
            self.msh[i][0]=3.0/16.0+(i-self.mp/2)/4.0
        for i in range(3*self.mp/4,self.mp):
            self.msh[i][0]=7.0/16.0+(i-self.mp/4*3)*(9.0/16.0)*4.0/self.mp
        return
    '''
    def mesh_test(self):
        self.msh[0][0]=0.0
        self.msh[0][1]=-1.0
        self.msh[1][0]=1.0
        self.msh[1][1]=-2.0
        self.msh[2][0]=2.0
        self.msh[2][1]=-1.0
        self.msh[3][0]=3.0
        self.msh[3][1]=0.0
        self.msh[4][0]=3.5
        self.msh[4][1]=1.0
        return
    def mesh_uniform(self):
        for i in range(0,self.mp):
            self.msh[i][0]=i*self.L/self.mp
            self.msh[i][1]=0.0
        return
    def update_b_l_2D(self,i):
        #for previous distance msh[i][2],msh[i][3],msh[i][4],msh[i][5]:b1, b2, l1, l2 of element i
        if i==0:
            #mirrow M_{i+1} to M_{i-1} about z-axis
            x0=(-self.msh[i+1][0])
            x1=self.msh[i+1][1]
            v1=[(self.msh[i+1][0]-x0),(self.msh[i+1][1]-x1)] #M_{i+1} - M_{i-1}
            v2=[(self.msh[i][0]  -x0),(self.msh[i][1]- x1)] #M_{i} - M_{i-1}
            v3=[(self.msh[i+1][0]-self.msh[i][0]),(self.msh[i+1][1]-self.msh[i][1])] #M_{i+1} - M_{i}
            self.msh[i][2]=np.dot(v2,v1)/np.linalg.norm(v1)
            self.msh[i][3]=np.dot(v3,v1)/np.linalg.norm(v1)

            self.msh[i][4]=0
            x0=self.msh[i+1][0]-self.msh[i][0]
            x1=self.msh[i+1][1]-self.msh[i][1]
            self.msh[i][5]=np.sqrt(x0**2+x1**2)
        elif i==(self.mp-1):
            #mirrow M_{i-1} to M_{i+1} about axis
            x0=self.L-self.msh[i-1][0]+self.L
            x1=self.msh[i-1][1]
            v1=[(x0 - self.msh[i-1][0]),(x1 - self.msh[i-1][1])] #M_{i+1} - M_{i-1}
            v2=[(self.msh[i][0]-self.msh[i-1][0]),(self.msh[i][1]-self.msh[i-1][1])] #M_{i} - M_{i-1}
            v3=[(x0 - self.msh[i][0]),(x1 - self.msh[i][1])] #M_{i+1} - M_{i}
            self.msh[i][2]=np.dot(v2,v1)/np.linalg.norm(v1)
            self.msh[i][3]=np.dot(v3,v1)/np.linalg.norm(v1)

            x0=self.msh[i][0]-self.msh[i-1][0]
            x1=self.msh[i][1]-self.msh[i-1][1]
            self.msh[i][4]=np.sqrt(x0**2+x1**2)
            self.msh[i][5]=0
            #check this ,for now set all to 1
            self.msh[i][2]=1
            self.msh[i][3]=1
            self.msh[i][4]=1
            self.msh[i][5]=1
        else:
            v1=[(self.msh[i+1][0]-self.msh[i-1][0]),(self.msh[i+1][1]-self.msh[i-1][1])] #M_{i+1} - M_{i-1}
            v2=[(self.msh[i][0]-self.msh[i-1][0]),(self.msh[i][1]-self.msh[i-1][1])] #M_{i} - M_{i-1}
            v3=[(self.msh[i+1][0]-self.msh[i][0]),(self.msh[i+1][1]-self.msh[i][1])] #M_{i+1} - M_{i}
            self.msh[i][2]=np.dot(v2,v1)/np.linalg.norm(v1)
            self.msh[i][3]=np.dot(v3,v1)/np.linalg.norm(v1)

            x0=self.msh[i][0]-self.msh[i-1][0]
            x1=self.msh[i][1]-self.msh[i-1][1]
            self.msh[i][4]=np.sqrt(x0**2+x1**2)
            x0=self.msh[i+1][0]-self.msh[i][0]
            x1=self.msh[i+1][1]-self.msh[i][1]
            self.msh[i][5]=np.sqrt(x0**2+x1**2)
        return
    def get_gm(self,i): # gamma surface energy
        if i == 0:
            gm1=0
            gm2=self.gamma # J/cm^2
            gm1p=0.0
            gm2p=0.0
        else:
            gm1=self.gamma # J/cm^2
            gm2=self.gamma # J/cm^2
            gm1p=0.0
            gm2p=0.0
        return gm1,gm2,gm1p,gm2p
    def get_h_i(self,i):
        if i==0:
            d=abs(self.msh[1][1]-self.msh[0][1])
        elif i == (self.mp-1):
            d=0.0
        else:
            p1=np.array([self.msh[i-1][0],self.msh[i-1][1]])
            p2=np.array([self.msh[i+1][0],self.msh[i+1][1]])
            p3=np.array([self.msh[i][0],self.msh[i][1]])
            d = np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1)
        return d
    def get_mu_i(self,i):
        mu_u0=0.0  #by default
        if i == 0 :
            gm   = self.get_gm(i)
            gm2  = gm[1]
            gm2p = gm[3]
            self.update_b_l_2D(i)
            b2=self.msh[i][3]
            l2=self.msh[i][5]
            h=self.get_h_i(i)
            E_surf=(h*gm2/l2)-(b2*gm2p/l2)
            mu_i =mu_u0 + 2.0/b2*(E_surf-self.gm_GB/2.0/self.m)
        else:
            self.update_b_l_2D(i)
            b1=self.msh[i][2]
            b2=self.msh[i][3]
            l1=self.msh[i][4]
            l2=self.msh[i][5]
            gm   = self.get_gm(i)
            gm1  = gm[0]
            gm2  = gm[1]
            gm1p = gm[2] 
            gm2p = gm[3] 
            h=self.get_h_i(i)
            mu_i = mu_u0 + 2.0/(b1+b2)*(h*(gm1/l1+gm2/l2)+(b2*gm2p/l2-b1*gm1p/l1))
        return mu_i
    def get_mu_i_all(self):
        for i in range(self.mp):
            self.msh[i][6]=self.get_mu_i(i)
        return
    def get_grad_s_mu_i(self,i):
        #i means the line segment between point i and i+1
        #gradient of mu in 2D: du/dx*i+du/dy*j (i,j standard unit vector)
        dx=self.msh[i+1][0]-self.msh[i][0]
        #dy=self.msh[i+1][1]-self.msh[i][1]
        dmudx=((self.msh[i+1][6]-self.msh[i][6])/(self.msh[i+1][0]-self.msh[i][0]))
        #dmudy=((self.msh[i+1][6]-self.msh[i][6])/(self.msh[i+1][1]-self.msh[i][1]))
        #grd=[dmudx,dmudy]
        #n=[dy,-dx]/np.linalg.norm([dy,-dx])
        #n_grd= np.dot(n,grd)
        #return grd-n_grd*n
        return [dmudx,0]
    def get_grad_s_mu_all(self):
        for i in range(self.mp-1):
            self.msh[i][7]= self.get_grad_s_mu_i(i)[0]
        self.msh[self.mp-1][7]=0.0
        return

    def get_flux_pt_i(self,i):
        if i==0:
            #grd_s_mu=self.get_grad_s_mu_i(i)[0]
            flux_n=-(self.D*self.nu*self.omega/self.kbT)*self.msh[0][7]
            #flux_n=-(self.D*self.nu/self.kbT)*grd_s_mu
            flux_x=flux_n
        elif i== (self.mp-1):
            flux_x=0
        elif i== (self.mp-2):
            flux_x=0
        else:
            #grd_s_mu=self.get_grad_s_mu_i(i)[0]
            flux_p=-(self.D*self.nu*self.omega/self.kbT)*(self.msh[i-1][7])
            #flux_p=-(self.D*self.nu/self.kbT)*grd_s_mu
            #grd_s_mu=self.get_grad_s_mu_i(i+1)[0]
            flux_n=-(self.D*self.nu*self.omega/self.kbT)*self.msh[i][7]
            #flux_n=-(self.D*self.nu/self.kbT)*grd_s_mu
            flux_x=flux_p + flux_n
        print flux_x
        return flux_x
    def update_flux_all(self):
        self.get_grad_s_mu_all()
        for i in range(self.mp):
            self.msh[i][8]=self.get_flux_pt_i(i)*self.dt
        return
    def update_coor_all(self):
        for i in range(self.mp):
            dz=self.omega*self.msh[i][8]/(self.msh[i][2]+self.msh[i][3])
            if i == 0:
                self.msh[i][1]+=dz
            elif i == (self.mp-1):
                self.msh[i][1]+=0
            elif i == (self.mp-2):
                self.msh[i][1]+=0
            else:
                x=self.msh[i+1][0]-self.msh[i-1][0]
                y=self.msh[i+1][1]-self.msh[i-1][1]
                n=[y,-x]/np.linalg.norm([y,-x])
                self.msh[i][0]+=n[0]*dz
                self.msh[i][1]+=n[1]*dz
        return
    def update_coor_all_z(self):
        for i in range(self.mp):
            dz=self.omega*self.msh[i][8]/(self.msh[i][2]+self.msh[i][3])
            if i == 0:
                self.msh[i][1]+=dz
            elif i == (self.mp-1):
                self.msh[i][1]+=0
            elif i == (self.mp-2):
                self.msh[i][1]+=0
            else:
                #x=self.msh[i+1][0]-self.msh[i-1][0]
                #y=self.msh[i+1][1]-self.msh[i-1][1]
                #n=[y,-x]/np.linalg.norm([y,-x])
                #self.msh[i][0]+=n[0]*dz
                self.msh[i][1]+=dz
        return

    def main_loop(self):
        self.mesh_uniform()
        #self.print_coor()
        for i in range(self.max_iter):
            #print ("#iteration %i" % i)
            self.get_mu_i_all()
            self.update_flux_all()
            self.update_coor_all()
            #self.print_coor()
            self.output(i)
            #print self.msh[0][0],self.msh[0][1]
        return
    def print_coor(self):
        for i in range(self.mp):
            print self.msh[i][0],self.msh[i][1]
        return
    def output(self,i):
        #if (i+1)%10 == 0:
        if (i+1)%1 == 0:
            filename=str(i+1)+'.out'
            with open(filename,'w')as fin:
                fin.write("#iteration %i\n" % (i+1))
                fin.write("#x(nm)          y(nm)   chemical_potential  chemical_potential_gradient matter_increment\n")
                for j in range(self.mp):
                    #fin.write('%12.9f    %12.9f\n'%(self.msh[j][0]*1.0e7,self.msh[j][1]*1.0e7))
                    fin.write('%12.9E    %12.9E    %12.9E    %12.9E    %12.9E\n'%(self.msh[j][0]*1.0e7,self.msh[j][1]*1.0e7,self.msh[j][6],self.msh[j][7],self.msh[j][8]))
        return

    def __init__(self):
        #distance in one-side
        self.L  = 5.0e-5  #1nm=1e-7cm  #2.5e-5cm=250nm 
        #numbers of meshing points in 1 dimension
        self.mp = 200
        self.max_iter = 5
        #self.D,self.gamma,self.kbT=2.5e-11,0.77*1e-4,1.38e-23*300.0
        self.D,self.gamma,self.kbT=4.3e-4,1.8*1.0e-4,1.38e-23*700.0
        #units: cm^2/s,J/cm^2,J
        self.omega=10.283/(6.02e23)
        self.nu= 1/(np.sqrt(3.0)/4.0*(4.08e-8)**2) #/cm^2
        #msh[i][0],msh[i][1]:x,y coordinates 
        #[2][3][4][5]:b1,b2,l1,l2  
        #[6]:chemical potential mu
        #[7]:chemical potential gradient
        #[8]: matter increment
        self.msh= np.zeros((self.mp,9))  
        #self.gm_GB= 0.154*1e-4  #J/cm^2
        self.gm_GB= 0.7*1e-4  #J/cm^2
        self.m=0.14 
        self.dt=1.0e-9# #unit second
        self.main_loop()
        #self.mesh_uniform()
        #self.mesh_test()
        #for i in range(10):
        #    print self.msh[i][0],self.msh[i][1]
        #    print self.get_h_i(i)
        return

if __name__ == "__main__":
    GBG_1d=grooving_iso()

