#!~/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=16)

class  solve_IBP(object):
    def __init__(self):
        #dimensionless benchmark
        self.max_iter=100
        self.L=5.0
        self.n=200
        self.h=self.L/self.n
        self.m=0.1
        self.B=2.0

        #diffusion rate and surface energy
        self.D0=1.0
        self.del_D=0.1
        self.gm0=1.0 
        self.del_gm=0.1
        self.eps=0.001


        self.Bt=[self.B]*self.n
        self.slp=[0.0]*self.n
        self.dt=0.001
        self.lmbd=self.dt/(self.h*self.h*self.h*self.h)
        self.A=np.full((self.n, self.n), 0.0, dtype=np.float)
        #self.show_spectral_radius()
        self.x=np.linspace(0,self.L,self.n)
        self.main_iso()

        self.main_aniso()
        #plt.legend(bbox_to_anchor=(0.85, 1), loc=2, borderaxespad=0.)
        #plt.show()
        return
    def show_spectral_radius(self):
        self.get_A_1()
        eign=np.linalg.eigvals(self.A)
        print eign
        print np.amax(eign)
        return
    def initial_yn(self):
        self.ynm1=np.full((self.n,1), 0.0, dtype=np.float)
        self.ynm1[0][0]=-self.h*self.m
        self.ynm1[1][0]=-self.h*self.m
        return
    def get_A_1(self):
        for i in range(self.n):
            if (i>1) and (i<(self.n-2)):
                self.A[i][i-2]=self.B*self.lmbd
                self.A[i][i-1]=-4.0*self.B*self.lmbd
                self.A[i][i]  =6.0*self.B*self.lmbd+1.0
                self.A[i][i+1]=-4.0*self.B*self.lmbd
                self.A[i][i+2]=self.B*self.lmbd
            elif i==0 :
                self.A[i][0]=1.0
                self.A[i][1]=-1.0
            elif i==1:
                self.A[i][1]=2.0
                self.A[i][2]=-3.0
                self.A[i][3]=1.0
            #elif i==2:
            #    self.A[i][2]=self.B*self.lmbd+1.0
            #    self.A[i][3]=-4.0*self.B*self.lmbd
            #    self.A[i][4]=6.0*self.B*self.lmbd
            #    self.A[i][5]=-4.0*self.B*self.lmbd
            #    self.A[i][6]=self.B*self.lmbd
            #elif i==3:
            #    self.A[i][2]=self.B*self.lmbd
            #    self.A[i][3]=-4.0*self.B*self.lmbd+1.0
            #    self.A[i][4]=6.0*self.B*self.lmbd
            #    self.A[i][5]=-4.0*self.B*self.lmbd
            #    self.A[i][6]=self.B*self.lmbd
            elif i==(self.n-2):
                self.A[i][i-3]=self.B*self.lmbd
                self.A[i][i-2]=-4.0*self.B*self.lmbd
                self.A[i][i-1]=6.0*self.B*self.lmbd
                self.A[i][i]  =-4.0*self.B*self.lmbd+1.0
                self.A[i][i+1]=self.B*self.lmbd
            elif i==self.n-1:
                self.A[i][i]=1.0
        return
    def get_A_t(self):
        for i in range(self.n):
            if (i>1) and (i<(self.n-2)):
                self.A[i][i-2]=self.Bt[i]*self.lmbd
                self.A[i][i-1]=-4.0*self.Bt[i]*self.lmbd
                self.A[i][i]  =6.0*self.Bt[i]*self.lmbd+1.0
                self.A[i][i+1]=-4.0*self.Bt[i]*self.lmbd
                self.A[i][i+2]=self.Bt[i]*self.lmbd
            elif i==0 :
                self.A[i][0]=1.0
                self.A[i][1]=-1.0
            elif i==1:
                self.A[i][1]=2.0
                self.A[i][2]=-3.0
                self.A[i][3]=1.0
            elif i==(self.n-2):
                self.A[i][i-3]=self.Bt[i]*self.lmbd
                self.A[i][i-2]=-4.0*self.Bt[i]*self.lmbd
                self.A[i][i-1]=6.0*self.Bt[i]*self.lmbd
                self.A[i][i]  =-4.0*self.Bt[i]*self.lmbd+1.0
                self.A[i][i+1]=self.Bt[i]*self.lmbd
            elif i==self.n-1:
                self.A[i][i]=1.0
        return
    def main_iso(self):
        self.get_A_1()
        self.initial_yn()
        for i in range(self.max_iter):
            self.yn = np.linalg.solve(self.A, self.ynm1)
            self.ynm1=np.copy(self.yn)
            self.ynm1[0][0]=-self.h*self.m
            self.ynm1[1][0]=-self.h*self.m
            self.ynm1[self.n-1][0]=0.0
            #print self.yn
            if (i+1)%100 == 0:
                self.print_profile(i+1,'.iso')
                #plt.plot(self.x,self.yn,label='iso')
        return
    def main_aniso(self):
        self.initial_yn()
        self.get_Bt()
        self.get_A_t()
        for i in range(self.max_iter):
            self.get_Bt()
            self.get_A_t()
            self.yn = np.linalg.solve(self.A, self.ynm1)
            self.ynm1=np.copy(self.yn)
            self.ynm1[0][0]=-self.h*self.m
            self.ynm1[1][0]=-self.h*self.m
            self.ynm1[self.n-1][0]=0.0
            #print self.yn
            if (i+1)%100 == 0:
                self.print_profile(i+1,'.ani')
                #plt.plot(self.x,self.yn,label='ani')
        return
    def get_Bt(self):
        self.Bt=[self.B]*self.n  #benchmark
        #calculate slope for each points
        self.get_slp()
        for i in range(self.n):
            self.Bt[i]=((self.D0+np.sin(self.slp[i])*self.del_D)) * (self.gm0+0.5*self.del_gm*(self.eps+np.sin(self.slp[i])**2/self.eps))
        #print np.linalg.eigvals(self.A)
        return
    def get_slp(self):
        for i in range(self.n):
            if (i>0) and (i<(self.n-1)):
                self.slp[i]=np.arctan((-1.0*self.ynm1[i-1][0]+1.0*self.ynm1[i+1][0])/(2.0*self.h))
            elif i==0 :
                self.slp[i]=np.arctan((-3.0*self.ynm1[i][0]+4.0*self.ynm1[i+1][0]-1.0*self.ynm1[i+2][0])/(2.0*self.h))
                #self.slp[i]=np.arctan((-11.0*self.ynm1[i][0]+18.0*self.ynm1[i+1][0]-9.0*self.ynm1[i+2][0]+2.0*self.ynm1[i+3][0])/(6.0*self.h))
            elif i==self.n-1:
                self.slp[i]=np.arctan((1.0*self.ynm1[i-2][0]-4.0*self.ynm1[i-1][0]+3.0*self.ynm1[i][0])/(2.0*self.h))
        return
    def print_profile(self,i,st):
        filename=str(i)+str(st)+'.out'
        with open(filename,'w') as fout:
            for i in range(self.n):
                fout.write('%12.9E   %12.9E\n'%(self.x[i],self.yn[i][0]))
        return



if __name__ == "__main__":
    sln=solve_IBP()
