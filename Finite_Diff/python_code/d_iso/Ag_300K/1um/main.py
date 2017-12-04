#!~/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=16)
class  solve_IBP(object):
    def __init__(self):
        #dimensionless benchmark
        #self.max_iter=16000
        self.max_iter=10000
        self.L=1.0e-5   # 50um
        self.n=1000
        self.h=self.L/self.n
        self.m=0.14
        self.dt=0.001
        self.lmbd=self.dt/(self.h*self.h*self.h*self.h)
        self.A=np.full((self.n, self.n), 0.0, dtype=np.float)
        self.x=np.linspace(0,self.L,self.n)

        omega=10.283/(6.02e23)
        nu= 1/(np.sqrt(3.0)/4.0*(4.08e-8)**2) #/cm^2
        T=300.0
        D,gamma,kbT=3.3e-13,0.77*1e-4,1.38e-23*T
        #self.B=1.0e-18   #cm4/s
        self.B=D*gamma*omega*omega*nu/kbT   #cm4/s
        self.main_iso()

        plt.legend(bbox_to_anchor=(0.85, 1), loc=2, borderaxespad=0.)
        plt.show()
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

    def print_profile(self,i):
        filename=str(i)+'.out'
        with open(filename,'w') as fout:
            for i in range(self.n):
                fout.write('%12.8f   %12.8f\n'%(self.x[i]*1e7,self.yn[i][0]*1e7))
        return

    def main_iso(self):
        self.get_A_1()
        self.initial_yn()
        self.x=np.linspace(0,self.L,self.n)
        for i in range(self.max_iter):
            self.yn = np.linalg.solve(self.A, self.ynm1)
            self.ynm1=np.copy(self.yn)
            self.ynm1[0][0]=-self.h*self.m
            self.ynm1[1][0]=-self.h*self.m
            self.ynm1[self.n-1][0]=0.0
            #print self.yn
            if (i+1)%2000 == 0:
                self.print_profile(i+1)
                plt.plot(self.x*1e7,self.yn*1e7,label=str(i+1))
        return
 
    
    

if __name__ == "__main__":
    sln=solve_IBP()
