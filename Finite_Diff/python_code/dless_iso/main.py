#!~/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=16)
class  solve_IBP(object):
    def __init__(self):
        #dimensionless benchmark
        self.max_iter=10
        self.L=10.0
        self.n=100
        self.h=self.L/self.n
        self.m=0.1
        self.B=2.5
        self.dt=0.01
        self.lmbd=self.dt/(self.h*self.h*self.h*self.h)
        self.A=np.full((self.n, self.n), 0.0, dtype=np.float)
        self.x=np.linspace(0,self.L,self.n)
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
            if i%1 == 0:
                plt.plot(self.x,self.yn,label=str(i))


        return
 
    
    

if __name__ == "__main__":
    sln=solve_IBP()
