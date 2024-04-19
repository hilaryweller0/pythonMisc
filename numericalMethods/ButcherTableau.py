import numpy as np
import sys
import matplotlib.pyplot as plt

# Calculating the order of accuracy of an RK scheme from the Butcher Tableau
class ButcherTableau:
    def __init__(self, name, A, w):
        self.name = name
        self.A = np.array(A)
        self.w = np.array(w).reshape(-1)
        self.n = len(A)
        self.c = sum(self.A.T)
        
        # Consistency checks
        if self.A.shape[0] != self.n or self.A.shape[1] != self.n:
            print(self.n)
            raise ValueError("ButcherTableau: A should be nxn not "
                          +str(self.A.shape[0])+'x'+str(self.A.shape[1]))
        if self.n != len(self.w):
            raise ValueError("ButcherTableau: Size of A and w should be the same, not "
                         +str(self.n)+' and '+str(len(self.w)))
    
    def order(self, orderConditions=None):
        first = sum(self.w)
        second = sum(self.c*self.w)
        wac = np.dot(self.w, np.dot(self.A, self.c))
        wcc = sum(self.w*self.c*self.c)
        
        o1 = max(1 - abs(first-1), 0)
        o2 = max(1 - abs(2*second - 1), 0)
        o3 = max(0.5 - abs(6*wac - 1), 0) + max(0.5 - abs(3*wcc - 1), 0)
        
        order =  o1 + np.floor(o1+sys.float_info.epsilon)*(o2 + np.floor(o2+sys.float_info.epsilon)*o3)

        if isinstance(orderConditions, list):
            orderConditions += [first]*(4-len(orderConditions))
            orderConditions[0] = first
            orderConditions[1] = second
            orderConditions[2] = wcc
            orderConditions[3] = wac
        else:
            print('First order condition', first, 'shoud be 1')
            print('Second order condition', second, 'should be 0.5')
            print('Third order condition wac is ', 6*wac, '/6 should be 1/6', sep='')
            print('Third order condition wcc is ', 3*wcc, '/3 should be 1/3', sep='')
        
        return order
    
    def Lstable(self):
        stable=False
        det = np.linalg.det(self.A)
        if det > sys.float_info.epsilon:
            L = np.dot(self.w, np.dot(np.linalg.inv(self.A), np.ones(self.n)))
            print('L-stable condition', L, 'should be 1')
            stable = abs(L-1) < sys.float_info.epsilon
        else:
            print('det(A) =', det)
        return stable
    
    def describe(self):
        print('RK scheme', self.name)
        print('A =', self.A)
        print('w =', self.w)
        print('c =', self.c)
        print('order', self.order())
        print('L-stable', self.Lstable())
    
    def stabilityFunction(self, z = np.linspace(0,10,101)):
        S = np.zeros_like(z)
        I = np.identity(self.n)
        for i in range(len(z)):
            M = I + z[i]*self.A
            if abs(np.linalg.det(M)) > sys.float_info.epsilon:
                S[i] = 1 - z[i]*np.dot(self.w, np.dot(np.linalg.inv(M), np.ones(self.n)))
        return S
    
    def stabilityFunctionImag(self, z = np.linspace(0,10,101)):
        S = np.zeros_like(z)
        I = np.identity(self.n)
        for i in range(len(z)):
            M = I - 1j*z[i]*self.A
            if abs(np.linalg.det(M)) > sys.float_info.epsilon:
                S[i] = abs(1 + 1j*z[i]*np.dot(self.w, np.dot(np.linalg.inv(M), np.ones(self.n))))
        return S
    
    def addZeroRow(self):
        A = np.concatenate(([np.zeros(self.n)], self.A))
        A = np.concatenate((A, np.array([np.zeros(self.n+1)]).T), axis=1)
        w = np.concatenate(([0], self.w))
        return ButcherTableau(self.name+str(0), A, w)

class doubleButcher:
    def __init__(self, B1, B2):
        self.B1 = B1
        self.B2 = B2
        
        # Consistency checks
        if self.B1.n != self.B2.n:
            raise ValueError("doubleButcher: should contain two ButcherTableau of the same size, not "
                          +str(self.B1.n)+' and '+str(self.B2.n))
    
    def order(self):
        B1 = self.B1
        B2 = self.B2
        o1 = B1.order()
        o2 = B2.order()
        
        oc2_12 = sum(B1.c*B2.w)
        oc2_21 = sum(B2.c*B1.w)
        print('Second order coupled conditions are', oc2_12, 'and', oc2_21, 'should be 0.5')
        
        wac112 = np.dot(B1.w, np.dot(B1.A, B2.c))
        wac121 = np.dot(B1.w, np.dot(B2.A, B1.c))
        wac122 = np.dot(B1.w, np.dot(B2.A, B2.c))
        wac212 = np.dot(B2.w, np.dot(B1.A, B2.c))
        wac221 = np.dot(B2.w, np.dot(B2.A, B1.c))
        wac211 = np.dot(B2.w, np.dot(B1.A, B1.c))
        
        wcc112 = sum(B1.w*B1.c*B2.c)
        wcc121 = sum(B1.w*B2.c*B1.c)
        wcc122 = sum(B1.w*B2.c*B2.c)
        wcc212 = sum(B2.w*B1.c*B2.c)
        wcc221 = sum(B2.w*B2.c*B1.c)
        wcc211 = sum(B2.w*B1.c*B1.c)
        
        print('Third order wac conditions are ', 6*wac112,'/6, ',6*wac121, '/6, ' ,6*wac122,  
              '/6, ', 6*wac212, '/6, ',6*wac221, '/6, ', 6*wac211, '/6. Should be 1/6', sep='')
        print('Third order wcc conditions are ', 3*wcc112, '/3, ',3*wcc121, '/3, ',3*wcc122, 
              '/3, ', 3*wcc212, '/3, ', 3*wcc221, '/3, ',3*wcc211, '/3. Should be 1/3', sep='')
        
    def describe(self):
        self.B1.describe()
        self.B2.describe()


class ShuOsher:
    def __init__(self, name, L0, L1, M0, M1):
        self.name = name
        self.L0 = np.array(L0)
        self.L1 = np.array(L1).reshape(-1)
        self.M0 = np.array(M0)
        self.M1 = np.array(M1).reshape(-1)
        self.s = len(L0)
        
        # Consistency checks
        if self.L0.shape[0] != self.s or self.L0.shape[1] != self.s:
            raise ValueError("ShuOsher: L0 should be nxn not "
                          +str(self.L0.shape[0])+'x'+str(self.L0.shape[1]))
        if self.s != len(self.L1):
            raise ValueError("ShuOsher: Size of L0 and L1 should be the same, not "
                         +str(self.s)+' and '+str(len(self.L1)))
        if self.M0.shape[0] != self.s or self.M0.shape[1] != self.s:
            raise ValueError("ShuOsher: M0 should be "
                            +str(self.s)+"x" +str(self.s) +", not "
                          +str(self.M0.shape[0])+'x'+str(self.M0.shape[1]))
        if self.s != len(self.M1):
            raise ValueError("ShuOsher: Size of M0 and M1 should be the same, not "
                         +str(self.M0.shape[0])+' and '+str(len(self.M1)))
    
    def describe(self):
        print('Shu-Osher RK scheme', self.name)
        print('L0 =', self.L0)
        print('L1 =', self.L1)
        print('M0 =', self.M0)
        print('M1 =', self.M1)
    
    def Butcher(self):
        s = self.s
        invImL0 = np.linalg.inv(np.eye(s) - self.L0)
        A = np.matmul(invImL0, self.M0)
        b = self.M1 + np.matmul(self.L1, A)
        print('A = ', A, ' b = ', b)
        return ButcherTableau(self.name, A, b)

def SDIRK3(s):
    m11 = 0.5*(1 - np.sqrt((s-1)/(s+1)))
    m21 = 0.5*(np.sqrt((s+1)/(s-1)) - 1)
    ms1 = (s+1)/(s*(s+1 + np.sqrt(s**2-1)))
    ls1 = (s+1)*(s-1 + np.sqrt(s**2-1))/(s*(s+1+np.sqrt(s**2-1)))
    M0 = m11*np.eye(s) + m21*np.eye(s, k=-1)
    L0 = np.eye(s, k=-1)
    M1 = np.zeros(s)
    M1[-1] = ms1
    L1 = np.zeros(s)
    L1[-1] = ls1
    return ShuOsher('SDIRK3_'+str(s), L0, L1, M0, M1)
