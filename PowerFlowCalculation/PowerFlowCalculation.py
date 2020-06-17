﻿# -*- coding: UTF-8 -*-
'''
Name:         PowerFlowCalculation.py
Func:         To calculate the power flow of power grid using
              Newton-Raphson method
              Derived from author Maples7's Github project
Author:       4thFever
Addr:         Beijing, China
Time:         2020-06-16 
Link:         https://github.com/4thfever
'''
import numpy as np
import pandas as pd
from math import sin, cos, fabs, pi

class Device(object):
    """docstring for Device"""
    def __init__(self, type_, *args):
        self.type_ = type_
        if type_ == 'load':
            self.i, self.a, self.b = args
        else:
            self.i, self.j, self.a, self.b, self.c = args

# const
 

class PowerFlow():
    def __init__(self):
        self.max_iter = 30               # max iter times
        self.max_error = 1.0e4   
        self.error_max = 1e-5
        self.Um, self.Ua, self.Jacob = None, None, None
        self.P, self.Q, self.Y = None, None, None
        self.num_node, self.num_line, self.num_tran, self.num_gene, self.num_load = 0, 0, 0, 0, 0
        self.Line, self.Tran, self.Gene, self.Load = Device, Device, Device, Device
        self.line, self.tran, self.gene, self.load = [], [], [], []

    def read_data(self):
        """
        Read data from input.csv
        """
        df_grid = pd.read_csv('input.csv')
        self.num_node = int(df_grid.loc[:,['node_i','node_j']].max().max())
        list_attr = [self.num_line, self.num_tran, self.num_gene, self.num_load]
        list_device = ['line', 'tran', 'gene', 'load']
        list_handle_device = [self.Line, self.Tran, self.Gene, self.Load]
        list_handles_device = [self.line, self.tran, self.gene, self.load]
        for i in range(4):
            list_attr[i] = int(df_grid[df_grid['type']==list_device[i]].shape[0])
        self.num_line, self.num_tran, self.num_gene, self.num_load = list_attr
        print(self.error_max, self.num_line, self.num_tran, self.num_gene, self.num_load) 

        for type_num in range(4):
            df_buf = df_grid[df_grid['type']==list_device[type_num]]
            df_buf = df_buf.dropna(axis=1, how='all')
            values = df_buf.iloc[:,1:].values
            list_handles_device[type_num].append(None)
            for i in range(values.shape[0]):
                if 'node_j' in df_buf.columns:
                    input_ = [list_device[type_num]] + [int(ele) if i <= 1 else ele for i, ele in enumerate(values[i])]
                else:
                    input_ = [list_device[type_num]] + [int(ele) if i == 0 else ele for i, ele in enumerate(values[i])]
                device_buf = list_handle_device[type_num](*input_)
                list_handles_device[type_num].append(device_buf)

    def admt_matrix(self):
        """
        Create admittance matrix
        """
        self.Y = np.zeros((self.num_node+1, self.num_node+1), dtype = complex)
        for lineNum in range(1, self.num_line+1):
            i = self.line[lineNum].i
            j = self.line[lineNum].j
            r = self.line[lineNum].a
            x = self.line[lineNum].b
            comp = 1/complex(r, x)
            if i==j:
                self.Y[i][i] += comp
            else:
                c = self.line[lineNum].c
                self.Y[i][j] -= comp
                self.Y[j][i] = self.Y[i][j]
                self.Y[i][i] += (comp + complex(0, c))
                self.Y[j][j] += (comp + complex(0, c))
        for tranNum in range(1, self.num_tran+1):
            i = self.tran[tranNum].i
            j = self.tran[tranNum].j
            r = self.tran[tranNum].a
            x = self.tran[tranNum].b
            c = self.tran[tranNum].c
            comp = 1/complex(r, x)
            self.Y[i][i] += comp
            self.Y[i][j] -= comp/c
            self.Y[j][i] = self.Y[i][j]
            self.Y[j][j] += comp/c/c


    def Um_and_Ua(self):
        """
        Set the amplitude and phase angle of voltage
        """
        self.Um = np.ones(self.num_node+1)
        self.Ua = np.zeros(self.num_node+1)
        for i in range(1, self.num_gene+1):
            if self.gene[i].j <= 0:
                self.Um[self.gene[i].i] = self.gene[i].c

    def form_Jacobian(self):
        """
        Form Jacobian Matrix & Calc the Power error
        """
        n2 = 2*self.num_node
        nu = n2 + 1
        for i in range(1, self.num_node+1):
            vi = self.Um[i]
            di = self.Ua[i]
            dp = 0.0
            dq = 0.0
            for j in range(1, self.num_node+1):
                if j != i:                  # when i <> j, off-diagonal elements
                    g = self.Y[i][j].real        # G        
                    b = self.Y[i][j].imag        # B
                    vj = self.Um[j]
                    dj = self.Ua[j]
                    dij = di - dj           # diff of Phase Angle
                    Hij = -self.Um[i] * self.Um[j] * (g*sin(dij) - b*cos(dij))
                    Lij = Hij
                    self.Jacob[i][j] = Hij
                    self.Jacob[i+self.num_node][j+self.num_node] = Lij
                    Nij = -self.Um[i]*self.Um[j]*(g*cos(dij)+b*sin(dij))
                    Mij = -Nij
                    self.Jacob[i][j+self.num_node]=Nij
                    self.Jacob[i+self.num_node][j] = Mij
                    p = self.Um[j]*(g*cos(dij)+b*sin(dij))
                    q = self.Um[j]*(g*sin(dij)-b*cos(dij))
                    dp += p
                    dq += q
            g = self.Y[i][i].real
            b = self.Y[i][i].imag
            Hii = vi*dq
            Nii = -vi*dp - 2*vi*vi*g
            Mii = -vi*dp
            Lii = -vi*dq + 2*vi*vi*b
            self.Jacob[i][i] = Hii
            self.Jacob[i][i+self.num_node] = Nii
            self.Jacob[i+self.num_node][i] = Mii
            self.Jacob[i+self.num_node][i+self.num_node] = Lii
            self.Jacob[i][nu] = -vi*(dp+vi*g)
            self.Jacob[i+self.num_node][nu] = -vi*(dq-vi*b)
            self.P[i] = vi * (dp+vi*g)
            self.Q[i] = vi * (dq-vi*b)
        for i in range(1, self.num_load+1):
            kk = self.load[i].i
            lp = self.load[i].a
            lq = self.load[i].b
            self.Jacob[kk][nu] += -lp
            self.Jacob[kk+self.num_node][nu] += -lq
        for i in range(1, self.num_gene+1):
            kk = self.gene[i].i
            gp = self.gene[i].a
            gq = self.gene[i].b
            self.Jacob[kk][nu] += gp
            self.Jacob[kk+self.num_node][nu] += gq
        for k in range(1, self.num_gene+1):
            ii = self.gene[k].i
            kk = self.gene[k].j
            if kk == 0:         # Balance nodes
                for j in range(1, n2+1):
                    self.Jacob[ii][j] = 0.0
                    self.Jacob[self.num_node+ii][j] = 0.0
                    self.Jacob[j][ii] = 0.0
                    self.Jacob[j][self.num_node+ii] = 0.0
                self.Jacob[ii][ii] = 1.0
                self.Jacob[self.num_node+ii][self.num_node+ii] = 1.0
                self.Jacob[ii][nu] = 0.0
                self.Jacob[self.num_node+ii][nu] = 0.0
            if kk < 0:          # PV nodes
                for j in range(1, n2+1):
                    self.Jacob[self.num_node+ii][j] = 0.0
                    self.Jacob[j][self.num_node+ii] = 0.0
                self.Jacob[self.num_node+ii][self.num_node+ii] = 1.0
                self.Jacob[self.num_node+ii][nu] = 0.0

    def node_flow(self):
        """
        output the power flow of nodes
        """
        print("\n\n\n\t\t* - * - * - Rasult of Power Flow Calculation * - * - * -")
        print("\n\n\t\t\t\t-------power flow of nodes-------")
        print("\n\n\tno.i\t Um\t\tUa\t\tPG\t\t  QG\t PL\tQL\n")
        for i in range(1, self.num_node+1):
            b1,b2,c1,c2 = 0.0, 0.0, 0.0, 0.0
            for j in range(1, self.num_gene+1):
                ii = self.gene[j].i
                kk = self.gene[j].j
                if i == ii and kk == 0:         # Balance nodes
                    b1 = self.P[ii]
                    b2 = self.Q[ii]
                    for k in range(1, self.num_load+1):
                        ii = self.load[k].i
                        if i == ii:
                            c1 = self.load[k].a
                            c2 = self.load[k].b
                            b1 += c1
                            b2 += c2
                    break
                if i == ii and kk == -1:            # PV nodes
                    b1 = self.gene[j].a
                    b2 = self.Q[ii]
                    for k in range(1, self.num_load+1):
                        ii = self.load[k].i
                        if i == ii:
                            c1 = self.load[k].a
                            c2 = self.load[k].b
                            b2 += c2
                    break 
            for j in range(1, self.num_load+1):
                ii = self.load[j].i
                if i == ii:
                    c1 = self.load[j].a
                    c2 = self.load[j].b
                    break
            print(" %6d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %(i, self.Um[i], self.Ua[i]*180.0/pi, b1, b2, c1, c2))

    def branch_flow(self):
        """
        output the power flow of branches
        """
        print("\n\n\t\t\t\t-------power flow of branches-------")
        print("\n\n\ti\t j\t\tPij\t\t   Qij\t\t  Pji\t\t Qji\t\t dP\t\t   dQ\n\n")
        ph, qh = 0.0, 0.0
        for p in self.line:
            if p == None:
                continue
            i = p.i
            j = p.j
            r = p.a
            x = p.b
            b = r*r + x*x
            if i == j:
                vi = self.Um[i]
                b = vi*vi/b
                pij = r*b
                qij = x*b
                pji = 0.0
                qji = 0.0
                dpb = pij
                ph += dpb
                dqb = qij
                qh += dqb
            else:
                r = r/b
                x = -x/b
                b = p.c
                dij = self.Ua[i] - self.Ua[j]
                vi = self.Um[i]
                vj = self.Um[j]
                vij = vi*vj
                vi *= vi
                vj *= vj
                cd = vij * cos(dij)
                sd = vij * sin(dij)
                pij = vi*r - r*cd - x*sd
                pji = vj*r - r*cd + x*sd
                dpb = pij + pji
                ph += dpb
                qij = -vi*(b+x) + x*cd - r*sd
                qji = -vj*(b+x) + x*cd + r*sd
                dqb = qij + qji
                qh += dqb
            print("  %3d  %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %(i, j, pij, qij, pji, qji, dpb, dqb))
        for p in self.tran:
            if p == None:
                continue
            i = p.i
            j = p.j
            r = p.a
            x = p.b
            t = p.c
            b = t*(r*r+x*x)
            r /= b
            x /= -b
            b = t - 1.0
            ri = r*b
            xi = x*b
            rj = -ri/t
            xj = -xi/t
            vi = self.Um[i]
            vj = self.Um[j]
            vij = vi*vj
            vi *= vi
            vj *= vj
            dij = self.Ua[i] - self.Ua[j]
            cd = vij * cos(dij)
            sd = vij * sin(dij)
            pij = vi*(ri+r) - r*cd - x*sd
            pji = vj*(rj+r) - r*cd + x*sd
            dpb = pij + pji
            ph += dpb
            qij = -vi*(xi+x) + x*cd - r*sd
            qji = -vj*(xj+x) + x*cd + r*sd
            dqb = qij + qji
            qh += dqb
            print("  %3d  %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" %(i, j, pij, qij, pji, qji, dpb, dqb))
        print("\n\n  The total loss of the system: - Active power:%8.5f\tReactive power:%8.5f" %(ph, qh))

    def solv_Eqn(self):
        """
        solve the Modified Equations
        """
        Jacob = self.Jacob
        n2 = 2*self.num_node
        nu = n2 + 1
        for i in range(1, n2+1):
            i1 = i+1
            d = 1.0/Jacob[i][i]
            for j in range(i1, nu+1):
                e = Jacob[i][j]
                if e != 0.0:
                    Jacob[i][j] = e*d
            if i != n2:
                for j in range(i1, n2+1):
                    e = Jacob[j][i]
                    if e != 0.0:
                        for k in range(i1, nu+1):
                            Jacob[j][k] -= Jacob[i][k]*e
        for k in range(2, n2+1):
            i = n2 - k + 1
            i1 = i + 1
            for j in range(i1, n2+1):
                Jacob[i][nu] = Jacob[i][nu] - Jacob[i][j]*Jacob[j][nu]

pf = PowerFlow()
# main 
pf.read_data()
pf.admt_matrix()
pf.Um_and_Ua()

pf.Jacob = np.zeros((2*pf.num_node+1, 2*pf.num_node+2))
pf.P = np.zeros(pf.num_node+1)
pf.Q = np.zeros(pf.num_node+1)

iter = 0  
x_axis = []             # for drawing graph
y_axis = []             # for drawing graph
while True:
    pf.form_Jacobian()
    error = 0.0
    for i in range(1, 2*pf.num_node+1):
        if fabs(pf.Jacob[i][2*pf.num_node+1]) > error:
            error = fabs(pf.Jacob[i][2*pf.num_node+1])
    print("Times of iteration: %2d\t\tThe maximum power error: %11.6f" %(iter+1, error))
    #print("%d %.6f\n" %(iter+1, error))
    x_axis.append(iter+1)           # for drawing graph   
    y_axis.append(error)            # for drawing graph
    if error < pf.error_max:
        pf.node_flow()
        pf.branch_flow()
        break
    if iter > pf.max_iter or error > pf.max_error:
        print("\n\n\t\tThe power flow is Divergence.")
        break
    pf.solv_Eqn()
    for i in range(1, pf.num_node+1):
        a = pf.Jacob[i][2*pf.num_node+1]
        pf.Ua[i] += -a
        a = pf.Jacob[pf.num_node+i][2*pf.num_node+1]
        pf.Um[i] *= 1-a
    iter += 1