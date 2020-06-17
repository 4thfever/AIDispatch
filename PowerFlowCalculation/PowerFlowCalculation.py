# -*- coding: UTF-8 -*-
'''
Name:         PowerFlowCalculation.py
Func:         To calc the power flow of elecrtcal grid using
              Newton-Raphson method, this is the main file of this project
              Derived from author:Maples7
Author:       4thFever
Addr:         Beijign, China
Time:         2020-06-16 
Link:         ?
'''
# import globalVariable as gv
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
MAX_ITER = 30               # max iter times
DIVERGENCE_ERROR = 1.0e4    

class PowerFlow():
    def __init__(self):
        self.num_node = 0
        self.num_line = 0
        self.num_tran = 0
        self.num_gene = 0
        self.num_load = 0
        self.Line, self.Tran, self.Gene, self.Load = Device, Device, Device, Device
        self.line, self.tran, self.gene, self.load = [], [], [], []

    def read_data(self):
        """
        Read Data from input.txt, the global variables are in the globalVariable.py
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
        self.error_max = 1e-5
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

    def output_file_ready(self):
        """
        make output.txt file ready to store output result
        """
        global fou
        try:
            fou = open("output.txt", "w+")
        except:
            print("Error: Can't create output.txt. (function output_file_ready() error)")
            exit()

    def admt_matrix(self):
        """
        Create admittance matrix
        """
        global Y
        Y = np.zeros((self.num_node+1, self.num_node+1), dtype = complex)
        for lineNum in range(1, self.num_line+1):
            i = self.line[lineNum].i
            j = self.line[lineNum].j
            r = self.line[lineNum].a
            x = self.line[lineNum].b
            comp = 1/complex(r, x)
            if i==j:
                Y[i][i] += comp
            else:
                c = self.line[lineNum].c
                Y[i][j] -= comp
                Y[j][i] = Y[i][j]
                Y[i][i] += (comp + complex(0, c))
                Y[j][j] += (comp + complex(0, c))
        for tranNum in range(1, self.num_tran+1):
            i = self.tran[tranNum].i
            j = self.tran[tranNum].j
            r = self.tran[tranNum].a
            x = self.tran[tranNum].b
            c = self.tran[tranNum].c
            comp = 1/complex(r, x)
            Y[i][i] += comp
            Y[i][j] -= comp/c
            Y[j][i] = Y[i][j]
            Y[j][j] += comp/c/c


    def Um_and_Ua(self):
        """
        Set the amplitude and phase angle of voltage
        """
        global Um
        global Ua
        Um = np.ones(self.num_node+1)
        Ua = np.zeros(self.num_node+1)
        for i in range(1, self.num_gene+1):
            if self.gene[i].j <= 0:
                Um[self.gene[i].i] = self.gene[i].c

    def form_Jacobian(self):
        """
        Form Jacobian Matrix & Calc the Power error
        """
        global Um
        global Ua
        global Jacob
        global P 
        global Q
        global Y

        n2 = 2*self.num_node
        nu = n2 + 1
        for i in range(1, self.num_node+1):
            vi = Um[i]
            di = Ua[i]
            dp = 0.0
            dq = 0.0
            for j in range(1, self.num_node+1):
                if j != i:                  # when i <> j, off-diagonal elements
                    g = Y[i][j].real        # G        
                    b = Y[i][j].imag        # B
                    vj = Um[j]
                    dj = Ua[j]
                    dij = di - dj           # diff of Phase Angle
                    Hij = -Um[i] * Um[j] * (g*sin(dij) - b*cos(dij))
                    Lij = Hij
                    Jacob[i][j] = Hij
                    Jacob[i+self.num_node][j+self.num_node] = Lij
                    Nij = -Um[i]*Um[j]*(g*cos(dij)+b*sin(dij))
                    Mij = -Nij
                    Jacob[i][j+self.num_node]=Nij
                    Jacob[i+self.num_node][j] = Mij
                    p = Um[j]*(g*cos(dij)+b*sin(dij))
                    q = Um[j]*(g*sin(dij)-b*cos(dij))
                    dp += p
                    dq += q
            g = Y[i][i].real
            b = Y[i][i].imag
            Hii = vi*dq
            Nii = -vi*dp - 2*vi*vi*g
            Mii = -vi*dp
            Lii = -vi*dq + 2*vi*vi*b
            Jacob[i][i] = Hii
            Jacob[i][i+self.num_node] = Nii
            Jacob[i+self.num_node][i] = Mii
            Jacob[i+self.num_node][i+self.num_node] = Lii
            Jacob[i][nu] = -vi*(dp+vi*g)
            Jacob[i+self.num_node][nu] = -vi*(dq-vi*b)
            P[i] = vi * (dp+vi*g)
            Q[i] = vi * (dq-vi*b)
        for i in range(1, self.num_load+1):
            kk = self.load[i].i
            lp = self.load[i].a
            lq = self.load[i].b
            Jacob[kk][nu] += -lp
            Jacob[kk+self.num_node][nu] += -lq
        for i in range(1, self.num_gene+1):
            kk = self.gene[i].i
            gp = self.gene[i].a
            gq = self.gene[i].b
            Jacob[kk][nu] += gp
            Jacob[kk+self.num_node][nu] += gq
        for k in range(1, self.num_gene+1):
            ii = self.gene[k].i
            kk = self.gene[k].j
            if kk == 0:         # Balance nodes
                for j in range(1, n2+1):
                    Jacob[ii][j] = 0.0
                    Jacob[self.num_node+ii][j] = 0.0
                    Jacob[j][ii] = 0.0
                    Jacob[j][self.num_node+ii] = 0.0
                Jacob[ii][ii] = 1.0
                Jacob[self.num_node+ii][self.num_node+ii] = 1.0
                Jacob[ii][nu] = 0.0
                Jacob[self.num_node+ii][nu] = 0.0
            if kk < 0:          # PV nodes
                for j in range(1, n2+1):
                    Jacob[self.num_node+ii][j] = 0.0
                    Jacob[j][self.num_node+ii] = 0.0
                Jacob[self.num_node+ii][self.num_node+ii] = 1.0
                Jacob[self.num_node+ii][nu] = 0.0

    def node_flow(self):
        """
        output the power flow of nodes
        """
        global fou
        global P
        global Q
        global Um
        global Ua
        fou.write("\n\n\n\t\t* - * - * - Rasult of Power Flow Calculation * - * - * -")
        fou.write("\n\n\t\t\t\t-------power flow of nodes-------")
        fou.write("\n\n\tno.i\t Um\t\t\tUa\t\t\tPG\t\t  QG\t\t PL\t\t\tQL\n\n")
        for i in range(1, self.num_node+1):
            b1,b2,c1,c2 = 0.0, 0.0, 0.0, 0.0
            for j in range(1, self.num_gene+1):
                ii = self.gene[j].i
                kk = self.gene[j].j
                if i == ii and kk == 0:         # Balance nodes
                    b1 = P[ii]
                    b2 = Q[ii]
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
                    b2 = Q[ii]
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
            fou.write(" %6d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %(i, Um[i], Ua[i]*180.0/pi, b1, b2, c1, c2))

    def branch_flow(self):
        """
        output the power flow of branches
        """
        global Um
        global Ua
        fou.write("\n\n\t\t\t\t-------power flow of branches-------")
        fou.write("\n\n\ti\t j\t\tPij\t\t   Qij\t\t  Pji\t\t Qji\t\t dP\t\t   dQ\n\n")
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
                vi = Um[i]
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
                dij = Ua[i] - Ua[j]
                vi = Um[i]
                vj = Um[j]
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
            fou.write("  %3d  %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %(i, j, pij, qij, pji, qji, dpb, dqb))
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
            vi = Um[i]
            vj = Um[j]
            vij = vi*vj
            vi *= vi
            vj *= vj
            dij = Ua[i] - Ua[j]
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
            fou.write("  %3d  %3d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %(i, j, pij, qij, pji, qji, dpb, dqb))
        fou.write("\n\n  The total loss of the system: - Active power:%8.5f\t\tReactive power:%8.5f" %(ph, qh))

    def solv_Eqn(self):
        """
        solve the Modified Equations
        """
        global Jacob
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
pf.output_file_ready()
pf.admt_matrix()
pf.Um_and_Ua()

Jacob = np.zeros((2*pf.num_node+1, 2*pf.num_node+2))
P = np.zeros(pf.num_node+1)
Q = np.zeros(pf.num_node+1)

iter = 0  
x_axis = []             # for drawing graph
y_axis = []             # for drawing graph
while True:
    pf.form_Jacobian()
    error = 0.0
    for i in range(1, 2*pf.num_node+1):
        if fabs(Jacob[i][2*pf.num_node+1]) > error:
            error = fabs(Jacob[i][2*pf.num_node+1])
    fou.write("Times of iteration: %2d\t\tThe maximum power error: %11.6f\n" %(iter+1, error))
    #fou.write("%d %.6f\n" %(iter+1, error))
    x_axis.append(iter+1)           # for drawing graph   
    y_axis.append(error)            # for drawing graph
    if error < pf.error_max:
        pf.node_flow()
        pf.branch_flow()
        break
    print(iter, MAX_ITER, error, DIVERGENCE_ERROR)
    if iter > MAX_ITER or error > DIVERGENCE_ERROR:
        fou.write("\n\n\t\tThe power flow is Divergence.")
        break
    pf.solv_Eqn()
    for i in range(1, pf.num_node+1):
        a = Jacob[i][2*pf.num_node+1]
        Ua[i] += -a
        a = Jacob[pf.num_node+i][2*pf.num_node+1]
        Um[i] *= 1-a
    iter += 1

# pf.close_file()

