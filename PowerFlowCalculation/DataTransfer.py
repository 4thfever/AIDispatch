# -*- coding: UTF-8 -*-
'''
Name:         DataTransfer.py
Func:         To transfer the topology or power flow result 
              to front end
Author:       4thFever
Addr:         Beijing, China
Time:         2020-06-16 
Link:         https://github.com/4thfever
'''
from PowerFlowCalculation import PowerFlow
import numpy as np

# 向前端发送拓扑信息
def export_topo(pf):
    topo = [ele.topo() for ele in pf.list_elements]
    return pf.num_node, np.array(topo)

# 向前端发送潮流结果
def export_pf_res(pf):
    return pf.df_branch, pf.df_node

def test():
    pf = PowerFlow()
    pf.run()
    res = export_topo(pf)
    print(res)
    res = export_pf_res(pf)
    print(res)

test()