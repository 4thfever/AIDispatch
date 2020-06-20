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
from PowerFlow import PowerFlow

from pyecharts import options as opts
import numpy as np

# 向前端发送拓扑信息
def export_topo(pf):
    # topo = [ele.topo() for ele in pf.list_elements]
    nodes = [opts.GraphNode(name=f"Node{num}", symbol_size=20) for num in range(1, pf.num_node+1)]
    links = []

    count_load = 1
    for ele in pf.list_elements:
        if ele.type_ == 'load':
            nodes.append(opts.GraphNode(name=f"Lode{count_load}", symbol_size=10))
            links.append(opts.GraphLink(source=f"Lode{count_load}", target=f"Node{ele.i}"))
            count_load += 1
        else:
            links.append(opts.GraphLink(source=f"Node{ele.i}", target=f"Node{ele.j}"))
    return nodes, links
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

if __name__ == '__main__':
    test()