# -*- coding: UTF-8 -*-
'''
Name:         data_transfer.py
Func:         To transfer the topology or power flow result 
              to front end
Author:       4thFever
Addr:         Beijing, China
Time:         2020-06-16 
Link:         https://github.com/4thfever
'''
from PowerFlow.power_flow import PowerFlow
from pyecharts import options as opts
from flask import url_for
import numpy as np

# 向前端发送拓扑信息
def export_topo(pf):
    links, nodes = [], []
    for num in range(1, pf.num_node+1):
        node_value = [list(map(
                      lambda x:round(x,3),
                      pf.df_node.loc[num, ['Um', 'Ua']].values
                    ))]
        node = opts.GraphNode(
                    name=f"Node{num}", 
                    symbol_size=20,
                    value=node_value
                    ) 
        nodes.append(node)

    count_load, count_tran, count_line = 1, 1, 1
    for ele in pf.list_elements:
        node_name = None
        # 默认value为空格，这样子在前端不显示
        # 如果不赋值，会显示undefined
        node_value, edge_value = " ", " "
        if ele.type_ == 'load':
            node_name=f"Lode{count_load}"
            node_symbol='diamond'
            # 用list形式存link信息，因为为tran时需两个link
            sources=[f"Lode{count_load}"]
            targets=[f"Node{ele.i}"]
            count_load += 1
        elif ele.type_ == 'gene':
            # 标记节点
            if ele.j == -1:
                value = 'PV节点'
            if ele.j == 0:
                value = '平衡节点'
            node_name=f"Gene{ele.j}"
            node_symbol='rect'
            sources=[f"Gene{ele.j}"]
            targets=[f"Node{ele.i}"]
        elif ele.type_ == 'tran':
            node_name=f"Tran{count_tran}"
            node_symbol='triangle'
            sources=[f"Tran{count_tran}"]*2
            targets=[f"Node{ele.i}", f"Node{ele.j}"]
            count_tran += 1
        elif ele.type_ == 'line':
            edge_value=[list(map(
                      lambda x:round(x,3),
                      pf.df_branch.loc[count_line,['dP', 'dQ']].values
                    ))]
            sources=[f"Node{ele.i}"] 
            targets=[f"Node{ele.j}"]
            count_line += 1
        if isinstance(node_name, str):
            node = opts.GraphNode(
                    name = node_name,
                    value = node_value,
                    symbol = node_symbol,
                    symbol_size = 20
                )
            nodes.append(node)
        for source, target in zip(sources, targets):
            link = opts.GraphLink(
                source=source,
                target=target,
                value=edge_value
            )
            links.append(link)
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