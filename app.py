# -*- coding: UTF-8 -*-
'''
Name:         app.py
Func:         The frontend showing power grid's topology
              and powerflow based on pyecharts
Author:       4thFever
Addr:         Beijing, China
Time:         2020-06-20
Link:         https://github.com/4thfever
'''
from PowerFlow.power_flow import PowerFlow
from data_transfer import export_topo, export_pf_res
from flask import Flask, render_template, redirect, url_for
from pyecharts import options as opts
from pyecharts.charts import Graph


app = Flask(__name__, static_folder="templates")


def make_graph(pf_data=False):
    pf = PowerFlow()
    pf.run()
    nodes, links = export_topo(pf)
    c = (
        Graph()
        .add("", 
            nodes, 
            links,
            edge_label=opts.LabelOpts(
                    is_show=True, 
                    position="middle", 
                    formatter="{c}"
            ), 
            label_opts=opts.LabelOpts(
                    is_show=True, 
                    position="middle", 
                    formatter="{c}"
            ), 
            repulsion=4000
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(),
            toolbox_opts=opts.ToolboxOpts(),
            legend_opts=opts.LegendOpts(
                is_show=True,
                legend_icon = 'rect',
            )
        )
    )
    return c


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/chart")
def get_chart():
    c = make_graph()
    return c.dump_options_with_quotes()

# 通过更新数据的方式，给拓扑形状输入具体的潮流值
# 这样做是为了之后动态训练的时候能够更新潮流情况
@app.route("/updateData")
def update_data():
    c = make_graph(pf_data=True)
    return c.dump_options_with_quotes()


if __name__ == "__main__":
    app.run(debug=True)
