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
from PowerFlow.data_transfer import export_topo, export_pf_res
from flask import Flask, render_template
from pyecharts import options as opts
from pyecharts.charts import Graph


app = Flask(__name__, static_folder="templates")


def make_graph():
    pf = PowerFlow()
    nodes, links = export_topo(pf)
    c = (
        Graph()
        .add("", nodes, links, repulsion=4000)
        .set_global_opts(title_opts=opts.TitleOpts(title="Graph-GraphNode-GraphLink"))
    )

    return c


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/barChart")
def get_bar_chart():
    c = make_graph()
    return c.dump_options_with_quotes()


if __name__ == "__main__":
    app.run()
