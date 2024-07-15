import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt 
import argparse
from matplotlib.ticker import FuncFormatter

def log_2_formatter(x, pos):
    return '$2^{{{}}}$'.format(int(x))

def log_10_formatter(x, pos):
    return '$10^{{{}}}$'.format(int(x))

result_folder = "./privGraph/results"
graph_folder = "./privGraph/graph"
target_methods = ["privGraph", "edgelist", "adjmat"]
cmap = plt.get_cmap("tab10")
format_dict = {
    "privGraph": {
        "linestyle": "-",
        "marker": "^",
        "color": cmap(3),
        "hatch": ".",
    },
    "edgelist": {
        "linestyle": "-.",
        "marker": "s",
        "color": cmap(0),
        "hatch": 'o',
    },
    "adjmat": {
        "linestyle": "--",
        "marker": "o",
        "color": cmap(1),
        "hatch": "/"
    }
}

query_dict = {
        "EdgeExistQuery": {
            "df_key": "EdgeExistQuery_mean",
            "label": "EdgeExist (ms)",
        },
        "OuttingEdgesCountQuery": {
            "df_key": "OuttingEdgesCountQuery_mean",
            "label": "NeighborsCount (ms)",
        },
        "NeighborsGetQuery": {
            "df_key": "NeighborsGetQuery_mean",
            "label": "NeighborsGet (ms)",   
        }
    }

df_dict = {method: None for method in target_methods}

for method in target_methods:
    target_file = result_folder + "/" + method + "_record.xlsx"
    df = pd.read_excel(target_file)
    df = df.fillna(method="ffill")
    
    if(method == "privGraph"):
        sub_df_privGraph = df
        # sub_df_privGraph = sub_df_privGraph.loc[((sub_df_privGraph["n"] != 32768) | ((sub_df_privGraph["k"] == 16) | (sub_df_privGraph["k"] == 32)))]
        df = sub_df_privGraph
    
    df_dict[method] = df


def query_performance_plot(query_key, df_dict, save_path, upper_legend_flag=False, lower_label_flag=False):
    
    plt.rcParams['font.size'] = 12

    gtypes = df_dict["privGraph"]["gtype"].unique()
    fig, axs= plt.subplots(1, len(gtypes), figsize=(15, 2))
    plt.subplots_adjust(wspace=0.3)

    for i in range(len(gtypes)):
        gtype = gtypes[i]
        ax = axs[i]
        ax2 = ax.twinx()
        n = df_dict["privGraph"].loc[df_dict["privGraph"]["gtype"] == gtype]["n"]
        e = df_dict["privGraph"].loc[df_dict["privGraph"]["gtype"] == gtype]["e"]
        
        query_df_key = query_dict[query_key]["df_key"]
        # ex_time_list = target_df[query_df_key]
        privGraph_time = df_dict["privGraph"].loc[df_dict["privGraph"]["gtype"] == gtype][query_df_key]
        for method in ["edgelist", "adjmat"]:
            base_time = df_dict[method].loc[df_dict[method]["gtype"] == gtype][query_df_key]
            # hatch = format_dict[method]["hatch"]
            offset = -0.25
            if(method == "adjmat"):
                offset = 0.25
            label = "(" + method + "/ privGraph)" + "$\\times$"
            valid_size = len(base_time)
            valid_n = df_dict[method].loc[df_dict[method]["gtype"] == gtype]["n"]

            speedup = base_time.to_numpy() / (privGraph_time[:valid_size]).to_numpy()
            
            if(i == 0):
                ax2.bar(np.log2(valid_n) + offset, np.log2(speedup), width=0.5, label=label, alpha=0.5)
            else:
                ax2.bar(np.log2(valid_n) + offset, np.log2(speedup), width=0.5, alpha=0.5)
            
            
            log2formatter = FuncFormatter(log_2_formatter)
            ax2.yaxis.set_major_formatter(log2formatter)
            ax2.xaxis.set_major_formatter(log2formatter)
            ax2.set_xlabel("|V|")
            
            if(i == (len(gtypes) - 1)):
                ax2.set_ylabel("SpeedUp ($\\times$)")
        
        for method in target_methods:
            linestyle = format_dict[method]["linestyle"]
            marker = format_dict[method]["marker"]
            color = format_dict[method]["color"]
            target_df = df_dict[method].loc[df_dict[method]["gtype"] == gtype]
            valid_n = target_df["n"]
            query_df_key = query_dict[query_key]["df_key"]
            ex_time_list = target_df[query_df_key]
            if(i == 0):
                ax.plot(np.log2(valid_n), np.log2(ex_time_list), label=method, linestyle=linestyle, marker = marker, color = color, linewidth=2.3, markerfacecolor='none', markersize=8)
            else:
                ax.plot(np.log2(valid_n), np.log2(ex_time_list), linestyle=linestyle, marker = marker, color = color, linewidth=2.3, markerfacecolor='none', markersize=8)
                
            log2formatter = FuncFormatter(log_2_formatter)
            ax.yaxis.set_major_formatter(log2formatter)
            if(i == 0):
                ax.set_ylabel(query_dict[query_key]["label"])
                
        ax.set_title(gtype)
    
    if(upper_legend_flag):
        fig.legend(bbox_to_anchor=(0.5, 1.03), loc='lower center', ncol=len(target_methods) + 2)
    
    if(lower_label_flag):
        fig.text(0.5, -0.05, 'Vertex Size $|V|$', ha='center', va='center', fontsize=13)
    
    plt.savefig(save_path, bbox_inches='tight')


if __name__ == "__main__":
    
    # get the test configs. 
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="performance", help="analysis_target")
    args = parser.parse_args()
    target = args.target
    
    # basic performance analysis.
    if target == "performance":
        query_performance_plot("EdgeExistQuery", df_dict, "./privGraph/graph/EdgeExistQuery.png", upper_legend_flag=True, lower_label_flag=True)
        query_performance_plot("OuttingEdgesCountQuery", df_dict, "./privGraph/graph/OuttingEdgesCountQuery.png", upper_legend_flag=False, lower_label_flag=False)
        query_performance_plot("NeighborsGetQuery", df_dict, "./privGraph/graph/NeighborsGetQuery.png", upper_legend_flag=False, lower_label_flag=False)
    