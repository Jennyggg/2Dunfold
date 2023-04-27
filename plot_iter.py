import numpy as np
import json
import yaml
from argparse import ArgumentParser
import os
import ROOT as rt
from Plotting_cfg import *
from scipy.stats import chi2 as CHI2
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from plot_utils import *
from unfold_utils import *
from arg_parsing import *
from plot_mplhep_utils import *
from dataclasses import dataclass, field, asdict
from typing import Union
from plot_configs import ResultPlotSettings, HistConfig, PlotConfig
from configs import ObsConfig
import pandas as pd
from ast import literal_eval

@dataclass
class HistConfigCollection:
  gen: HistConfig = None
  reco: HistConfig = None
  eff: HistConfig = None
  acc: HistConfig = None
  mig: HistConfig = None


@dataclass
class Chi2Collection:
  gen_both_error: HistConfig = None
  gen_both_error_rel: HistConfig = None
  gen_target_error: HistConfig = None
  gen_target_error_rel: HistConfig = None
  reco_both_error: HistConfig = None
  reco_both_error_rel: HistConfig = None
  reco_target_error: HistConfig = None
  reco_target_error_rel: HistConfig = None

@dataclass
class KSCollection:
  gen_distance: HistConfig = None
  gen_distance_rel: HistConfig = None
  reco_distance: HistConfig = None
  reco_distance_rel: HistConfig = None
  gen_p: HistConfig = None
  gen_p_rel: HistConfig = None
  reco_p: HistConfig = None 
  reco_p_rel: HistConfig = None



def get_df_entries(df, **kwargs_selections):
  for ikey,key in enumerate(kwargs_selections.keys()):
    if ikey == 0:
      selection = (df[key]==kwargs_selections[key])
    else:
      selection = selection & (df[key]==kwargs_selections[key])
    
  return df[selection]

def get_histconfigcollection(df,color_list,style_list,legend_list):
  gen_entry = df[(df['histtype']=='inclusive')&df['dim1_isgen']&df['dim2_isgen']]
  reco_entry = df[(df['histtype']=='inclusive')&~df['dim1_isgen']&~df['dim2_isgen']] 
  #Plan: process the results from bootstraps -> add 'histtype' options 'inclusive_sysunc', 'inclusive_statunc','inclusive_totalunc'
  eff_entry = df[(df['histtype']=='efficiency')&df['dim1_isgen']&df['dim2_isgen']]
  acc_entry = df[(df['histtype']=='acceptance')&~df['dim1_isgen']&~df['dim2_isgen']]
  mig_entry = df[(df['histtype']=='migration')&df['dim1_isgen']&~df['dim2_isgen']] 

  if not isinstance(color_list,list):
    color_list = [color_list]*4+[None]

  if not isinstance(style_list,list):
    style_list = [style_list]*4+[None]

  if not isinstance(legend_list,list):
    legend_list = [legend_list,legend_list,legend_list+" Eff.",legend_list+" Acc.",legend_list+" Mig."]

  histarray_statonly = [0]*5 # Replace with the HistArray from stat-unc-only unfolding results

  list_HistConfig = []
  for i,entry in enumerate([gen_entry,reco_entry,eff_entry,acc_entry,mig_entry]):
    if len(entry)>0:
      histarray = HistArray.from_df_entry(entry.iloc[0])
      if i<2:
        histarray.divide_by_bin_width()
      list_HistConfig.append(HistConfig(histarray,
                                        histarray_statonly[i],
                                        color_list[i], style_list[i],legend_list[i]))
    else:
      list_HistConfig.append(None)

  return HistConfigCollection(*list_HistConfig)

def get_histconfig_gof(df,color,legend,histtype_list):
  #gen_both_error_entry = df[(df['histtype']=='gen_both_error')]
  #gen_target_error_entry = df[(df['histtype']=='gen_target_error')]
  #reco_both_error_entry = df[(df['histtype']=='reco_both_error')]
  #reco_target_error_entry = df[(df['histtype']=='reco_target_error')]

  list_gof = []
  #for i,entry in enumerate([gen_both_error_entry,gen_target_error_entry,reco_both_error_entry,reco_target_error_entry]):
  for i,histtype in enumerate(histtype_list):
    entry = df[(df['histtype']==histtype)]
    if len(entry)>0:
      histarray = HistArray.from_df_entry(entry.iloc[0])
      list_gof.append(HistConfig(histarray,0,color,"hist",legend))
      list_gof.append(HistConfig(histarray/histarray.nested_value[0][0],0,color,"hist",legend))
    else:
      list_gof.append(None)
      list_gof.append(None)

  return Chi2Collection(*list_gof),KSCollection(*list_gof)



def plot_wrapper(plt_list,**kwargs):
  input_args = {
    "hist_ref_stat": plt_list.ref.stat,
    "style_ref": plt_list.ref.style,
    "color_ref": plt_list.ref.color,
    "list_style_compare": plt_list.styles,
    "list_color_compare": plt_list.colors,
    "label_ratio":plt_list.ratio
  }
  use_root  = kwargs.pop("use_root")
  input_args.update(kwargs)
  if use_root:
    plot_flat_hists_root(plt_list.ref.hist, plt_list.hists, plt_list.ref.legend, plt_list.legends, **input_args)
  else:
    plot_flat_hists_mpl(plt_list.ref.hist, plt_list.hists, plt_list.ref.legend, plt_list.legends, **input_args)

def draw_plot(plt_list, plotdir, var1_nm, var2_nm, obs2, txt_list,use_root=True):
    if plt_list.ref is None: return
    path = f'{plotdir}/{plt_list.name}_{var1_nm}_{var2_nm}'
    os.system(f"mkdir -p {plotdir}")
    axis_title = obs2.reco.name if plt_list.is_reco else obs2.gen.name

    plot_args = {
      "title": axis_title,
      "is_logY": 0,
      "do_ratio": 1,
      "output_path": path,
      "text_list": txt_list,
      "labelY":'Normalized Events/Bin Width',
      "range_ratio": 0.1 if "eff" in plt_list.name or "acc" in plt_list.name else 0.3,
      "use_root": use_root
    }
    plot_wrapper(plt_list,**plot_args)

    if not("eff" in plt_list.name or "acc" in plt_list.name):
      plot_args["is_logY"] = 1
      plot_wrapper(plt_list,**plot_args)

def draw_gof(plt_list, plotdir):
    if plt_list.ref is None: return
    path = f'{plotdir}/{plt_list.name}'
    os.system(f"mkdir -p {plotdir}")
    axis_title = "Iteration"

    plot_args = {
      "title": axis_title,
      "is_logY": 0,
      "do_ratio": 0,
      "output_path": path,
      "text_list": [""],
      "labelY":'Goodness of fit',
      "range_ratio": 0.,
      "use_root": 0
    }
    plot_wrapper(plt_list,**plot_args)


def get_sel_legend_color(config,config_style,datatype):
  sel = {'datatype': datatype,
         'dataset': config[f"{datatype}_name"]}
  legend = datatype + ": " + config[f"{datatype}_name"]
  color = config_style['mpl'][f"{datatype}_color"]
  return sel,legend,color

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--config',default="config/plot_1d_v7_MCEPOS_unfoldCP1.json",help="The configration file including the unfolding setup")
    parser.add_argument('--config-style', type=str, default = 'config/results_style.yml')
    parser.add_argument('--obs', type=obs_pair,nargs = '+', help="Pair of observables separated by a comma")
    args = parser.parse_args()

    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    with open(args.config, 'r') as configjson:
        config = json.load(configjson)

    style_file = args.config_style
    with open(style_file, 'r') as style_json:
        config_style = yaml.safe_load(style_json)

    df = pd.read_csv(config["input"])

    for method in config["methods"]:
      result_settings = ResultPlotSettings.from_yaml(style_file, ['mpl', method])
      result_settings.gen_reweight = (config['workflow']=='genreweight')
      result_settings.sys_reweight = (config['workflow']=='sysreweight')

      chi2_collection_list = []
      ks_collection_list = []
      for icolor,obs in enumerate(args.obs):
      
        (obs1_name, obs1_binning), (obs2_name, obs2_binning) = parse_obs_args( config, obs )

        obs1 = ObsConfig.from_yaml( config["varunfold"], [obs1_name], binning=obs1_binning )
        obs2 = ObsConfig.from_yaml( config["varunfold"], [obs2_name], binning=obs2_binning )
        TextListReco = [f'{obs1.reco.edges[i]} #leq {obs1.reco.shortname} < {obs1.reco.edges[i+1]}' for i in range(obs1.reco.nbins)]
        TextListGen = [f'{obs1.gen.edges[i]} #leq {obs1.gen.shortname} < {obs1.gen.edges[i+1]}' for i in range(obs1.gen.nbins)]+(["background"] if config["addbkg"] else [])

        base_sel = {'obs1': obs[0],
                    'obs2': obs[1],
                    'method': method
                   }
        if "pseudodata" in config.keys() and config["pseudodata"]:
          data_sel, data_legend, data_color = get_sel_legend_color(config,config_style,"pseudodata")
          chi2_sel = {'datatype': f"chi2_{config['workflow']}_pseudodata",
                      'dataset':f"{config['workflow']}_{config['pseudodata_name']}"}
          ks_sel = {'datatype': f"ks_{config['workflow']}_pseudodata",
                      'dataset':f"{config['workflow']}_{config['pseudodata_name']}"}
        else:
          data_sel, data_legend, data_color = get_sel_legend_color(config,config_style,"data")
          chi2_sel = {'datatype': f"chi2_{config['workflow']}_data",
                      'dataset':f"{config['workflow']}_{config['data_name']}"}
          ks_sel = {'datatype': f"ks_{config['workflow']}_data",
                      'dataset':f"{config['workflow']}_{config['data_name']}"}


        MC_sel, MC_legend, MC_color = get_sel_legend_color(config,config_style,"MC")



        records_data = get_df_entries(df, **base_sel|data_sel)
        records_MC = get_df_entries(df, **base_sel|MC_sel)
        records_chi2 = get_df_entries(df, **base_sel|chi2_sel)
        records_ks = get_df_entries(df, **base_sel|ks_sel)


        hcc_data = get_histconfigcollection(records_data,data_color,['triangle','cross','triangle','triangle',None],data_legend)

        hcc_MC = get_histconfigcollection(records_MC,MC_color,['fillederror','fillederror','cross','cross',None],MC_legend)

        chi2_collection,_ = get_histconfig_gof(records_chi2,default_colors[icolor],f"{obs[0]}:{obs[1]}",["gen_both_error","gen_target_error","reco_both_error","reco_target_error"])
        chi2_collection_list.append(chi2_collection)

        _,ks_collection = get_histconfig_gof(records_ks,default_colors[icolor],f"{obs[0]}:{obs[1]}",["gen_distance","reco_distance","gen_p","reco_p"])
        ks_collection_list.append(ks_collection)

        iters = np.sort(np.unique(df[~np.isnan(df['iter'])]['iter']))

        for it in iters:
          iter_sel = {'datatype': config['workflow'],
                      'dataset': config['MC_name'],
                      'iter': it}
          records_it = get_df_entries(df, **base_sel|iter_sel)
          it_color = result_settings.color
          it_color_unfold = result_settings.color_unfold
          it_legend = result_settings.legend_refold
          it_legend_unfold = result_settings.legend_unfold
          tag = result_settings.tag
          hcc_it = get_histconfigcollection(records_it,[it_color_unfold,it_color,it_color_unfold,it_color_unfold,None],["marker","filled","cross","cross",None],[it_legend_unfold,it_legend,it_legend_unfold+" Eff.",it_legend_unfold+" Acc.",it_legend_unfold+" Mig."])

          if method in config['evaluate-sys-bootstrap']:
            iter_sysbs_sel = {'datatype': config['workflow']+"_sys",
                      'dataset': config['MC_name'],
                      'iter': it}
            records_it_sysbs = get_df_entries(df, **base_sel|iter_sysbs_sel)
            it_sysbs_color = result_settings.color
            it_sysbs_color_unfold = result_settings.color_unfold
            it_sysbs_legend = result_settings.legend_refold+" sys."
            it_sysbs_legend_unfold = result_settings.legend_unfold+" sys."
            hcc_it_sysbs = get_histconfigcollection(records_it_sysbs,[it_sysbs_color_unfold,it_sysbs_color,it_sysbs_color_unfold,it_sysbs_color_unfold,None],["hatch","hatch","hatch","hatch",None],[it_sysbs_legend_unfold,it_sysbs_legend,it_sysbs_legend_unfold+" Eff.",it_sysbs_legend_unfold+" Acc.",it_sysbs_legend_unfold+" Mig."])
          else:
            hcc_it_sysbs = None

          pltLists = {}
          ratio = 'Refold / data' if config['workflow']=="unfold" else "Reweight / sys. var."
          pltLists[f"Refoldcompare_iter{it}"] = PlotConfig(
                                                           hcc_data.reco,
                                                           [hcc_it.reco,hcc_it_sysbs.reco,hcc_MC.reco],
                                                           f'data_refold_{tag}_iter{it}',
                                                           ratio)

          ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          pltLists[f"Unfoldcompare_iter{it}"] = PlotConfig(
                                                           hcc_MC.gen,
                                                           [hcc_it.gen,hcc_it_sysbs.gen],
                                                           f'MC_unfold_{tag}_iter{it}',
                                                           ratio )

          ratio = 'Unfold / Truth' if config['workflow']=="unfold" else "Reweight / sys. var."
          pltLists[f"Unfoldcomparepseudodata_iter{it}"] = PlotConfig(
                                                                     hcc_data.gen,
                                                                     [hcc_it.gen,hcc_it_sysbs.gen,hcc_MC.gen],
                                                                     f'pseudodata_truth_unfold_{tag}_iter{it}',
                                                                     ratio)

          ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          pltLists[f"Unfoldcompareeff_iter{it}"]= PlotConfig(
                                                              hcc_MC.eff,
                                                              [hcc_it.eff,hcc_it_sysbs.eff],
                                                              f'MC_unfoldedeff_{tag}_iter{it}',
                                                              ratio )

          pltLists[f"Unfoldcompareacc_iter{it}"]= PlotConfig( hcc_MC.acc,
                                                                [hcc_it.acc,hcc_it_sysbs.acc],
                                                                f'MC_unfoldacc_{tag}_iter{it}',
                                                                ratio )

          ratio = 'Unfold / Truth' if config['workflow']=="unfold" else "Reweight / sys. var"
          pltLists[f"Unfoldcomparepseudodataeff_iter{it}"]= PlotConfig(
                                                              hcc_data.eff,
                                                              [hcc_it.eff,hcc_it_sysbs.eff, hcc_MC.eff],
                                                              f'pseudodatatruth_unfoldedeff_{tag}_iter{it}',
                                                              ratio )

          pltLists[f"Unfoldcomparepseudodataacc_iter{it}"]= PlotConfig( hcc_data.acc,
                                                                [hcc_it.acc,hcc_it_sysbs.acc,hcc_MC.acc],
                                                                f'pseudodatatruth_unfoldacc_{tag}_iter{it}',
                                                                ratio )
          to_plot = [f'Refoldcompare_iter{it}', f'Unfoldcompare_iter{it}', f'Unfoldcompareeff_iter{it}', f'Unfoldcompareacc_iter{it}',f'Unfoldcomparepseudodata_iter{it}', f"Unfoldcomparepseudodataeff_iter{it}", f"Unfoldcomparepseudodataacc_iter{it}" ]

          print(f"plotting iteration {it}")
          for plt_type in to_plot:
            txt_list = TextListReco if pltLists[plt_type].is_reco else TextListGen
            draw_plot( pltLists[plt_type], config["outputdir"], obs1_name, obs2_name, obs2, txt_list, use_root = False)
        chi2_plot = ["gen_target_error","gen_target_error_rel","reco_target_error","reco_target_error_rel"]
        ks_plot = ["gen_distance","gen_distance_rel","reco_distance","reco_distance_rel","gen_p","gen_p_rel","reco_p","reco_p_rel"]

        for (gof,gof_plot,gof_collection_list) in zip(["chi2","ks"],[chi2_plot,ks_plot],[chi2_collection_list,ks_collection_list]):
          for name in gof_plot:
            gof_list = [getattr(gof_coll,name) for gof_coll in gof_collection_list]
            gof_list = [histconfig for histconfig in gof_list if histconfig is not None]
            if len(gof_list)>0:
              plotconfig = PlotConfig( gof_list[0],
                                       gof_list[1:] if len(gof_list)>1 else [],
                                       f"{gof}_{name}_{method}",
                                       "" )
              draw_gof( plotconfig, config["outputdir"] )




