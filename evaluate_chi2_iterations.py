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


def GOF(HistList1,HistList2):
  assert len(HistList1) == len(HistList2)
  wchi2 = 0.
  ndof = 0.
  wchi2_hist2unc = 0.

  for (hist1,hist2) in zip(HistList1,HistList2):
    assert hist1.GetNbinsX() == hist2.GetNbinsX()
    wchi2 += sum([(hist1.GetBinContent(ibin+1)-hist2.GetBinContent(ibin+1))**2/(hist1.GetBinError(ibin+1)**2+hist2.GetBinError(ibin+1)**2) for ibin in range(hist1.GetNbinsX())])
    wchi2_hist2unc += sum([(hist1.GetBinContent(ibin+1)-hist2.GetBinContent(ibin+1))**2/(hist2.GetBinError(ibin+1)**2) for ibin in range(hist1.GetNbinsX())])
    ndof += hist1.GetNbinsX()

  ndof -= 1
  p = CHI2.sf(wchi2,ndof)
  p_hist2unc = CHI2.sf(wchi2_hist2unc,ndof)
  return wchi2,wchi2_hist2unc,ndof,wchi2/ndof,wchi2_hist2unc/ndof,p,p_hist2unc

def GetRefoldIter(File):
  names = [key.GetName() for key in File.GetListOfKeys() if ("HistRecoInclusive_MC" in key.GetName() and "iter" in key.GetName())]

  #extract heads and remove duplicates
  heads= list(dict.fromkeys([name.split("iter")[0] for name in names]))
  print(heads)
  heads_sort =sorted(heads,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  print(heads_sort)

  if any("iter0" in s for s in names):
    start_iter = 0
  else:
    start_iter = 1
  names_rank = [[head+"iter"+str(i) for head in heads_sort] for i in range(int(start_iter),len(names)//len(heads_sort)+start_iter)]
  print(names_rank)
  names_data = [key.GetName() for key in File.GetListOfKeys() if "HistRecoInclusive_Pseudodata" in key.GetName()]

  if len(names_data)>0:
    names_data_sort =sorted(names_data,key = lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  else:
    names_data = [key.GetName() for key in File.GetListOfKeys() if "HistRecoInclusive_Data" in key.GetName()]
    names_data_sort =sorted(names_data,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])

  names_MC = [key.GetName() for key in File.GetListOfKeys() if "HistRecoInclusive_MC_" in key.GetName() and "iter" not in key.GetName()]
  names_MC_sort =sorted(names_MC,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])

  return names_rank,names_data_sort,names_MC_sort

def GetUnfoldIter(File):

  names = [key.GetName() for key in File.GetListOfKeys() if ("HistGenInclusive_MC" in key.GetName() and "iter" in key.GetName())]
  heads = list(dict.fromkeys([name.split("iter")[0] for name in names]))
  heads_sort = sorted(heads,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])

  if any("iter0" in s for s in names):
    start_iter = 0
  else:
    start_iter = 1

  names_rank = [[head+"iter"+str(i) for head in heads_sort] for i in range(int(start_iter),len(names)//len(heads_sort)+start_iter)]
  names_MC = [key.GetName() for key in File.GetListOfKeys() if "HistGenInclusive_MC" in key.GetName() and "iter" not in key.GetName()]
  names_MC_sort =sorted(names_MC,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])

  return names_rank,names_MC_sort

def GetPseudoDataTruth(File):
  names = [key.GetName() for key in File.GetListOfKeys() if "HistGenInclusive_Pseudodata" in key.GetName()]
  names_sort = None
  if len(names)>0:
    names_sort = sorted(names,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  return names_sort

def make_iter_hist( folding, unc, names):
    hist = rt.TH1F(f'{folding}_chi2_{unc}', f'{folding} #chi^{2}/n.d.o.f', len(names), 0, len(names) )
    hist.GetXaxis().SetTitle("iterations")
    hist.GetYaxis().SetTitle('#chi^{2}/n.d.o.f')
    return hist

def prepare_histlist( list_name, hist_name, hist_file):
    hlist = HistList(list_name)
    hlist.root_hists_name = hist_name
    hlist.read_hist_from_file(hist_file)
    hlist.divide_by_bin_width()
    hlist.flatten_hist()

    return hlist

def draw_initial_line( chi2_val, n_iter, hist_name ):
    line = rt.TLine(0, chi2_val, n_iter, chi2_val)
    line.SetLineColor(rt.kRed)
    line.Write(hist_name)

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


def GetEffAcc(f,label,names_gen,names_reco):

  hist_list_geneff = HistList(label+"GenEff")
  hist_list_geneff.root_hists_name=[name_gen.replace("Inclusive","Eff") for name_gen in names_gen]
  if hist_list_geneff.read_hist_from_file(f)==-1:
    hist_list_geneff=None
  else:
    hist_list_geneff.flatten_hist()

  hist_list_recoacc= HistList(label+"RecoAcc")
  hist_list_recoacc.root_hists_name=[name_reco.replace("Inclusive","Acc") for name_reco in names_reco]
  if hist_list_recoacc.read_hist_from_file(f)==-1:
    hist_list_recoacc=None
  else:
    hist_list_recoacc.flatten_hist()

  return {"Eff":hist_list_geneff, "Acc":hist_list_recoacc}


if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--input',default='results_finebin_v7_MCCP1ES_CP5sys_trksys_1d_optimize/unfold_nparticle_eta2p4pt05_pur_1d_spherocity_eta2p4pt05_pur_1d_nominal_optimize_omnifold.root',help='The input root file containing the results in iterations')
    parser.add_argument('--output',default='results_finebin_v7_MCCP1ES_CP5sys_trksys_1d_optimize/iter_nparticle_eta2p4pt05_pur_1d_spherocity_eta2p4pt05_pur_1d_nominal_omnifold.root',help='The output root file containing the refolfing chi^2 w.r.t. iterations')
    parser.add_argument('--config',default="config/plot_1d_v7_MCEPOS_unfoldCP1.json",help="The configration file including the unfolding setup")
    parser.add_argument('--plot',action="store_true",default=False)
    parser.add_argument('--plotdir',default="results_finebin_v7_MCCP1ES_CP5sys_trksys_1d_optimize/plots_optimize")
    parser.add_argument('--sysreweight',action="store_true",default=False)
    parser.add_argument('--plot-software', type=str, choices=["root", "mpl"], default="mpl")
    parser.add_argument('--config-style', type=str, default = 'config/results_style.yml')
    parser.add_argument('--obs', type=obs_pair, help="Pair of observables separated by a comma")
    args = parser.parse_args()

    with open(args.config, 'r') as configjson:
        config = json.load(configjson)

    style_file = args.config_style
    with open(style_file, 'r') as style_json: 
        config_style = yaml.safe_load(style_json)

    (obs1_name, obs1_binning), (obs2_name, obs2_binning) = parse_obs_args( config, args.obs )

    obs1 = ObsConfig.from_yaml( config["varunfold"], [obs1_name], binning=obs1_binning )
    obs2 = ObsConfig.from_yaml( config["varunfold"], [obs2_name], binning=obs2_binning )

    TextListReco = [f'{obs1.reco.edges[i]} #leq {obs1.reco.shortname} < {obs1.reco.edges[i+1]}' for i in range(obs1.reco.nbins)]
    TextListGen = [f'{obs1.gen.edges[i]} #leq {obs1.gen.shortname} < {obs1.gen.edges[i+1]}' for i in range(obs1.gen.nbins)]+(["background"] if config["addbkg"] else [])

    f = rt.TFile(args.input,"READ")

    names_refold, name_data, name_MCreco = GetRefoldIter(f)
    names_unfold, name_MC = GetUnfoldIter(f)
    names_pseudodata_truth = GetPseudoDataTruth(f)

    fout = rt.TFile(args.output,"RECREATE")


    hist_chi2_iter_dataMCunc = make_iter_hist( "Refold", "dataMCunc", names_refold)
    hist_chi2_iter_dataunc = make_iter_hist( "Refold", "dataunc", names_refold)
    hist_chi2_iter_MC = make_iter_hist( "Unfold", "MC", names_unfold)
    hist_chi2_iter_MC_MCunfoldunc = make_iter_hist("Unfold", "MC_MCunfoldunc", names_unfold)

    hist_p_iter_MC_MCunfoldunc = rt.TH1F("Unfold_p_MC_MCunfoldunc", "Unfold p", len(names_unfold), 0, len(names_unfold))
    hist_p_iter_MC_MCunfoldunc.GetXaxis().SetTitle("iterations")
    hist_p_iter_MC_MCunfoldunc.GetYaxis().SetTitle("p-value")

    if not(names_pseudodata_truth is None):

      hist_chi2_iter_truth_dataunc = make_iter_hist("Unfold", "pseudodatatruth_dataunc", names_unfold)
      hist_chi2_iter_truth_dataMCunc = make_iter_hist("Unfold", "pseudodatatruth_dataMCunc", names_unfold)


    _, _, _, wchi2_per_ndof_MCdatareco, wchi2_hist2unc_per_ndof_MCdatareco, p_MCdatareco, _  = GOF([f.Get(name) for name in name_MCreco],[f.Get(name) for name in name_data])
    if not(names_pseudodata_truth is None):
      _, _, _, wchi2_per_ndof_MCdatagen, wchi2_hist2unc_per_ndof_MCdatagen, _, _ = GOF([f.Get(name) for name in name_MC],[f.Get(name) for name in names_pseudodata_truth])

    hist_list_data = prepare_histlist("Data", name_data, f)
    hist_list_MCgeninclusive = prepare_histlist("MCGenInclusive", name_MC, f)
    hist_list_MCrecoinclusive = prepare_histlist("MCRecoInclusive", name_MCreco, f)
    hist_list_MC_eff_acc = GetEffAcc(f, "MC", name_MC,name_MCreco)

    if not(names_pseudodata_truth is None):
      hist_list_pseudodatatruthinclusive = prepare_histlist("PseduodataTruthInclusive", names_pseudodata_truth, f)
      hist_list_pseudodata_eff_acc = GetEffAcc(f, "Pseudodata", names_pseudodata_truth, name_data)

    else:
      hist_list_pseudodatatruthinclusive = None
      hist_list_pseudodata_eff_acc = {"Eff":None, "Acc":None}

    histCfg = {}
    if args.plot:

      data_legend = "Data" if names_pseudodata_truth is None else "Pseudo-data"
      if args.sysreweight:
          data_legend = f"sys variation: {config['syslegend'][0]}/Plot"
      data_color = config_style[args.plot_software]["data_color"]
      MC_color = config_style[args.plot_software]["MC_color"]
      pseudodata_color = config_style[args.plot_software]["pseudodata_color"]

      histCfg["Data"] =            HistConfig(hist_list_data, 0, data_color, "cross", data_legend)
      histCfg["MCGenInclusive"] =  HistConfig(hist_list_MCgeninclusive, 0, MC_color, "fillederror", config["MClegend"])
      histCfg["MCRecoInclusive"] = HistConfig( hist_list_MCrecoinclusive, 0, MC_color, "fillederror", config["MClegend"])

      ps_legend = "Pseudo-data truth" if not args.sysreweight else f"sys variation: {config['pseudodata_name'][0]}" 
      histCfg["PseudodataTruthInclusive"] = HistConfig(hist_list_pseudodatatruthinclusive, 0, pseudodata_color, "triangle", ps_legend)

      histCfg["MCGenEff"] =  HistConfig( hist_list_MC_eff_acc["Eff"], 0, MC_color, "cross", f'{config["MClegend"]} Eff.')
      histCfg["MCRecoAcc"] = HistConfig( hist_list_MC_eff_acc["Acc"], 0, MC_color, "cross", f'{config["MClegend"]} Acc.')

      ps_legend = "Pseudo-data truth" if not args.sysreweight else f"sys variation: {config['pseudodatalegend']} Eff."
      histCfg["PseudodataTruthEff"] = HistConfig( hist_list_pseudodata_eff_acc["Eff"], 0, pseudodata_color, "triangle", ps_legend )

      ps_legend = "Pseudo-data truth" if not args.sysreweight else f"sys variation: {config['pseudodatalegend']} Acc."
      histCfg["PseudodataTruthAcc"] = HistConfig(hist_list_pseudodata_eff_acc["Acc"], 0, pseudodata_color, "triangle", ps_legend )


    pltLists= {}


    print(name_data)
    for i, name_refold in enumerate(names_refold):
      print(name_refold)
      _, _, _, wchi2_per_ndof, wchi2_hist2unc_per_ndof, _, _ = GOF([f.Get(name) for name in name_refold],[f.Get(name) for name in name_data])
      hist_chi2_iter_dataMCunc.SetBinContent(i+1,wchi2_per_ndof)
      hist_chi2_iter_dataMCunc.SetBinError(i+1,0)
      hist_chi2_iter_dataunc.SetBinContent(i+1,wchi2_hist2unc_per_ndof)
      hist_chi2_iter_dataunc.SetBinError(i+1,0)

    for i, name_unfold in enumerate(names_unfold):
      _, _, _, wchi2_per_ndof, wchi2_hist2unc_per_ndof, p, _ = GOF([f.Get(name) for name in name_unfold],[f.Get(name) for name in name_MC])
      iter_index = re.search("iter(\d+)",name_unfold[0]).group(1)
      hist_chi2_iter_MC.SetBinContent(i+1,wchi2_hist2unc_per_ndof)
      hist_chi2_iter_MC.SetBinError(i+1,0)
      hist_chi2_iter_MC_MCunfoldunc.SetBinContent(i+1,wchi2_per_ndof)
      hist_chi2_iter_MC_MCunfoldunc.SetBinError(i+1,0)
      hist_p_iter_MC_MCunfoldunc.SetBinContent(i+1,p)
      hist_p_iter_MC_MCunfoldunc.SetBinError(i+1,0)

      if not(names_pseudodata_truth is None):
        _, _, _, wchi2_per_ndof, wchi2_hist2unc_per_ndof, _, _ = GOF([f.Get(name) for name in name_unfold],[f.Get(name) for name in names_pseudodata_truth])
        hist_chi2_iter_truth_dataunc.SetBinContent(i+1,wchi2_hist2unc_per_ndof)
        hist_chi2_iter_truth_dataunc.SetBinError(i+1,0)
        hist_chi2_iter_truth_dataMCunc.SetBinContent(i+1,wchi2_per_ndof)
        hist_chi2_iter_truth_dataMCunc.SetBinError(i+1,0)


      if args.plot:

        default_rps   = ResultPlotSettings.from_yaml(style_file, [args.plot_software, "MLE"])
        omnifold_rps  = ResultPlotSettings.from_yaml(style_file, [args.plot_software, "omnifold"])
        multifold_rps = ResultPlotSettings.from_yaml(style_file, [args.plot_software, "multifold"])
        unifold_rps   = ResultPlotSettings.from_yaml(style_file, [args.plot_software, "omnifold"])

        result_settings = default_rps
        if "omnifold" in args.input:
          result_settings = omnifold_rps
        elif "multifold" in args.input:
          result_settings = multifold_rps
        elif "unifold" in args.input:
          result_settings = unifold_rps

        result_settings.sys_reweight = args.sysreweight
        color         = result_settings.color
        legend        = result_settings.legend_refold
        color_unfold  = result_settings.color_unfold
        legend_unfold = result_settings.legend_unfold
        tag           = result_settings.tag

        hist_list_refold = prepare_histlist(f"Refold_iter{iter_index}", names_refold[i], f)
        hist_list_unfold = prepare_histlist(f"Unfolder_iter{iter_index}", names_unfold[i], f)

        hist_list_unfold_eff_acc = GetEffAcc(f,"Unfold_iter"+iter_index,names_unfold[i],names_refold[i])


        histCfg["Refold_iter"+iter_index] = HistConfig(hist_list_refold, 0, color, "filled", legend)
        histCfg["Unfold_iter"+iter_index] =  HistConfig(hist_list_unfold, 0, color_unfold, "marker", legend_unfold)
        histCfg[f"Unfold_iter{iter_index}Eff"] = HistConfig( hist_list_unfold_eff_acc["Eff"], 0, color_unfold, "cross", f'{legend_unfold} Eff.')
        histCfg[f"Unfold_iter{iter_index}Acc"] = HistConfig( hist_list_unfold_eff_acc["Acc"], 0, color_unfold, "cross", f'{legend_unfold} Acc.')

        ratio = 'Refold / data' if not args.sysreweight else "Reweight / sys. var."
        pltLists["Refoldcompare_iter"+iter_index] = PlotConfig(
                                                              histCfg["Data"],
                                                              [histCfg[f"Refold_iter{iter_index}"],histCfg["MCRecoInclusive"]],
                                                              f'data_refold_{tag}_iter{iter_index}', 
                                                              ratio) 

        ratio = 'Unfold / MC' if not args.sysreweight else "Reweight / MC"
        pltLists["Unfoldcompare_iter"+iter_index] = PlotConfig( 
                                                              histCfg["MCGenInclusive"], 
                                                              [histCfg[f"Unfold_iter{iter_index}"]], 
                                                              f'MC_unfold_{tag}_iter{iter_index}', 
                                                              ratio )

        ratio = 'Unfold / Truth' if not args.sysreweight else "Reweight / sys. var."
        pltLists["Unfoldcomparepseudodata_iter"+iter_index] = PlotConfig( 
                                                                        histCfg["PseudodataTruthInclusive"], 
                                                                        [histCfg[f"Unfold_iter{iter_index}"],histCfg["MCGenInclusive"]], 
                                                                        f'pseudodata_truth_unfold_{tag}_iter{iter_index}', 
                                                                        ratio) 

        ratio = 'Unfold / MC' if not args.sysreweight else "Reweight / MC"
        pltLists["Unfoldcompareeff_iter"+iter_index]= PlotConfig( 
                                                              histCfg["MCGenEff"],
                                                              [histCfg[f"Unfold_iter{iter_index}Eff"]],
                                                              f'MC_unfoldedeff_{tag}_iter{iter_index}',
                                                              ratio )

        pltLists["Unfoldcompareacc_iter"+iter_index]= PlotConfig( histCfg["MCRecoAcc"],
                                                                [histCfg[f"Unfold_iter{iter_index}Acc"]],
                                                                f'MC_unfoldacc_{tag}_iter{iter_index}',
                                                                ratio )

        ratio = 'Unfold / Truth' if not args.sysreweight else "Reweight / sys. var"
        pltLists["Unfoldcomparepseudodataeff_iter"+iter_index]= PlotConfig( 
                                                              histCfg["PseudodataTruthEff"],
                                                              [histCfg[f"Unfold_iter{iter_index}Eff"], histCfg["MCGenEff"]],
                                                              f'pseudodatatruth_unfoldedeff_{tag}_iter{iter_index}',
                                                              ratio )

        pltLists["Unfoldcomparepseudodataacc_iter"+iter_index]= PlotConfig( histCfg["PseudodataTruthAcc"],
                                                                [histCfg[f"Unfold_iter{iter_index}Acc"],histCfg["MCRecoAcc"]],
                                                                f'pseudodatatruth_unfoldacc_{tag}_iter{iter_index}',
                                                                ratio )

        to_plot = [f'Refoldcompare_iter{iter_index}', f'Unfoldcompare_iter{iter_index}', f'Unfoldcompareeff_iter{iter_index}', f'Unfoldcompareacc_iter{iter_index}' ]
        if names_pseudodata_truth is not None:
            to_plot += [ f'Unfoldcomparepseudodata_iter{iter_index}', f"Unfoldcomparepseudodataeff_iter{iter_index}", f"Unfoldcomparepseudodataacc_iter{iter_index}" ]

        print(f"plotting iteration {iter_index}")
        for plt_type in to_plot:
            txt_list = TextListReco if pltLists[plt_type].is_reco else TextListGen
            draw_plot( pltLists[plt_type], args.plotdir, obs1_name, obs2_name, obs2, txt_list, use_root = (args.plot_software == "root"))


    hist_chi2_iter_dataMCunc.Write()
    hist_chi2_iter_dataunc.Write()
    hist_chi2_iter_MC.Write()
    hist_chi2_iter_MC_MCunfoldunc.Write()

    hist_p_iter_MC_MCunfoldunc.Write()

    line_reco_dataMCunc = draw_initial_line( wchi2_per_ndof_MCdatareco, len(names_refold), "Chi2dataMC_reco_dataMCunc")
    line_reco_dataunc = draw_initial_line( wchi2_hist2unc_per_ndof_MCdatareco, len(names_refold), "Chi2dataMC_reco_dataunc")
    line_reco_dataMCunc_p = draw_initial_line(p_MCdatareco,len(names_refold), "p_dataMC_reco_dataMCunc")

    if not(names_pseudodata_truth is None):
      hist_chi2_iter_truth_dataunc.Write()
      hist_chi2_iter_truth_dataMCunc.Write()

      line_gen_dataunc = draw_initial_line( wchi2_hist2unc_per_ndof_MCdatagen, len(names_unfold), "Chi2pseudodataMC_gen_pseudodataunc")
      line_gen_dataMCunc = draw_initial_line( wchi2_per_ndof_MCdatagen, len(names_unfold), "Chi2pseudodataMC_gen_pseudodataMCunc")

    fout.Close()

