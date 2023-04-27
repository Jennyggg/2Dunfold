import numpy as np
import json
import yaml
from argparse import ArgumentParser
import os
import ROOT
import re
from unfold_utils import *
from arg_parsing import *
from configs import ObsConfig
import pandas as pd
from GOF.binned import *
import uproot
import awkward as ak
import subprocess
from multiprocessing.pool import Pool

def get_metadata_config(dataset,tag="",**dic_config):
  workflow = dic_config.pop("workflow")
  MC_name = dic_config.pop("MC_name")
  data_name = dic_config.pop("data_name")
  pseudodata_name = dic_config.pop("pseudodata_name")
  from_step1 = dic_config["from_step1"]
  if 'chi2' in dataset:
    gof = 'chi2'
  elif 'ks' in dataset:
    gof = 'ks'
  else:
    gof = ''
  if dataset == "MC":
    datatype = "MC"
    df_dataset = MC_name
    from_step1 = False
  elif dataset == "Data":
    datatype = "data"
    df_dataset = data_name
    from_step1 = False
  elif dataset == "Pseudodata":
    datatype = "pseudodata"
    df_dataset = pseudodata_name
    from_step1 = False
  elif dataset == f"{gof}_MC_data":
    datatype = f"{gof}_MC_data"
    df_dataset = MC_name+"_"+data_name
    from_step1 = False
  elif dataset == f"{gof}_MC_pseudodata":
    datatype = f"{gof}_MC_pseudodata"
    df_dataset = MC_name+"_"+pseudodata_name
    from_step1 = False
  elif dataset == f"{gof}_unfold_pseudodata":
    datatype = f"{gof}_"+workflow+"_pseudodata"
    df_dataset = workflow+"_"+pseudodata_name
  elif dataset == f"{gof}_unfold_data":
    datatype = f"{gof}_"+workflow+"_data"
    df_dataset = workflow+"_"+data_name
  elif "sysbs" in dataset:
    ibs = re.findall(r'\d+', dataset)[-1]
    datatype = workflow+"_sysbs_"+ibs
    print(datatype)
    df_dataset = MC_name
  else:
    datatype = workflow
    df_dataset = MC_name
  if "iter" in tag:
    i_iter = int(re.findall(r'\d+',tag)[0])
  else:
    i_iter = np.nan
  dic_config["datatype"] = datatype
  dic_config["dataset"] = df_dataset
  dic_config["iter"] = i_iter
  dic_config["from_step1"] = from_step1
  return dic_config



def fill_hist_lists(dataset,o1,o2,edges_gen,edges_reco,source,genWeight="",from_root=True,weight_array=None,store_mig=False,tag="", reco_only=False,**kwargs):
  hists = {}

  metadata_config = get_metadata_config(dataset,tag,**kwargs)

  reco_inclusive_args = {"histtype":"inclusive"}
  reco_inclusive_args.update(metadata_config)
  reco_inclusive=HistList("HistRecoInclusive_"+dataset,tag,**reco_inclusive_args)
  reco_inclusive.read_settings_from_config_dim1(o1,isgen=False)
  reco_inclusive.read_settings_from_config_dim2(o2,isgen=False)
  reco_inclusive.bin_edges_dim2 = edges_reco
  reco_inclusive.cut = cuts[CutType.PassReco]
  hists["reco_inclusive"] = reco_inclusive

  if not reco_only:
    gen_passreco_args = {"histtype":"pass_gen_reco"}
    gen_passreco_args.update(metadata_config)
    gen_passreco=HistList("HistGen_"+dataset,tag,**gen_passreco_args)
    gen_passreco.read_settings_from_config_dim1(o1,isgen=True)
    gen_passreco.read_settings_from_config_dim2(o2,isgen=True)
    gen_passreco.bin_edges_dim2 = edges_gen
    gen_passreco.cut = cuts[CutType.PassReco_PassGen]
    hists["gen_passreco"] = gen_passreco

    gen_inclusive_args = {"histtype":"inclusive"}
    gen_inclusive_args.update(metadata_config)
    gen_inclusive=HistList("HistGenInclusive_"+dataset,tag,**gen_inclusive_args)
    gen_inclusive.read_settings_from_config_dim1(o1,isgen=True)
    gen_inclusive.read_settings_from_config_dim2(o2,isgen=True)
    gen_inclusive.bin_edges_dim2 = edges_gen
    gen_inclusive.cut = cuts[CutType.PassGen]
    hists["gen_inclusive"] = gen_inclusive

    reco_passgen_args = {"histtype":"pass_gen_reco"}
    reco_passgen_args.update(metadata_config)
    reco_passgen=HistList("HistReco_"+dataset,tag,**reco_passgen_args)
    reco_passgen.read_settings_from_config_dim1(o1,isgen=False)
    reco_passgen.read_settings_from_config_dim2(o2,isgen=False)
    reco_passgen.bin_edges_dim2 = edges_reco
    reco_passgen.cut = cuts[CutType.PassReco_PassGen]
    hists["reco_passgen"] = reco_passgen

    gen_eff_args = {"histtype":"efficiency"}
    gen_eff_args.update(metadata_config)
    gen_eff = HistList("HistGenEff_"+dataset,tag,**gen_eff_args)
    gen_eff.read_settings_from_config_dim1(o1,isgen=True)
    gen_eff.read_settings_from_config_dim2(o2,isgen=True)
    gen_eff.bin_edges_dim2 = edges_gen
    gen_eff.fill_root_hists_name()

    reco_acc_args = {"histtype":"acceptance"}
    reco_acc_args.update(metadata_config)
    reco_acc = HistList("HistRecoAcc_"+dataset,tag,**reco_acc_args)
    reco_acc.read_settings_from_config_dim1(o1,isgen=False)
    reco_acc.read_settings_from_config_dim2(o2,isgen=False)
    reco_acc.bin_edges_dim2 = edges_reco
    reco_acc.fill_root_hists_name()


  for _, hist in hists.items():
    hist.fill_root_hists_name()
    print("histograms:", hist.root_hists_name)
    hist.fill_hist(source, from_root, weightarray=weight_array,genWeight=genWeight)
  print("filled")

  if not reco_only:
    gen_eff.get_hist_from_division(gen_passreco,gen_inclusive)
    reco_acc.get_hist_from_division(reco_passgen,reco_inclusive)
    hists["eff"] = gen_eff
    hists["acc"] = reco_acc

  mig = None
  if store_mig:
    mig = []
    for i, edge_i in enumerate(o1.gen.edges[:-1]):
      mig.append([])
      for j, edge_j in enumerate(o1.reco.edges[:-1]):
        edge_ip1 = o1.gen.edges[i+1]
        edge_jp1 = o1.reco.edges[j+1]

        #mig_ij
        mig_args = {"histtype":"migration",
                    "mig_o1_index":{"gen": i,
                                    "reco": j},
                    "mig_o1_range":{"gen":[edge_i,edge_ip1],
                                    "reco": [edge_j,edge_jp1]}
                   }
        mig_args.update(metadata_config)
        this_mig = HistList(f'HistMig_{dataset}_{o1.gen.np_var}{edge_i}-{edge_ip1}_{o1.reco.np_var}{edge_j}-{edge_jp1}{tag}',**mig_args)
        mig[i].append(this_mig)
        mig[i][j].read_settings_from_config_dim1(o2,isgen=True)
        mig[i][j].read_settings_from_config_dim2(o2,isgen=False)

        gen_cuts_root = f" ({o1.gen.root_var}>={edge_i}) * ({o1.gen.root_var}<{edge_ip1}) "
        reco_cuts_root = f" ({o1.reco.root_var}>={edge_j}) * ({o1.reco.root_var}<{edge_jp1}) "
        all_cuts_root = root_cuts[CutType.PassReco_PassGen] + "*" + gen_cuts_root + "*" + reco_cuts_root
        mig[i][j].root_cut = all_cuts_root

        gen_cuts_np = [[o1.gen.np_var,">=",str(edge_i)], [o1.gen.np_var,"<",str(edge_ip1)]]
        reco_cuts_np = [[o1.reco.np_var, ">=", str(edge_j)], [o1.reco.np_var,"<",str(edge_jp1)]]
        all_cuts_np = np_cuts[CutType.PassReco_PassGen] + gen_cuts_np + reco_cuts_np
        mig[i][j].npy_cut = all_cuts_np

        mig[i][j].fill_2Dhist(source, from_root=from_root, weightarray=weight_array, genWeight=genWeight)
        print("filled histogram:",mig[i][j].name)

    for i in range(len(o1.gen.edges[:-1])):
      bin_norm_gen = np.sum([np.array(hist.bin_norm) for hist in mig[i]], axis = 0 )
      underflow = np.sum([hist.dim1_underflow])
      overflow = np.sum([hist.dim1_overflow])
      for j in range(len(o1.reco.edges[:-1])):
        mig[i][j].binwise_multiply_2Dhist(scalar_list = np.nan_to_num(1./bin_norm_gen),scalar_underflow=np.nan_to_num(1./underflow), scalar_overflow=np.nan_to_num(1./overflow))
  hists["mig"] = mig
  return hists

def fill_gof_histdata(chi2collections,**kwargs):
  dataset = chi2collections.name
  metadata_config = get_metadata_config(dataset,'',**kwargs)
  return HistData.diclist_from_gofcoll(chi2collections,**metadata_config)


def get_input_file_from_config_info( input_val, var ): 
    if isinstance(input_val, dict):
        if var in input_val.keys():
            return input_val[var]
        else:
            return input_val["default"]

    return input_val

def get_sys_variations( config, obs_nm ):
    tree_sys_list = []
    for sys in config["inputfilesim_sys"].keys():
      for vari in ["up", "down"]:
        inputfile_sys = get_input_file_from_config_info( config["inputfilesim_sys"][sys][vari], obs_nm )
        fin_sys_vari = ROOT.TFile(inputfile_sys)
        default = fin_sys_vari.Get("ntuplizer/tree")
        if default:
            tree = default
        else:
            tree = fin_sys_vari.Get("tree")
        if tree is not None:
            tree_sys_list.append(tree)
    return tree_sys_list

def get_bin_edges( conf, var1, var2, trees):

    cut = root_cuts[CutType.PassReco_PassGen]
    bin_edges = []
    for lvl in ["reco","gen"]:
        do_merge = conf[f"merge{lvl}bin"]
        if do_merge:
            edges = merge_bins(obs=[var1[lvl].root_var , var2[lvl].root_var], trees=trees, root_cut=cut, threshold=conf[f"mergethreshold{lvl}"], bin_edges_dim1_1d=var1[lvl].edges, bin_edges_dim2_1d=var2[lvl].edges)
        else:
            edges = [var2[lvl].edges * var1[lvl].nbins ]
        bin_edges.append(edges)

    return bin_edges


def write_all_hists( hist_dict ):
    for key, hist in hist_dict.items():
        if hist == None:
            continue

        if key == 'mig':
          for column in hist:
            for row in column:
              row.write_2Dhist()
        else:
          hist.write_hist_list()

def get_all_dicts( hist_dict ):
    list_dicts = []
    for key, hist in hist_dict.items():
        if hist == None:
            continue

        if key == 'mig':
          for column in hist:
            for row in column:
              list_dicts.append(asdict(HistData.from_HistList(histlist=row)))
        else:
          list_dicts.append(asdict(HistData.from_HistList(histlist=hist)))
    return list_dicts

def load_events( source ):
  if isinstance(source,list):
    f_list = [ np.load(str(file_one), allow_pickle=True) for file_one in source ]
    events = {}
    for key in list(f_list[0].keys()):
      if key != 'tracks' and key != 'charged':
        events[key] = np.concatenate([ f[key] for f in f_list ], axis=0) 
    from_root = False
  elif '.npz' in source:
    events = np.load(source,allow_pickle=True)
    from_root = False
  else:
    f = uproot.open(source)
    for obj_key in f.keys():
      if 'tree' in obj_key:
        tree = f[obj_key]
        break
    branch_list = []
    for branch in root_var_list:
      if branch in tree.keys():
        branch_list.append(branch)
    events = tree.arrays(branch_list)
    from_root = True

  return events,from_root

def filter_weight(events,filter,weight = None):
  if not isinstance(weight, list):
    weight = [weight]
  weight_array = np.ones(len(filter))
  for w in weight:
    if w is not None:
      if isinstance(w,str):
        weight_array *= np.array(events[w])
      else:
        weight_array *= np.array(w)
  try:
    weight_array = np.array(weight_array[filter])
  except ValueError:
    weight_array = np.array(weight_array[ak.any(filter,axis=1)])
    pass
  return weight_array


def load_var_from_ak( events, cut, var_list, weight = None):
  '''
  events: dictionary of awkward arrays
  cut: CutType object
  var_list: list of VarConfig instances
  weight: the name of the weight branch or an array or a list of weight names or arrays
  '''
  filter = filter_ak_cut(events,cut)
  var_arrays = []
  for var in var_list:
    try:
      var_arrays.append(np.array(ak.flatten(events[var.root_var][filter])))
    except ValueError:
      var_arrays.append(np.array(events[var.root_var][filter]))
  weight_array = filter_weight(events,filter,weight)
  return np.vstack(var_arrays), weight_array

def load_var_from_np( events, cut, var_list, weight = None):
  '''
  events: dictionary of numpy arrays
  cut: CutType object
  var_list: list of VarConfig instances (to prepare for multi-dim ks-test)
  weight: name of weight in the events.keys() or an array or a list of either types
  '''
  filter = filter_np_cut(events,np_cuts[cut])
  var_arrays = []
  for var in var_list:
    var_arrays.append(events[var.np_var][filter])
  weight_array = filter_weight(events,filter,weight)
  return np.vstack(var_arrays), weight_array


def files_in_path( path, sel):
  ls_output = subprocess.getoutput("ls "+os.path.join(path,sel))
  ls_output_list = ls_output.split('\n')
  return ls_output_list





def load_obs_array( events, from_root, obs_list, weight = None ):
  '''
  events: dictionary of numpy arrays or awkward arrays
  from_root: bool
  obs_list: list of ObsConfig
  weight: name of weight or an array or a list of either types
  '''
  gen_var_list = []
  reco_var_list = []
  if isinstance(events,dict):
    var_available = events.keys()
  else:
    var_available = events.fields
  for obs in obs_list:
    if obs.reco.np_var in var_available or obs.reco.root_var in var_available:
      reco_var_list.append(obs.reco)
    if obs.gen.np_var in var_available or obs.gen.root_var in var_available:
      gen_var_list.append(obs.gen)
  obs_array = ObsArray()
  if from_root:
    if len(reco_var_list)>0:
      obs_array.reco,obs_array.weight_reco = load_var_from_ak(events,CutType.PassReco,reco_var_list,weight)
    if len(gen_var_list)>0:
      obs_array.gen,obs_array.weight_gen = load_var_from_ak(events,CutType.PassGen,gen_var_list,weight)
  else:
    if len(reco_var_list)>0:
      obs_array.reco,obs_array.weight_reco = load_var_from_np(events,CutType.PassReco,reco_var_list,weight)
    if len(gen_var_list)>0:
      obs_array.gen,obs_array.weight_gen = load_var_from_np(events,CutType.PassGen,gen_var_list,weight)
  return obs_array

def histlist_dic_applyMCeff(histlist_dic,eff_from_nominal,gen_inveff):
  if eff_from_nominal:
    histlist_dic["gen_inclusive"].get_hist_from_multiplication(histlist_dic["gen_passreco"],gen_inveff)

def histlist_dic_normalize_to_data(histlist_dic,normalization_histlist_dic):
  norm_factor = normalization_histlist_dic.norm / histlist_dic["reco_inclusive"].norm
  for key, hist in histlist_dic.items():
    if not (key == 'mig' or key == "eff" or key == "acc"):
      hist.multiply(norm_factor)



if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--config', default="config/dataset_v7_MCEPOS_unfoldCP1_optimize_1p3M.json", help="The configration file including the unfolding setup")
    parser.add_argument('--method', default="multifold", help="omnifold/multfold/unifold")
    parser.add_argument('--migiter', type=int,nargs='+', help="Which iterations to plot the migration matrices")
    parser.add_argument('--step1', action="store_true", default=False, help="Process the histograms from the step1, otherwise process the step 2 results")
    parser.add_argument('--eff-acc', action="store_true", default=False, help="Consider the efficiency and acceptance in the unfolding")
    parser.add_argument('--eff-from-nominal', action="store_true", default=False, help="Apply the reconstruction efficiency of the nominal MC to the unfolded one, otherwise use the efficiency given by the unfolding algorithm.")
    parser.add_argument('--obs', type=obs_pair, help="Which pair of observables to run.")
    parser.add_argument('--df-tag',type = str, default="debug",help="The tag added to the csv file")
    parser.add_argument('--df-overwrite',action = "store_true", default = False, help = "Overwrite the file storing the dataframe")
    args = parser.parse_args()


    with open(args.config, 'r') as configjson:
        config = json.load(configjson)

    (obs1_name, obs1_bins), (obs2_name, obs2_bins) = parse_obs_args( config, args.obs )

    obs1 = ObsConfig.from_yaml( config["varunfold"], [obs1_name], binning=obs1_bins )
    obs2 = ObsConfig.from_yaml( config["varunfold"], [obs2_name], binning=obs2_bins )

    if not os.path.exists(config["outputdir"]):
      os.makedirs(config["outputdir"])
    mc_infile = get_input_file_from_config_info( config["inputfilesim"], obs2_name )
    fin = ROOT.TFile(mc_infile,"READ")
    tree = fin.Get("ntuplizer/tree") if fin.Get("ntuplizer/tree") else fin.Get("tree")

    weightname=config["reweight"] if ("reweight" in config.keys() and config["reweight"]!="") else "genWeight"

    tree_sys_list=[]
    if config["addsys"]:
        tree_syst_list = get_sys_variations( config, obs2_name )

    all_trees = [tree] + tree_sys_list
    bin_edges_reco, bin_edges_gen = get_bin_edges( config, obs1, obs2, all_trees)
    df_config = { "obs1": args.obs[0],
                  "obs2": args.obs[1],
                  "workflow": config["workflow"],
                  "MC_name": config["MC_name"],
                  "data_name": config["data_name"],
                  "pseudodata_name": config["pseudodata_name"] if config["pseudodata"] else None,
                  "method": args.method,
                  "from_step1": args.step1,
                  "do_eff_acc": args.eff_acc
    }
    mc_hists = fill_hist_lists("MC", obs1, obs2, bin_edges_gen, bin_edges_reco, tree, genWeight=weightname, store_mig=True, **df_config)
    mc_events, mc_from_root = load_events(mc_infile)
    mc_obsarray = load_obs_array( mc_events, mc_from_root, [obs2], weightname )
    #efficiency=reconstructed and generated / generated
    gen_inveff = HistList("HistGenInvEff")
    gen_inveff.read_settings_from_config_dim1(obs1,isgen=True)
    gen_inveff.read_settings_from_config_dim2(obs2,isgen=True)
    gen_inveff.bin_edges_dim2 = bin_edges_gen
    gen_inveff.fill_root_hists_name()
    gen_inveff.get_hist_from_division(mc_hists["gen_inclusive"],mc_hists["gen_passreco"])

    infile_pseudodata = get_input_file_from_config_info( config["inputfilepseudodata"], obs2_name )
    pseudodata_NPZ =  config["pseudodata"] and (isinstance(infile_pseudodata,list) or '.npz' in infile_pseudodata)
    #tree_data, tree_refdata = get_tree_data( config, pseudodata_NPZ)

    infile_refdata = get_input_file_from_config_info( config["inputfiledata"], obs2_name )
    data_NPZ = isinstance(infile_refdata,list) or '.npz' in infile_refdata
    if data_NPZ:
      events_refdata = infile_refdata
    else: 
      fin_refdata=ROOT.TFile(infile_refdata,"READ")
      tree_refdata=fin_refdata.Get("ntuplizer/tree") if fin_refdata.Get("ntuplizer/tree") else fin_refdata.Get("tree")

    if config["pseudodata"]:
      weight_pseudodata = None
      if pseudodata_NPZ:
        if config.get("pseudodataweight",None) is not None:
          file_weight_pseudodata=np.load(config["pseudodataweight"],allow_pickle=True) 
          weight_pseudodata=file_weight_pseudodata[-1]
        events_pseudodata = infile_pseudodata
        from_root = False
      else:
        fin_pseudodata=ROOT.TFile(infile_pseudodata,"READ")
        events_pseudodata = fin_pseudodata.Get("ntuplizer/tree") if  fin_pseudodata.Get("ntuplizer/tree") else  fin_pseudodata.Get("tree")
        from_root = True
      pseudo_hists = fill_hist_lists("Pseudodata", obs1, obs2, bin_edges_gen, bin_edges_reco, events_pseudodata, genWeight=weightname, from_root=from_root, weight_array=weight_pseudodata, store_mig=True,**df_config)
      tag = "_Ref"
      chi2_mc_pseudo = Chi2Collections.from_source(mc_hists,pseudo_hists,"MC","pseudodata")
      pseudo_events,pseudo_from_root = load_events(infile_pseudodata)
      pseudo_obsarray = load_obs_array( pseudo_events, pseudo_from_root, [obs2], [weightname,weight_pseudodata] )
      ks_mc_pseudo = KSDistanceCollections.from_source(mc_obsarray, pseudo_obsarray, "MC","pseudodata")
      target_hists_tuple = ("pseudodata",pseudo_hists, pseudo_obsarray)
      chi2_mc_target = chi2_mc_pseudo
      ks_mc_target = ks_mc_pseudo
    else:
      tag = ""

    weight_data = None
    if config.get("dataweight",None) is not None:
      file_weight_data=np.load(config["dataweight"],allow_pickle=True)
      weight_data=file_weight_data[-1]
    if data_NPZ:
      reco_data_events = events_refdata
      from_root = False
    else:
      reco_data_events = tree_refdata
      from_root = True
      
    data_hists = fill_hist_lists("Data", obs1, obs2, None, bin_edges_reco, reco_data_events,from_root=from_root,weight_array=weight_data, tag=tag, reco_only=True,**df_config)
    chi2_mc_data = Chi2Collections.from_source(mc_hists,data_hists,"MC","data")
    data_events,data_from_root = load_events(infile_refdata)
    data_obsarray = load_obs_array( data_events, data_from_root, [obs2])
    ks_mc_data = KSDistanceCollections.from_source(mc_obsarray, data_obsarray, "MC","data")


    if not config["pseudodata"]:
      chi2_mc_target = chi2_mc_data
      target_hists_tuple = ("data",data_hists,data_obsarray)
      ks_mc_target = ks_mc_data


    normalization_hist = pseudo_hists["reco_inclusive"] if config["pseudodata"] else data_hists["reco_inclusive"]
    histlist_dic_normalize_to_data(mc_hists,normalization_hist)

    if not os.path.exists(config[args.method]["weight"]):
      print("Cannot find weight file ",config[args.method]["weight"])
      exit(0)

    weights=np.load(config[args.method]["weight"],allow_pickle=True)
    if config[args.method]["addsysbs"]:
      #weights_sysbs = []
      ibs_list = []
      sysbs_name_list = files_in_path( config[args.method]['sysbsweightdir'], "*MCbs*npy")
      for sysbs_name in sysbs_name_list:
        ibs_list.append(re.findall(r'\d+', sysbs_name)[-1])
      #  print(f"opening {sysbs_name}")
      #  with open(sysbs_name,'rb') as sysbs_file:
      #    weights_sysbs.append(np.load(sysbs_file,allow_pickle=True))



    weights_per_iter = 4 if args.eff_acc else 2

    step1_tag = "_step1" if args.step1 else ""
    fout = ROOT.TFile(f'{config["outputdir"]}/unfold_{obs1_name}_{obs2_name}_{config["MCtag"]}_optimize_{args.method}{step1_tag}.root', "recreate")
    csv_out = f'{config["outputdir"]}/{config["MCtag"]}_{args.df_tag}.csv'
    list_dict_out = []

    niter = len(weights)//weights_per_iter
    list_chi2_unfold_target = []
    list_ks_distance_unfold_target = []
    for i in range(0,niter+1):
      offset = weights_per_iter // 2 if (args.step1 and i > 0 ) else 0
      weight_index = weights_per_iter*i - offset
      weight_iter = weights[weight_index ]

      store_mig = i in args.migiter

      unfold_hists = fill_hist_lists("MC_"+args.method, obs1, obs2, bin_edges_gen,bin_edges_reco,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_iter,store_mig=store_mig,tag="_iter"+str(i),**df_config)
      unfold_events,unfold_from_root = load_events(config[args.method]["sim"])
      unfold_obsarray = load_obs_array( unfold_events,unfold_from_root, [obs2], [weightname, weight_iter] )
      histlist_dic_applyMCeff(unfold_hists,args.eff_from_nominal,gen_inveff)

      histlist_dic_normalize_to_data(unfold_hists,normalization_hist)

      list_chi2_unfold_target.append(Chi2Collections.from_source(unfold_hists,target_hists_tuple[1],"unfold",target_hists_tuple[0]))
      list_ks_distance_unfold_target.append(KSDistanceCollections.from_source(unfold_obsarray,target_hists_tuple[2],"unfold",target_hists_tuple[0]))

      write_all_hists(unfold_hists)
      list_dict_out += get_all_dicts(unfold_hists)
      if config[args.method]["addsysbs"]:
        #weight_sysbs_iter = [oneweight_sysbs[weight_index] for oneweight_sysbs in weights_sysbs if weight_index<len(oneweight_sysbs)  ]
        #unfold_hists_sysbs = fill_hist_lists("MC_"+args.method+"_sys", obs1, obs2, bin_edges_gen,bin_edges_reco,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_sysbs_iter,store_mig=False,tag="_iter"+str(i),**df_config)
        #histlist_dic_applyMCeff(unfold_hists_sysbs,args.eff_from_nominal,gen_inveff)
        #histlist_dic_normalize_to_data(unfold_hists_sysbs,normalization_hist)
        #write_all_hists(unfold_hists_sysbs)
        #list_dict_out += get_all_dicts(unfold_hists_sysbs)


        def process_bs(ibs):
          with open(sysbs_name_list[ibs],'rb') as sysbs_file:
            weight_sysbs = np.load(sysbs_file,allow_pickle=True)
            if weight_index<len(weight_sysbs):
              weight_sysbs_thisiter = weight_sysbs[weight_index]
            else:
              return None
          unfold_hists_sysbs_i = fill_hist_lists("MC_"+args.method+"_sysbs_"+ibs_list[ibs], obs1, obs2, bin_edges_gen,bin_edges_reco,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_sysbs_thisiter,store_mig=False,tag="_iter"+str(i),**df_config)
          histlist_dic_applyMCeff(unfold_hists_sysbs_i,args.eff_from_nominal,gen_inveff)
          histlist_dic_normalize_to_data(unfold_hists_sysbs_i,normalization_hist)
          write_all_hists(unfold_hists_sysbs_i)
          return get_all_dicts(unfold_hists_sysbs_i)

        with Pool() as pool:
          for result in pool.map(process_bs,range(len(sysbs_name_list))):
            if result is not None:
              list_dict_out += result

        #unfold_hists_sysbs_list = [fill_hist_lists("MC_"+args.method+"_sysbs_"+ibs_list[ibs], obs1, obs2, bin_edges_gen,bin_edges_reco,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_sysbs_iter[ibs],store_mig=False,tag="_iter"+str(i),**df_config) for ibs in range(len(weight_sysbs_iter))]

        #for ibs in range(len(weight_sysbs_iter)):
        #  histlist_dic_applyMCeff(unfold_hists_sysbs_list[ibs],args.eff_from_nominal,gen_inveff)
        #  histlist_dic_normalize_to_data(unfold_hists_sysbs_list[ibs],normalization_hist)
        #for ibs in range(len(weight_sysbs_iter)):
        #  write_all_hists(unfold_hists_sysbs_list[ibs])
        #  list_dict_out += get_all_dicts(unfold_hists_sysbs_list[ibs])

    chi2_unfold_target = Chi2Collections.merge(list_chi2_unfold_target)
    list_dict_out += fill_gof_histdata(chi2_unfold_target,**df_config)

    ks_distance_unfold_target = KSDistanceCollections.merge(list_ks_distance_unfold_target)
    list_dict_out += fill_gof_histdata(ks_distance_unfold_target,**df_config)

    write_all_hists(mc_hists)
    list_dict_out += get_all_dicts(mc_hists)
    list_dict_out += fill_gof_histdata(chi2_mc_target,**df_config)
    list_dict_out += fill_gof_histdata(ks_mc_target,**df_config)


    gen_inveff.write_hist_list()
    write_all_hists(target_hists_tuple[1])
    list_dict_out += get_all_dicts(target_hists_tuple[1])

    if args.df_overwrite:
      print("Overwriting the file",csv_out)
      if os.path.exists(csv_out):
        os.system(f"rm {csv_out}")
      mode = 'w'
    else:
      print("Appending to the file",csv_out)
      mode = 'a'
    pd.DataFrame(data=list_dict_out).to_csv(csv_out,mode=mode,index=False,header=not os.path.exists(csv_out))
