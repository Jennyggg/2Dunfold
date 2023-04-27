import numpy as np
from argparse import ArgumentParser
import os
import matplotlib as mpl
mpl.use('Agg')
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import MaxNLocator
hep.style.use("CMS")
import ast
from unfold_utils import HistList

def lumi_text(lumi):
    if np.log10(lumi) >= 0 and np.log10(lumi) < 3:
        lumi_print = lumi
        lumi_unit = '$\mathrm{fb}^{-1}$'
    elif np.log10(lumi) < 6 and np.log10(lumi) >= 3:
        lumi_print = lumi * 0.001
        lumi_unit = '$\mathrm{ab}^{-1}$'
    elif np.log10(lumi) < 0 and np.log10(lumi) >= -3:
        lumi_print = lumi * 1000
        lumi_unit = '$\mathrm{pb}^{-1}$'
    elif np.log10(lumi) < -3 and np.log10(lumi) >= -6:
        lumi_print = lumi * 1000000.0
        lumi_unit = '$\mathrm{nb}^{-1}$'
    elif np.log10(lumi) <= -6 and np.log10(lumi) >= -9:
        lumi_print = lumi * 1000000000.0
        lumi_unit = '$\mathrm{\mu b}^{-1}$'
    elif np.log10(lumi) <= -9 and np.log10(lumi) >= -12:
        lumi_print = lumi * 1000000000000.0
        lumi_unit = '$\mathrm{mb}^{-1}$'
    elif np.log10(lumi) <= -12 and np.log10(lumi) >= -15:
        lumi_print = lumi * 1000000000000000.0
        lumi_unit = '$\mathrm{b}^{-1}$'
    return ('{0:.1f}').format(lumi_print) + ' ' + lumi_unit + ' 2018 (13 TeV)'

def latex_root_to_mpl(text):
  text=text.replace('#','\\')
  text=text.replace(' ','\ ')
  text = "$\mathrm{" + text + "}$"
  return text


#temporary class from root histograms to arrays, will be replaced by boost hisograms
class HistArray:

  def __init__(self, histlist = None, bins = None):
    self._nested_value = [np.array([hist.GetBinContent(ibin+1) for ibin in range(hist.GetNbinsX())]) for hist in histlist.root_hists] if histlist is not None else None
    self._nested_error = [np.array([hist.GetBinError(ibin+1) for ibin in range(hist.GetNbinsX())]) for hist in histlist.root_hists] if histlist is not None else None
    if  histlist is not None:
      self._nested_bins = []
      for hist in histlist.root_hists:
        bins = []
        for ibin in range(hist.GetNbinsX()):
          bins.append(hist.GetXaxis().GetBinLowEdge(ibin+1))
        bins.append(hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX()))
        self._nested_bins.append(np.array(bins))
    else:
      self._nested_bins = None
    self._bins = bins
  @property
  def nested_value(self):
    return self._nested_value

  @property
  def nested_error(self):
    return self._nested_error

  @property
  def nested_bins(self):
    return self._nested_bins

  def bins(self):
    return self._bins

  @classmethod
  def from_df_entry(cls,entry):
    bin_values = entry['bin_values']
    bin_errors = entry['bin_errors']
    bins = entry['bin_edges_dim1']
    nested_bins = entry['bin_edges_dim2']
    if isinstance(bin_values,str):
      bin_values = ast.literal_eval(bin_values.replace('inf','2e308'))
    if not isinstance(bin_values[0],list):
      bin_values = [bin_values]
    else:
      bin_values = [value[1:-1] for value in bin_values]
    if isinstance(bin_errors,str):
      bin_errors = ast.literal_eval(bin_errors.replace('inf','2e308'))
    if not isinstance(bin_errors[0],list):
      bin_errors = [bin_errors]
    else:
      bin_errors = [error[1:-1] for error in bin_errors]
    if isinstance(bins,str):
      bins = ast.literal_eval(bins.replace('inf','2e308'))
    if isinstance(nested_bins,str):
      nested_bins = ast.literal_eval(nested_bins)
    if not isinstance(nested_bins[0],list):
      nested_bins = [nested_bins]

    histarray = cls()
    histarray._nested_value = [np.array(value) for value in bin_values]
    histarray._nested_error = [np.array(error) for error in bin_errors]
    histarray._bins = np.array(bins)
    histarray._nested_bins = [np.array(edges) for edges in nested_bins]
    return histarray


  def __len__(self):
    return len(self._nested_value)

  def __truediv__(self,other):
    result = HistArray()
    result._nested_bins = self._nested_bins
    if isinstance(other,HistArray):
      result._nested_value = [self._nested_value[i] / other._nested_value[i] for i in range(len(self))]
      result._nested_error = [self._nested_error[i] / other._nested_value[i] for i in range(len(self))]
    else:
      result._nested_value = [self._nested_value[i] / other for i in range(len(self))]
      result._nested_error = [self._nested_error[i] / other for i in range(len(self))]
    return result

  def divide_by_bin_width(self):
    for ibin1 in range(len(self._nested_bins)):
      if self._bins is not None:
        self._nested_value[ibin1] = self._nested_value[ibin1] / (self._bins[ibin1+1]-self._bins[ibin1])
        self._nested_error[ibin1] = self._nested_error[ibin1] / (self._bins[ibin1+1]-self._bins[ibin1])
      for ibin2 in range(len(self._nested_bins[ibin1])-1):
        self._nested_value[ibin1][ibin2] = self._nested_value[ibin1][ibin2] / (self._nested_bins[ibin1][ibin2+1]-self._nested_bins[ibin1][ibin2])
        self._nested_error[ibin1][ibin2] = self._nested_error[ibin1][ibin2] / (self._nested_bins[ibin1][ibin2+1]-self._nested_bins[ibin1][ibin2])

def map_marker_style(style):
  marker_dict = { 'marker' :'o', 
                  'cross': '.',
                  'triangle': '^',
                  'triangle_down': 'v' }
  return marker_dict.get( style, None)

def draw_array(value,error,bins,style,ax,color,legend):

  marker = map_marker_style(style)
  if marker is not None:
    center = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
    xerrs = [(bins[i+1] - bins[i])/2 for i in range(len(bins)-1)]
    ax.errorbar(
      center,
      list(value),
      xerr = xerrs,
      yerr = list(error),
      color = color,
      marker = marker,
      label = legend,
      ls = 'none'
    )
  elif style == 'fillederror':
    error_up = np.append(value,value[-1]) + np.append(error,error[-1])
    error_down = np.append(value,value[-1]) - np.append(error,error[-1])
    ax.fill_between(
      x = bins,
      y1 = error_up,
      y2 = error_down,
      color = color,
      step = "post",
      linewidth = 0,
      alpha = 0.8,
      label = legend
    )
    ax.fill_between(
      x = bins,
      y1 = np.append(value,value[-1]),
      y2 = np.append(value,value[-1]),
      color = color,
      step = "post"
    )
  elif style == 'filled':
    error_up = np.append(value,value[-1]) + np.append(error,error[-1])
    error_down = np.append(value,value[-1]) - np.append(error,error[-1])
    ax.fill_between(
      x = bins,
      y1 = error_up,
      y2 = error_down,
      color = color,
      step = "post",
      linewidth = 0,
      alpha = 0.8,
      label = legend
    )
    hep.histplot(
      value-error,
      bins=bins,
      histtype = 'fill',
      color = color,
      alpha=0.5,
      edgecolor=None,
      linewidth = 0,
      ax = ax
    )
  elif style == 'hist':
    hep.histplot(
      value,
      bins=bins,
      histtype = 'step',
      color = color,
      label = legend,
      ax = ax
    )
  elif style == 'line':
    hep.histplot(
      value,
      bins=bins,
      yerr = error,
      histtype = 'step',
      color = color,
      label = legend,
      ax = ax
    )
  elif style == 'hatch':
    error_up = np.append(value,value[-1]) + np.append(error,error[-1])
    error_down = np.append(value,value[-1]) - np.append(error,error[-1])
    ax.fill_between(
      x = bins,
      y1 = error_up,
      y2 = error_down,
      facecolor = "none",
      edgecolor = color,
      step = "post",
      linewidth = 0,
      hatch = "///",
      label = legend
    )
  else:
    print(f"style {style} not recognized")

def get_histarray(source):
  if isinstance(source,HistList):
    return HistArray(source)
  elif isinstance(source,HistArray):
    return source
  else:
    return None

def plot_flat_hists_mpl(hist_ref, list_hist_compare, legend_ref, list_legend_compare, title="", is_logY=0, do_ratio=1, output_path="./", hist_ref_stat=0, text_list=[], style_ref='marker', color_ref="black", list_style_compare=[], list_color_compare=[], labelY='Normalized Events/Bin Width', label_ratio='Data/MC',range_ratio=0.3):

  hist_ref_arrays = get_histarray(hist_ref)
  hist_ref_stat_arrays = get_histarray(hist_ref_stat)
  list_hist_compare_arrays = [get_histarray(hist_compare) for hist_compare in list_hist_compare]
  if do_ratio:
    f, axs = plt.subplots(2,len(hist_ref_arrays),sharex=True,sharey='row',gridspec_kw={"height_ratios": (2,1) },figsize=(12.0+(len(hist_ref_arrays)-1)*2,10.0))
  else:
    f, axs = plt.subplots(1,len(hist_ref_arrays),sharex=True,sharey='row')
  if do_ratio:
    axs = axs.flatten()
  elif len(hist_ref_arrays) == 1:
    axs = [axs]
  hep.cms.label('Preliminary', data=False, rlabel="", loc=0, ax = axs[0],fontsize = 20)
  hep.cms.lumitext(lumi_text(1.4817568788812e-08), ax = axs[len(hist_ref_arrays)-1])

  order_hist_array = list_hist_compare_arrays.copy()
  order_style = list_style_compare.copy()
  order_color = list_color_compare.copy()
  order_legend = list_legend_compare.copy()
  if style_ref == "filled":
     order_hist_array.insert(0,hist_ref_arrays)
     order_style.insert(0,style_ref)
     order_color.insert(0,color_ref)
     order_legend.insert(0,legend_ref)
  else:
     order_hist_array.append(hist_ref_arrays)
     order_style.append(style_ref)
     order_color.append(color_ref)
     order_legend.append(legend_ref)

  for ihist in range(len(hist_ref_arrays)):
    for (hist_array,style,color,legend) in zip(order_hist_array,order_style,order_color,order_legend):
      draw_array(hist_array.nested_value[ihist], hist_array.nested_error[ihist],hist_array.nested_bins[ihist],style, axs[ihist],color,legend) 
      axs[ihist].text(.5,.65,latex_root_to_mpl(text_list[ihist]),horizontalalignment='center',transform=axs[ihist].transAxes,fontsize = 'x-small')
      if ihist != len(hist_ref_arrays)-1 and not do_ratio:
        x_ticks = axs[ihist].xaxis.get_major_ticks()
        x_ticks[-2].label1.set_visible(False)
  axs[0].set_ylabel(labelY)
  axs[0].ticklabel_format(axis='y',style='sci',useMathText=True)
  axs[ihist].legend(facecolor='white', framealpha=0.95)
  bottom, top = axs[ihist].get_ylim()
  axs[ihist].set_yscale('log' if is_logY else 'linear')
  if is_logY:
    axs[ihist].set_ylim(top = top*50)
  else:
    axs[ihist].set_ylim(top = top*1.5, bottom = max(bottom,0))
  if do_ratio:
    order_hist_array_ratio = [hist_compare_arrays / hist_ref_arrays for hist_compare_arrays in list_hist_compare_arrays]
    order_style_ratio = ["fillederror" if style_compare == "filled" else style_compare for style_compare in list_style_compare ]
    order_color_ratio = list_color_compare.copy()
    order_hist_array_ratio.insert(0,hist_ref_arrays / hist_ref_arrays)
    order_style_ratio.insert(0,"fillederror")
    order_color_ratio.insert(0,"darkgrey")
    if hist_ref_stat_arrays is not None:
      order_hist_array_ratio.insert(1,hist_ref_stat_arrays / hist_ref_arrays)
      order_style_ratio.insert(1,"fillederror")
      order_color_ratio.insert(1,"dimgray")
    for iratio in range(len(hist_ref_arrays)):
      for (hist_array,style,color) in zip(order_hist_array_ratio,order_style_ratio,order_color_ratio):
        draw_array(hist_array.nested_value[iratio], hist_array.nested_error[iratio],hist_array.nested_bins[iratio],style, axs[iratio+len(hist_ref_arrays)],color,None)
        if iratio != len(hist_ref_arrays)-1:
          x_ticks = axs[iratio+len(hist_ref_arrays)].xaxis.get_major_ticks()
          x_ticks[-2].label1.set_visible(False)
    axs[len(hist_ref_arrays)].set_ylabel(label_ratio)
  axs[len(axs)-1].set_xlabel(latex_root_to_mpl(title))
  axs[len(axs)-1].set_xlim(left = hist_ref_arrays.nested_bins[0][0], right = hist_ref_arrays.nested_bins[0][-1])
  ratio_bottom, ratio_top = axs[len(axs)-1].get_ylim()
  ratio_range_sym = max(1.0-ratio_bottom,ratio_top-1.0)
  ratio_bottom_sym, ratio_top_sym = 1.0-ratio_range_sym, 1.0+ratio_range_sym
  axs[len(axs)-1].set_ylim(bottom = max(0,ratio_bottom_sym), top = ratio_top_sym)
  plt.subplots_adjust(wspace=0, hspace=0)
  if is_logY:
    output_path += "_logy"
  plt.savefig(output_path + '.pdf')
  plt.savefig(output_path + '.png')
  plt.close()

