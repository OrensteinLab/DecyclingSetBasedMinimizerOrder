import sys,os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


labels_dict = {'random':'Random', 'decycling': 'Decycling', 'double' : 'Double decycling', 'pasha' : 'UHS (PASHA)', 'Miniception' : 'Miniception', 'docks' : 'UHS (DOCKS)', 'lower_bound' : 'Lower bound'}
# 'decycling_set' : 'Pre-computed decycling',

color_dict = {'random' : 'tab:red', 'docks' : 'black', 'pasha' : 'tab:orange', 'decycling' : 'blue', 'double' : 'lime', 'lower_bound' : 'grey', 'Miniception' : 'tab:cyan' }
#markers_dict = {'random' : 'p', 'docks' : 's', 'pasha' : 'o', 'decycling' : 'D', 'double' : '|' }
lines_dict = {'random' : '-', 'docks' : ':', 'pasha' : (0, (5, 10)), 'decycling' : '-.', 'double' : (0, (5, 1)), 'lower_bound' : (0, (3, 1, 1, 1, 1, 1)), 'Miniception' : (0, (3, 5, 1, 5, 1, 5)) }

ks = [11,20,50,100]

Ls = 100


def plot_k(k,df,seq='expected',lowerbound=False,stds=False):
  plt.figure(figsize=(10, 7), dpi=600)
  dfk = df[df['k'] == k]
  if len(dfk)==0: return
  methods = dfk['method'].unique()
  methods = [m for m in methods if m in labels_dict.keys()]
  legend = []
  for method in methods:
    method_df = dfk[dfk['method']==method]
    Ls = np.array(sorted(method_df['L'].unique()))
    dens = np.array([method_df[method_df['L']==L]['density'].values[0] for L in Ls])
    stddevs = np.array([method_df[method_df['L']==L]['density_std'].values[0] for L in Ls])
    factors = (Ls-k+2)
    if not stds:
      plt.plot(Ls,dens*factors,color=color_dict[method],linestyle=lines_dict[method],linewidth=2.5)
    else:
      plt.errorbar(Ls,dens*factors,yerr=stddevs*factors, color=color_dict[method], linestyle=lines_dict[method],linewidth=2.5)
    legend.append(labels_dict[method])

  if seq == 'expected':
    title = "k={}".format(k)#"Expected density factor for k={}".format(k)
    plt.ylabel("Expected density factor",fontsize=22)
  else:
    title = "k={}, {}".format(k,seq)#"Density factor on {} for k={}".format(seq,k)
    plt.ylabel("Particular density factor",fontsize=22)
  plt.title(title,fontsize=24)
  plt.xlabel("L",fontsize=22)
  if lowerbound:
      Ls = np.array(sorted(dfk['L'].unique()))
      method = 'lower_bound'
      ws = Ls-k
      maxes = np.array([max(0,(k-w)/(w)) for w  in ws])
      bounds = (ws+1)*(1.5+maxes+1/(2*(ws)))/Ls
      if not stds:
          plt.plot(Ls,bounds,color = color_dict[method],linestyle=lines_dict[method],linewidth=2.5)
      else:
          plt.errorbar(Ls,bounds,color = color_dict[method],linestyle=lines_dict[method],linewidth=2.5)
      legend.append(labels_dict[method])
   
  xtic = np.arange((min(dfk['L'])//10)*10,max(dfk['L'])+1,10)
  if len(xtic) > 20:
      xtic = np.arange((min(Ls)//10)*10,max(Ls)+1,20)
  plt.xticks(xtic,rotation=45,fontsize=18)
  plt.yticks(fontsize=18)

  plt.legend(legend,loc='lower right',fontsize=16)#, bbox_to_anchor=(1, 0.5))
  plt.tight_layout()
  outname = "/home/gaga/data-scratch/dpellow/uhsnn/decycling/outputs/rerun_figs/{}_{}.pdf".format(title,stds)
  plt.savefig(outname)
  plt.clf()


def plot_L(L,df,seq='expected',lowerbound=False,stds=False,max_k=15):
  plt.figure(figsize=(10, 7), dpi=600)
  dfL = df[df['L'] == L]
  methods = dfL['method'].unique()
  methods = [m for m in methods if m in labels_dict.keys()]
  legend = []
  for method in methods:
    method_df = dfL[dfL['method']==method]
    ks = np.array(sorted(method_df['k'].unique()))
    ks = np.array([k for k in ks if k<=max_k])
    dens = np.array([method_df[method_df['k']==k]['density'].values[0] for k in ks])
    stddevs = np.array([method_df[method_df['k']==k]['density_std'].values[0] for k in ks])
    factors = (L-ks+2)
    if not stds:
      plt.plot(ks,dens*factors,color=color_dict[method],linestyle=lines_dict[method],linewidth=2.5)
    else:
      plt.errorbar(ks,dens*factors,yerr=stddevs*factors, color=color_dict[method], linestyle=lines_dict[method],linewidth=2.5)
    legend.append(labels_dict[method])

  if seq == 'expected':
    title = "L={}".format(L)#"Expected density factor for L={}".format(L)
    plt.ylabel("Expected density factor",fontsize=22)
  else:
    title = "L={}, {}".format(L,seq)#"Density factor on {} for L={}".format(seq,L)
    plt.ylabel("Particular density factor",fontsize=22)
  plt.title(title,fontsize=24)
  plt.xlabel("k",fontsize=22)

  if lowerbound:
      ks = np.array(sorted(dfL['k'].unique()))
      ks = np.array([k for k in ks if k<=max_k])  
      method = 'lower_bound'
      ws = L-ks
      maxes = np.array([max(0,(ks[i]-w)/(w)) for i,w in enumerate(ws)])
      bounds = (ws+1)*(1.5+maxes+1/(2*(ws)))/L
      if not stds:
          plt.plot(ks,bounds,color = color_dict[method],linestyle=lines_dict[method],linewidth=2.5)
      else:
          plt.errorbar(ks,bounds,color = color_dict[method],linestyle=lines_dict[method],linewidth=2.5)
      legend.append(labels_dict[method])
   
  xtic = np.arange(min(dfL['k']),min(max(dfL['k']),max_k)+1,1)
  if len(xtic) > 20:
        xtic = np.arange(min(dfL['k']),min(max(dfL['k']),max_k)+1,2)
  plt.xticks(xtic,rotation=45,fontsize=18)
  plt.yticks(fontsize=18)

  plt.legend(legend, fontsize=16, loc='upper right')#, bbox_to_anchor=(1, .35))# loc='lower left')#'center right', bbox_to_anchor=(1, .35))# ,loc='lower right'
  plt.tight_layout()
  outname = "/home/gaga/data-scratch/dpellow/uhsnn/decycling/outputs/rerun_figs/{}_{}.pdf".format(title,stds)
  plt.savefig(outname)
  plt.clf()    



if __name__ == '__main__':


  infile = sys.argv[1]
  inseq = 'expected'

  if len(sys.argv) > 2:
    inseq = sys.argv[2]
      
  df = pd.read_csv(infile,sep='\t')
  for k in ks:
    plot_k(k,df,seq=inseq,lowerbound=True)
    plot_k(k,df,stds=True,seq=inseq)
  plot_L(100,df,seq=inseq)
  plot_L(100,df,stds=True,seq=inseq)




