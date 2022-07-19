#!/usr/bin/env python
# coding: utf-8

# In[22]:


import pandas as pd
import numpy as np
import pickle
import os 

#F_path=os.path.abspath("D:\Euhydrated heating\VGLUTDTRM1\Heating April 22 2022\suite2p\plane0\F.npy")
#Fneu_path=os.path.abspath("D:\Euhydrated heating\VGLUTDTRM1\Heating April 22 2022\suite2p\plane0\Fneu.npy")
#iscell_path=os.path.abspath("D:\Euhydrated heating\VGLUTDTRM1\Heating April 22 2022\suite2p\plane0\iscell.npy")

F=np.load('NaCl infusion-VGLUTDTRM1/F.npy', allow_pickle=True) ##input local file path, or variable name if file not locally stored in jupyter
Fneu=np.load('NaCl infusion-VGLUTDTRM1/Fneu.npy', allow_pickle=True)
iscell=np.load('NaCl infusion-VGLUTDTRM1/iscell.npy', allow_pickle=True)
spks=np.load('NaCl infusion-VGLUTDTRM1/spks.npy', allow_pickle=True)

iscell=pd.DataFrame(iscell)
F=pd.DataFrame(F)
Fneu=pd.DataFrame(Fneu)

F_cells=F.loc[(iscell[0] == 1)]
Fneu_cells=Fneu.loc[(iscell[0] == 1)]

#corrected fluorescence
F_corr=F_cells-0.7*Fneu_cells
F_corr=F_corr.T
F_corr.index=F_corr.index/60/5

#zscored
p=(F_corr-(F_corr.loc[1:15].mean(axis=0)))/(F_corr.loc[1:15].std())


# In[85]:


class sorting_zscore():
    def __init__(self, df_name, stim_start, stim_end, inhib_thresh, act_thresh):
        self.df_name=df_name
        self.stim_start=stim_start
        self.stim_end=stim_end
        self.inhib_thresh=inhib_thresh
        self.act_thresh=act_thresh
        
    def sort_df_by_dff(self):
        x=self.df_name.loc[self.stim_start:self.stim_end]
        y=x.mean(axis=0)
        dff=pd.DataFrame(y)
        dff.columns=['dff']
        self.dff=dff.sort_values(by=['dff'])
        y=self.df_name.T.reindex(self.dff.index)
        return y, self.dff
    
    def sort_dff_by_cat(self):
        self.dff_inhib=self.df_name.T.loc[self.dff['dff'] <=self.inhib_thresh]
        self.dff_nada=self.df_name.T.loc[(self.dff['dff'] >= self.inhib_thresh) & (self.dff['dff']<=self.act_thresh)]
        self.dff_act=self.df_name.T.loc[self.dff['dff']>=self.act_thresh]
        return self.dff_inhib, self.dff_nada, self.dff_act
    
    def sort_df_by_cat(self):
        inhib=self.df_name.T.loc[self.dff_inhib.index]
        nada=self.df_name.T.loc[self.dff_nada.index]
        act=self.df_name.T.loc[self.dff_act.index]
        return inhib.T, nada.T, act.T


        
    
new_df=sorting_zscore(p, 20,45,(-0.5),1)

[S_df,dff_df]=new_df.sort_df_by_dff()
[dff_i, dff_nc, dff_a]=new_df.sort_dff_by_cat()
[i, nc, a]=new_df.sort_df_by_cat()

lst=[S_df,dff_df,dff_i, dff_nc, dff_a, i, nc, a]

with open('NaCl_VGDTR_sorted.p', 'wb') as f:
    pickle.dump(lst,f)
    

