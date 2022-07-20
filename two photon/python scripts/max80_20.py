#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import DivergingNorm
import matplotlib.colors as colors
import numpy as np
from sklearn.cluster import KMeans
from tslearn.clustering import TimeSeriesKMeans, KShape, KernelKMeans
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
import peakutils
from scipy.ndimage.filters import uniform_filter1d
import pickle

with open('NaCl_VGDTR_sorted.p', 'rb') as f:
    test = pickle.load(f)
'''
The dataframes in the pickle file correspond to

0:S_df, 1:dff_df, 2:dff_i, 3:dff_nc, 4:dff_a, 5:i, 6:nc, 7:a

You care about 0, which is your sorted dataframe, 5 which is inhibitory cluster, 6 is no change clust, and 7 which is activated cluster
'''    
df=test[7]





# In[5]:


##elbow method for clustering 
wcss = []

for i in range(1, 11):
    model = KMeans(n_clusters = i,     
                    init = 'k-means++',                 # Initialization method for kmeans
                    max_iter = 300,                     # Maximum number of iterations 
                    n_init = 10,                        # Choose how often algorithm will run with different centroid 
                    random_state = 0)                   # Choose random state for reproducibility
    model.fit(df)                              
    wcss.append(model.inertia_)
    
# Show Elbow plot
plt.plot(range(1, 11), wcss)
plt.title('Elbow Method')                               # Set plot title
plt.xlabel('Number of clusters')                        # Set x axis name
plt.ylabel('Within Cluster Sum of Squares (WCSS)')      # Set y axis name
plt.show()


# In[2]:



preprocessing_meanvar = False # True to use TimeSeriesScalerMeanVariance preprocessing
smooth_n = 5 # n observations to smooth over
smooth_func = 'mean' # one of ['mean','min','max','sum']
norm = True # normalize the data to 0-1 range
model = ''
n_clusters = 2 # number of clusters to fit


if smooth_n > 0:
    if smooth_func == 'mean':
        df = df.rolling(smooth_n).mean().dropna(how='all')
    elif smooth_func == 'max':
        df = df.rolling(smooth_n).max().dropna(how='all')
    elif smooth_func == 'min':
        df = df.rolling(smooth_n).min().dropna(how='all')
    elif smooth_func == 'sum':
        df = df.rolling(smooth_n).sum().dropna(how='all')
    else:
        df = df.rolling(smooth_n).mean().dropna(how='all')

# normalize the data if specified
if norm:
    df = (df-df.min())/(df.max()-df.min())
    

# get values to cluster on
X = df.transpose().values
if preprocessing_meanvar:
    X = TimeSeriesScalerMeanVariance().fit_transform(X)
    df = pd.DataFrame(X.reshape(df.shape), columns=df.columns, index=df.index)
if model == 'kshape':
    model = KShape(n_clusters=n_clusters, max_iter=10, n_init=2).fit(X)
elif model == 'kmeans':
    model = TimeSeriesKMeans(n_clusters=n_clusters, metric="euclidean", max_iter=10, n_init=2).fit(X)
elif model == 'dtw':
    model = TimeSeriesKMeans(n_clusters=n_clusters, metric="dtw", max_iter=5, n_init=2).fit(X)
elif model == 'kernelkmeans':
    model = KernelKMeans(n_clusters=n_clusters, kernel="gak", max_iter=5, n_init=2).fit(X)
else:
    model = TimeSeriesKMeans(n_clusters=n_clusters, metric="euclidean", max_iter=10, n_init=2).fit(X)
    

# build helper df to map metrics to their cluster labels
df_cluster = pd.DataFrame(list(zip(df.columns, model.labels_)), columns=['metric', 'cluster'])

# make some helper dictionaries and lists
cluster_metrics_dict = df_cluster.groupby(['cluster'])['metric'].apply(lambda x: [x for x in x]).to_dict()
cluster_len_dict = df_cluster['cluster'].value_counts().to_dict()
clusters_dropped = [cluster for cluster in cluster_len_dict if cluster_len_dict[cluster]==1]
clusters_final = [cluster for cluster in cluster_len_dict if cluster_len_dict[cluster]>1]
clusters_final.sort()

df_cluster=df_cluster.set_index('metric')


# In[3]:


act_0=df.T.loc[df_cluster['cluster'] == 0].T
act_1=df.T.loc[df_cluster['cluster'] == 1].T


# In[21]:


class peaks():
    def __init__(self, df, stim_0, stim_3):
        self.df = df
        self.normie()
        self.stim_0 = stim_0
        self.stim_3 = stim_3
        self.peak_calc()
        self.risetime()
    
    def normie(self):
        self.df_norm = (self.df - self.df.min() ) / (self.df.max() - self.df.min() )
        self.df_norm_roll = self.df_norm.rolling(30, min_periods = 1).mean()
        self.df_norm_roll_mean = self.df_norm_roll.mean(axis = 1)
        return self.df_norm, self.df_norm_roll, self.df_norm_roll_mean
    
    def plot_mean(self):
        ax = self.df_norm_roll_mean.plot( figsize = (20,10), fontsize = 16)
        ax.set_xlabel('Time (min)', fontsize = 16)
        ax.set_ylabel('Zscore', fontsize = 16)
        
    def peak_calc(self):
        
        ##First off, you need to find the max val, the index of the max val, and the estimated 80% and 20% of that max val
        self.list_attributes = {}
        self.list_calculated80_20 = {}
        
        for i in self.df_norm_roll.columns:
            self.max_index = self.df_norm_roll[i].loc[ self.stim_0 : self.stim_3 ].idxmax()
            self.max_val = self.df_norm_roll[i].loc[ self.stim_0 : self.stim_3 ].max()
            self.max_80 = 0.8*self.max_val
            self.max_20 = 0.2*self.max_val
            name = int(i)
            self.list_attributes[name] = self.max_index, self.max_val
            self.list_calculated80_20[name] = self.max_80, self.max_20
            
        self.max_index_df=pd.DataFrame(self.list_attributes)[:1] # here is your dataframe with the index values you want to stay below to find your 80%
        self.max_80=pd.DataFrame(self.list_calculated80_20)[:1] # dataframe to isolate max80 estimated values
        
        ## make a dictionary so you can easily access each column and value separately 
        bob={}
        for i in self.max_80.columns:
            new = self.max_80[i].values
            bob[i] = new
        
        # get the value for each column that is closest to the estimated max 80 and get the index for that value
        self.list_real80idx={}

        for i in self.max_index_df.columns:
            self.idx80 = abs(self.df_norm_roll[i].loc[ self.stim_0 : (float(self.max_index_df[i])) ] - (bob[i])).idxmin()
            self.real_80 = self.df_norm_roll[i].loc[self.idx80]
            name=int(i)
            self.list_real80idx[name] = self.idx80, self.real_80
        
        
        self.idx80_df = pd.DataFrame(self.list_real80idx)[:1] #make dataframe for indexes of the max 80 values so that you can stay to the left for max20%
        self.max_20 = pd.DataFrame(self.list_calculated80_20)[1:] #dataframe to isolate max20 estimated values
        
        ## make a dictionary so you can easily access each column and value separately 
        boob={}
        for i in self.max_20.columns:
            new = self.max_20[i].values
            boob[i] = new
        
        # get the value for each column that is closest to the estimated max 20 and get the index for that value
        self.list_real20idx = {}
        for i in self.idx80_df.columns:
            self.idx20 = abs(self.df_norm_roll[i].loc[self.stim_0:(float(self.idx80_df[i].values))]-(boob[i])).idxmin()
            self.real_20=self.df_norm_roll[i].loc[self.idx20]
            name=int(i)
            self.list_real20idx[name] = self.idx20, self.real_20

        self.max_peak_wholetrace = {}
        for i in self.df_norm_roll.columns:
            self.max_peak_idx = self.df_norm_roll[i].idxmax()
            self.max_peak_val = self.df_norm_roll[i].max()
            name=int(i)
            self.max_peak_wholetrace[name] = self.max_peak_idx, self.max_peak_val
        
        return self.list_attributes[name], self.list_calculated80_20[name], self.list_real80idx[name], self.list_real20idx[name], self.max_peak_wholetrace[name]
    
    def plot_peaks(self):
        
        plt.rcParams["figure.figsize"] = [5,5]
        plt.rcParams["font.size"] = 18

        for i in self.df_norm_roll.columns: #[0:len(act_0_norm_roll.columns)]:
    
            plt.plot(self.df_norm_roll[i], 'k')
            plt.xlim([0,45])
            plt.plot(self.list_attributes[i][0],self.list_attributes[i][1] , 'X', color='orange' )
            plt.plot(self.list_real80idx[i][0], self.list_real80idx[i][1], "<", color="r")
            plt.plot(self.list_real20idx[i][0], self.list_real20idx[i][1], "<", color="m")
            plt.plot(self.max_peak_wholetrace[i][0], self.max_peak_wholetrace[i][1], "*", color='royalblue')
            plt.xlabel("Time (m)")
            plt.ylabel("Z-Score")
            #plt.plot(list_real20idx[i], ">")
            plt.show()
            
    def risetime(self):
        self.RT=(pd.DataFrame(self.list_real80idx)-pd.DataFrame(self.list_real20idx))[:1]
        return self.RT
    
    def slewtime(self):
        self.slew_time = {}
        for i in self.RT.columns:
            new = self.RT[i]
            st = ((self.list_real80idx[i][1] - self.list_real20idx[i][1]) / (new))[0]# (self.list_real80idx[i][1] - self.list_real20idx[i][1]) 
            self.slew_time[i]=st
        return self.slew_time
    

    def peak_per(self):
        peak_per=((pd.DataFrame(self.list_real80idx))[1:]/pd.DataFrame(self.max_peak_wholetrace)[1:])*100
        return peak_per    
            
            

    
    
nacl=peaks(act_0, 19, 40)


# In[ ]:





# In[ ]:




