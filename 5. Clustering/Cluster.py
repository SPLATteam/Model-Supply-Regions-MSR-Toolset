#Li Yunshu, Oct 2021, for IRENA
import pandas as pd
from tslearn.clustering import TimeSeriesKMeans
import matplotlib.pyplot as plt
import numpy as np
import re
import os

# USER INPUTS

ctrl_path = "ctrl.csv"
ctrl = pd.read_csv(ctrl_path).set_index('ctrl')
f_in = ctrl.loc['f_in']['val']
it = int(ctrl.loc['iterations']['val'])     #number of iterations in the clustering
ctrl_country_in = ctrl.loc['ctrl_country_in']['val'] #input path for countries file
tech = ctrl.loc['tech']['val']  #tech the run is about
metric = ctrl.loc['cluster_metric']['val']  #default is euclidean
ctrl_param_in = ctrl.loc['ctrl_param_in']['val']  #parameters' mode of aggregation
ctryf_writeout = ctrl.loc['ctryf_writeout']['val']  #if write to separate country folders

raw = pd.read_excel(f_in)
raw = raw.set_index('Zone_ID')
#print(raw)

#params_to_output = ['Area sq.Km', 'Max Capacity MW', 'DistFromTL Km'] #not used yet

ctrl_ctry = pd.read_csv(ctrl_country_in).set_index('country').to_dict() # load number of clusters for each country from file
nclust_all = ctrl_ctry['n_cluster'] # number of clusters

ctrl_param = pd.read_csv(ctrl_param_in)
sum_param = ctrl_param[ctrl_param['mode'] == 'sum']['param']
mean_param = ctrl_param[ctrl_param['mode'] == 'mean']['param']
wmean_param = ctrl_param[ctrl_param['mode'] == 'wmean']['param']
max_param = ctrl_param[ctrl_param['mode'] == 'max']['param']
wmean_IEC_param = ctrl_param[ctrl_param['mode'] == 'wmean_IEC']['param']

debug = pd.DataFrame()
out_dir = tech
os.mkdir(out_dir)  # make output directory

for ctry in nclust_all:

    #number of clusters
    nclust = nclust_all[ctry]

    #load country data
    raw_ctry = raw[raw['CtryName'] == ctry]
    if not raw_ctry.empty:
        raw_trans = raw_ctry.T
        raw_trans_ind = raw_trans.index.tolist()
        raw_trans['H'] = [int(re.findall('\d+', i)[0]) if re.findall('H\d+', i) else 'nan' for i in raw_trans_ind]
        raw_trunc = raw_trans[raw_trans['H'] != 'nan'].set_index('H')
        raw_trunc_inv = raw_trans[raw_trans['H'] == 'nan']

        #K-Means Clustering
        model = TimeSeriesKMeans(n_clusters=nclust, metric=metric, max_iter=it, verbose=1).fit(raw_trunc.T)
        labels = model.labels_

        map = dict(zip(raw_trunc.columns,labels)) #not used, write out only

        #Make write out folder
        saveas = '{} {} zones {}clust_{}it'.format(ctry, tech, nclust,it)

        #Write out original parameters and map to clusters
        all_table_ = raw_trunc_inv.T.copy()
        all_table_ = all_table_[all_table_.index != 'H']
        all_table_['cluster'] = labels
        try:
            all_table_all = pd.concat([all_table_all, all_table_])
        except:
            all_table_all = all_table_
    #    debug = debug.append(all_table_)

        all_table = all_table_.groupby(['cluster'])
#
#         #Aggregate parameters
        clust_table = pd.DataFrame(index=range(nclust))
        clust_table['Zones count'] = all_table['Area sq.Km'].count()
#
        for p in sum_param:
            clust_table[p] = all_table[p].sum()
        for p in mean_param:
            clust_table[p] = all_table[p].mean()
        for p in max_param:
            clust_table[p] = all_table[p].max()
#
        all_table_['weights'] = len(all_table_)
        for c in range(nclust):
            all_table_.loc[all_table_['cluster']==c,'weights'] = all_table_[all_table_['cluster']==c]['Max Capacity MW']/clust_table.loc[c]['Max Capacity MW']
#
        for p in wmean_param:
            p_ = p + '_'
            all_table_[p + '_'] = all_table_[p] * all_table_['weights']
            clust_table[p] = (all_table_.groupby(['cluster']))[p+'_'].sum()

        try:

            for p in wmean_IEC_param:
                p_ = p + '_'
                all_table_[p + '_'] = all_table_[p].str.replace(r'\D', '').astype(int) * all_table_['weights']
                clust_table[p] = 'Class-' + round((all_table_.groupby(['cluster']))[p+'_'].sum()).astype(int).astype(str)
        except:
            print('WARNING: no parameter using wmean_IEC, if running for wind, check if this is defined. Perhaps you are running for solar?')

#
#        if ctryf_writeout: clust_table.to_csv(out_dir +'/'+ saveas+'_params.csv')
#
#         #Aggregate capacity factors
        all_table_profiles = raw_trunc.copy().T.astype(float)
        all_table_profiles = all_table_profiles.multiply(all_table_['weights'], axis=0)
        all_table_profiles['cluster'] = labels
        all_table_profiles = all_table_profiles.groupby(['cluster']).sum()

#        if ctryf_writeout: all_table_profiles.to_csv(out_dir +'/'+saveas + '_profiles.csv')

        #Plot results and save
        plot_count = nclust
        frac = [0, 1/4, 2/4, 3/4]
        plot_hours = [int(8760* i) for i in frac]
        #
        fig, axs = plt.subplots(4, plot_count, figsize=(40, 25))
        fig.suptitle('Clusters')
        # For each label there is,
        # plots every series with that label
        for label in set(labels):
            column_j = label
            cluster = []
            for i in range(len(labels)):
                if (labels[i] == label):
                    for row_i in range(4):
                        axs[row_i, column_j].plot(raw_trunc.iloc[plot_hours[row_i]:plot_hours[row_i]+72,i], c="gray", alpha=0.4)
                        cluster.append(raw_trunc.iloc[:,i])
            if len(cluster) > 0:
                for row_i in range(4):
                    axs[row_i, column_j].plot(range(1+plot_hours[row_i],1+plot_hours[row_i]+72),all_table_profiles.loc[label][plot_hours[row_i]:plot_hours[row_i]+72], c="red")
        #
        fig.savefig(out_dir +'/'+ saveas)
        # #plt.show()
        #
        # #concat
        #
        smallT = pd.concat([clust_table, all_table_profiles], axis=1)
        smallT['CtryName'] = ctry
        try:
            bigT = pd.concat([bigT, smallT], axis=0)
        except:
            bigT = smallT



bigT = bigT.rename_axis('Zone_ID')
bigT.index = bigT.index + 1
cols = bigT.columns.tolist()
cols = cols[-1:] + cols[:-1] #move countryname to first position
bigT = bigT[cols]
bigT.to_csv('bigT.csv')
all_table_all.to_csv('unclustered_zones_all.csv')
