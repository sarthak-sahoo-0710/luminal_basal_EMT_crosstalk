#%%

# import packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.special import entr


plt.rcParams.update({'font.size': 18, 'font.family': 'arial'})
relevant_columns = ['HALLMARK_ESTROGEN_RESPONSE_EARLY', 'basal_breast_cancer', 'luminal_breast_cancer', "pEMT", "Epithelial_tumour_signature", "Mesenchymal_tumour_signature"]


#%%
# Load data
data1 = pd.read_csv('./../CID4290_metadata.txt', sep='\t', header=0)
data1 = data1[relevant_columns]
# normalise each column by min-max scaling
data1 = (data1 - data1.min()) / (data1.max() - data1.min())

data2 = pd.read_csv('./../CID4535_metadata.txt', sep='\t', header=0)
data2 = data2[relevant_columns]
# normalise each column by min-max scaling
data2 = (data2 - data2.min()) / (data2.max() - data2.min())

data3 = pd.read_csv('./../CID4465_metadata.txt', sep='\t', header=0)
data3 = data3[relevant_columns]
# normalise each column by min-max scaling
data3 = (data3 - data3.min()) / (data3.max() - data3.min())

data4 = pd.read_csv('./../CID44971_metadata.txt', sep='\t', header=0)
data4 = data4[relevant_columns]
# normalise each column by min-max scaling
data4 = (data4 - data4.min()) / (data4.max() - data4.min())

data5 = pd.read_csv('./../1160920F_metadata.txt', sep='\t', header=0)
data5 = data5[relevant_columns]
# normalise each column by min-max scaling
data5 = (data5 - data5.min()) / (data5.max() - data5.min())

data6 = pd.read_csv('./../1142243F_metadata.txt', sep='\t', header=0)
data6 = data6[relevant_columns]
# normalise each column by min-max scaling
data6 = (data6 - data6.min()) / (data6.max() - data6.min())


# %%

# remove the y-axis labels
sns.clustermap(data6, cmap='RdBu_r', xticklabels=True, yticklabels=False)
#plt.savefig('heatmap_1142243F.png', dpi=500, bbox_inches='tight')

#%%
# plot all the luminal_breast_cancer data columns as one violin plot for the 6 different datasets

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['luminal_breast_cancer'], data2['luminal_breast_cancer'], data3['luminal_breast_cancer'], data4['luminal_breast_cancer'], data5['luminal_breast_cancer'], data6['luminal_breast_cancer']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('Luminal breast cancer')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()
# %%
# plot for basal_breast_cancer
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['basal_breast_cancer'], data2['basal_breast_cancer'], data3['basal_breast_cancer'], data4['basal_breast_cancer'], data5['basal_breast_cancer'], data6['basal_breast_cancer']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('Basal breast cancer')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()

# %%
# plot for epithelial_tumour_signature
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.boxplot(data=[data1['Epithelial_tumour_signature'], data2['Epithelial_tumour_signature'], data3['Epithelial_tumour_signature'], data4['Epithelial_tumour_signature'], data5['Epithelial_tumour_signature'], data6['Epithelial_tumour_signature']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('Epithelial tumour signature')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()
# %%
# plot for mesenchymal_tumour_signature
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['Mesenchymal_tumour_signature'], data2['Mesenchymal_tumour_signature'], data3['Mesenchymal_tumour_signature'], data4['Mesenchymal_tumour_signature'], data5['Mesenchymal_tumour_signature'], data6['Mesenchymal_tumour_signature']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('Mesenchymal tumour signature')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()
# %%
# plot for pEMT
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['pEMT'], data2['pEMT'], data3['pEMT'], data4['pEMT'], data5['pEMT'], data6['pEMT']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('pEMT')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()

# %%
# plot for HALLMARK_ESTROGEN_RESPONSE_EARLY
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['HALLMARK_ESTROGEN_RESPONSE_EARLY'], data2['HALLMARK_ESTROGEN_RESPONSE_EARLY'], data3['HALLMARK_ESTROGEN_RESPONSE_EARLY'], data4['HALLMARK_ESTROGEN_RESPONSE_EARLY'], data5['HALLMARK_ESTROGEN_RESPONSE_EARLY'], data6['HALLMARK_ESTROGEN_RESPONSE_EARLY']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('HALLMARK_ESTROGEN_RESPONSE_EARLY')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()
# %%
# plot violin plot for difference between mesenchymal and epithelial tumour signature
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['Mesenchymal_tumour_signature'] - data1['Epithelial_tumour_signature'], data2['Mesenchymal_tumour_signature'] - data2['Epithelial_tumour_signature'], data3['Mesenchymal_tumour_signature'] - data3['Epithelial_tumour_signature'], data4['Mesenchymal_tumour_signature'] - data4['Epithelial_tumour_signature'], data5['Mesenchymal_tumour_signature'] - data5['Epithelial_tumour_signature'], data6['Mesenchymal_tumour_signature'] - data6['Epithelial_tumour_signature']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('Mesenchymal - Epithelial tumour signature')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()

# %%
# plot violin plot for difference between basal and luminal breast cancer
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.violinplot(data=[data1['basal_breast_cancer'] - data1['luminal_breast_cancer'], data2['basal_breast_cancer'] - data2['luminal_breast_cancer'], data3['basal_breast_cancer'] - data3['luminal_breast_cancer'], data4['basal_breast_cancer'] - data4['luminal_breast_cancer'], data5['basal_breast_cancer'] - data5['luminal_breast_cancer'], data6['basal_breast_cancer'] - data6['luminal_breast_cancer']], ax=ax, palette=['lightblue', 'lightblue', '#FF7F7F', '#FF7F7F', '#FF7F7F', '#FF7F7F'], width=1.0)
ax.set_ylabel('Basal - Luminal breast cancer')
ax.set_xticklabels(['CID4290', 'CID4535', 'CID4465', 'CID44971', '1160920F', '1142243F'])
plt.show()

# %%
# write a fuction to return the probability distribution of an array of continuous data
def get_prob_dist(data):
    # get the number of bins
    n_bins = 20
    # get the histogram
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    # get the bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return hist, bin_centers

# merge prob dists from all the datasets
pathway2 = 'Epithelial_tumour_signature'
pathway1 = 'Mesenchymal_tumour_signature'
hist1, bin_centers1 = get_prob_dist(data1[pathway1]-data1[pathway2])
hist2, bin_centers2 = get_prob_dist(data2[pathway1]-data2[pathway2])
hist3, bin_centers3 = get_prob_dist(data3[pathway1]-data3[pathway2])
hist4, bin_centers4 = get_prob_dist(data4[pathway1]-data4[pathway2])
hist5, bin_centers5 = get_prob_dist(data5[pathway1]-data5[pathway2])
hist6, bin_centers6 = get_prob_dist(data6[pathway1]-data6[pathway2])

# merge the hist arrays into one numpy array
hist = np.array([hist1, hist2, hist3, hist4, hist5, hist6])
hist = hist/np.sum(hist, axis=1)[:, None]

print(hist.sum(axis=1))

print(entr(hist).sum(axis=1))

array = entr(hist).sum(axis=1)
labels = ['ER+', 'ER+', 'TNBC', 'TNBC', 'TNBC', 'TNBC']
# create dataframe
df_ent = pd.DataFrame({'ER status': labels, 'Entropy': array})
# plot barplot with error bars
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.barplot(data=df_ent, x='ER status', y='Entropy', ax=ax, palette=['lightblue', '#FF7F7F'])
ax.set_ylabel('Entropy')
plt.show()
# %%
pathway1 = 'Epithelial_tumour_signature'

hist1, bin_centers1 = get_prob_dist(data1[pathway1])
hist2, bin_centers2 = get_prob_dist(data2[pathway1])
hist3, bin_centers3 = get_prob_dist(data3[pathway1])
hist4, bin_centers4 = get_prob_dist(data4[pathway1])
hist5, bin_centers5 = get_prob_dist(data5[pathway1])
hist6, bin_centers6 = get_prob_dist(data6[pathway1])

# merge the hist arrays into one numpy array
hist = np.array([hist1, hist2, hist3, hist4, hist5, hist6])
hist = hist/np.sum(hist, axis=1)[:, None]

print(hist.sum(axis=1))

print(entr(hist).sum(axis=1))
array = entr(hist).sum(axis=1)
labels = ['ER+', 'ER+', 'TNBC', 'TNBC', 'TNBC', 'TNBC']
# create dataframe
df_ent = pd.DataFrame({'ER status': labels, 'Entropy': array})
# plot barplot with error bars
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.barplot(data=df_ent, x='ER status', y='Entropy', ax=ax, palette=['lightblue', '#FF7F7F'])
ax.set_ylabel('Entropy')
plt.show()

# %%
