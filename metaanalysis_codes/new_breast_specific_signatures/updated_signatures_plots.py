#%%
import os
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams

rcParams.update({'font.family':'Arial','font.size':18})
# %%
f = plt.figure()
f.set_figwidth(5)
f.set_figheight(5)

filename = 'basal_vs_clusterC_t.txt'

# read "luminal_vs_clusterA_t.txt" to a dataframe without quotes
df = pd.read_csv('./data_luminal_basal_mam_sp/'+filename, sep=' ', quotechar='"', index_col=0, header=0)

corr1 = list(df['Correlation'])
corr2 = list(df['Pvalue'])

c0 = c1 = c2 = 0
for i,j in enumerate(corr1):
    c0 += 1
    if j >= 0.3 and corr2[i] <= 0.05:
        c1 += 1
        plt.scatter([j],-1*np.log10([corr2[i]]),c='red',s=20)
    elif j <= -0.3 and corr2[i] <= 0.05:
        c2 += 1
        plt.scatter([j],-1*np.log10([corr2[i]]),c='blue',s=20)
    else:
        plt.scatter([j],-1*np.log10([corr2[i]]),c='black',s=20)

print(c1,c2,c0)
plt.xlim([-1.1,1.1])
plt.axhline(y = 1.301,c='grey')
plt.axvline(x = 0.3,c='grey')
plt.axvline(x = 0,c='grey')
plt.axvline(x = -0.3,c='grey')

# save the plot
plt.savefig('./figures_luminal_basal_mam_sp/'+filename[:-4]+'.png',bbox_inches='tight',dpi=800)
# %%

filename = 'basal_vs_clusterC_t.txt'
df = pd.read_csv('./data_luminal_basal_mam_sp/'+filename, sep=' ', quotechar='"', index_col=0, header=0)

corr1 = list(df['Correlation'])
corr2 = list(df['Pvalue'])

filename = 'basal_vs_clusterB_t.txt'
df = pd.read_csv('./data_luminal_basal_mam_sp/'+filename, sep=' ', quotechar='"', index_col=0, header=0)

corr3 = list(df['Correlation'])
corr4 = list(df['Pvalue'])

afters = []
mids = []
for i,j in enumerate(corr2):
    if j < 0.05 or corr4[i] < 0.05:
        afters.append(corr1[i])
        mids.append(corr3[i])

print(ss.ttest_rel(mids, afters))


df = pd.DataFrame({'condition 1': afters,'condition 2': mids})

jitter = 0.05
df_x_jitter = pd.DataFrame(np.random.normal(loc=0, scale=jitter, size=df.values.shape), columns=df.columns)
df_x_jitter += np.arange(len(df.columns))

fig, ax = plt.subplots()
for col in df:
    if col[-1] == '1':
        ax.plot(df_x_jitter[col], df[col], 'o', alpha=.9, zorder=1, ms=6, mew=1, color='red')
    if col[-1] == '2':
        ax.plot(df_x_jitter[col], df[col], 'o', alpha=.9, zorder=1, ms=6, mew=1, color='orange')
ax.set_xticks(range(len(df.columns)))
ax.set_xticklabels(df.columns)
ax.set_xlim(-0.5,len(df.columns)-0.5)

for idx in df.index:
    ax.plot(df_x_jitter.loc[idx,['condition 1','condition 2']], df.loc[idx,['condition 1','condition 2']], color = 'grey', linewidth = 0.5, linestyle = '-', zorder=-1)

plt.savefig('./figures_luminal_basal_mam_sp/'+filename[:-4]+'_comp.png',bbox_inches='tight',dpi=800)

# %%
