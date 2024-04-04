#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18,'font.family': 'arial'})
# %%
df_tgfb1 = pd.read_csv("data_tgfb1.txt",sep="\t")
# %%
sns.barplot(x="x",y="y",data=df_tgfb1,hue="data",palette=["grey","red"],order=["Lum","Lum_H","Bas_epi","Bas_H","Bas_mes"],capsize=.1)
plt.legend([],frameon=False)
plt.ylim([0,0.5])
plt.savefig("tgfb.png",dpi=800)
# %%
df_tgfb1 = pd.read_csv("data_il1b.txt",sep="\t")
sns.barplot(x="x",y="y",data=df_tgfb1,hue="data",palette=["grey","orange"],order=["Lum","Lum_H","Bas_epi","Bas_H","Bas_mes"],capsize=.1)
plt.legend([],frameon=False)
plt.ylim([0,0.4])
plt.savefig("il1b.png",dpi=800)

# %%
df_tgfb1 = pd.read_csv("data_percent.txt",sep="\t")
sns.barplot(x="x",y="y",data=df_tgfb1,hue="data",palette=["orange","red"],order=["Lum","Lum_H","Bas_epi","Bas_H","Bas_mes"],capsize=.1)
plt.legend([],frameon=False)
plt.ylim([-55,100])
plt.savefig("percent.png",dpi=800)
# %%
