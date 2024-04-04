#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import product, starmap
# %%
def jaccard_similarity(dictionary_of_lists):
    keys = list(dictionary_of_lists.keys())
    jaccard_indices = {}
    for key1, list1 in dictionary_of_lists.items():
        for key2, list2 in dictionary_of_lists.items():
            if key1 != key2:
                intersection = len(set(list1).intersection(list2))
                union = len(set(list1).union(list2))
                jaccard_indices[(key1, key2)] = intersection / union

    matrix_data = []
    for key1 in dictionary_of_lists.keys():
        row = []
        for key2 in dictionary_of_lists.keys():
            if key1 == key2:
                row.append(1.0)  # Jaccard index of a list with itself is 1
            else:
                row.append(jaccard_indices[(key1, key2)])
        matrix_data.append(row)

    sns.set(font_scale=1.2)
    sns.clustermap(matrix_data, cmap="Blues", linewidths=0.5, figsize=(8, 8),vmin=0,vmax=1,annot=True, xticklabels=keys, yticklabels=keys)
    plt.show()
#%%
def common_genes_across_all_patients(data):
    first_list = list(data.values())[0]

    common_elements = []
    for element in first_list:
        is_common = all(element in sublist for sublist in data.values())
        if is_common:
            common_elements.append(element)

    return common_elements
# %%
patient_ids = ["CID3941","CID4067","CID4290A","CID4463","CID4530N"]
#patient_ids = ["CID3963","CID4465","CID4495","CID44971","CID4513","CID4515","CID4523"]
source_cell_type = ["T cells CD8+"]
target_cell_type = ["Cancer LumA SC"]
ligand_lists = {}
receptor_lists = {}
patients_considered = []
for patient_id in patient_ids:

    df_metadata = pd.read_csv("./patient_data/"+patient_id+"/metadata.csv",index_col=0)
    df_liana = pd.read_csv("./patient_data/"+patient_id+"/output_data_liana.tsv",sep="\t",index_col=0)
    #source_frq = df_metadata.groupby("celltype_minor").count()["celltype_major"][source_cell_type[0]]
    target_frq = df_metadata.groupby("celltype_minor").count()["celltype_major"][target_cell_type[0]]

    if target_frq < 100: #source_frq < 100 or
        continue
    df_liana_subset = df_liana[(df_liana["target"].isin(target_cell_type)) ]  #& (df_liana["source"].isin(source_cell_type))
    df_selected = df_liana_subset[(df_liana_subset["cellphonedb.pvalue"] < 0.05) | (df_liana_subset["sca.LRscore"] > 0.8)]
    ligand_lists[patient_id] = set(df_selected["ligand.complex"])
    receptor_lists[patient_id] = set(df_selected["receptor.complex"])
    patients_considered.append(patient_id)


# %%
jaccard_similarity(receptor_lists)
# %%
common_receptors = common_genes_across_all_patients(receptor_lists)
print(len(common_receptors))
df_targets = set([])
for i in common_receptors:
    if "_" in i:
        df_targets.add(i.split("_")[0])
        df_targets.add(i.split("_")[1])
    else:
        df_targets.add(i)
#%%
with open("luminal_receptors.txt","w") as f:
    for gene in df_targets:
        f.write(gene+"\n")

# %%
for patient_id in patient_ids:
    df_metadata = pd.read_csv("./patient_data/"+patient_id+"/metadata.csv",index_col=0)

    target_frq = df_metadata.groupby("celltype_minor").count()["celltype_major"][target_cell_type[0]]

    if target_frq < 100:
        continue

    df_AUCell_imputed = pd.read_csv("./patient_data/"+patient_id+"/output_AUCellscores_imputed.tsv",sep="\t",index_col=0).T
    df_AUCell_not_imputed = pd.read_csv("./patient_data/"+patient_id+"/output_AUCellscores_not_imputed.tsv",sep="\t",index_col=0).T
    df_liana = pd.read_csv("./patient_data/"+patient_id+"/output_data_liana.tsv",sep="\t",index_col=0)
    df_imputed = pd.read_csv("./patient_data/"+patient_id+"/output_imputed_values.tsv",sep="\t",index_col=0)

    df_concat = pd.concat([df_imputed,df_AUCell_imputed,df_metadata],axis=1)

    df_subset = df_concat[df_concat["celltype_minor"].isin(target_cell_type)]
    sns.set(font_scale=0.5)
    sns.clustermap(df_subset[list(df_targets)+['Epithelial_tumour_signature','Mesenchymal_tumour_signature','pEMT','luminal_breast_cancer','basal_breast_cancer']].corr(method="spearman"),vmin=-1,vmax=1,cmap="RdBu_r", figsize=(10, 10))
    plt.show()
    plt.close()
# %%
df_metaanalysis = pd.read_csv("meta_analysis_results.txt",index_col=0)
df_labels_genes = pd.read_csv("all_ligand_receptors_combined.txt",sep='\t',index_col=0)
# %%
lut = dict(zip(df_labels_genes['status'].unique(), ['red','green','blue','pink','black','yellow']))
row_colors = df_labels_genes['status'].map(lut)
sns.clustermap(df_metaanalysis,cmap="RdBu_r",row_colors=row_colors,col_cluster=False)
# %%
