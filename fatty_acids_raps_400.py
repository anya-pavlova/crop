import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%matplotlib inline

data_cam_df = pd.read_csv('xcms_centroid_rape_400_minfrac_0.5_camera.csv', sep = ',', index_col=0)
data_cam_df[data_cam_df==0]=np.nan

#Palm

new_mz_palm = np.abs((data_cam_df['mz'] - 271.2866)) / 271.2866 
new_mz_palm = new_mz_palm * 1000000
sunf_ppm_palm = data_cam_df[new_mz_palm<10]

sunf_ppm_palm.head()
sunf_ppm_palm.index

# Oleic compound

new_mz_oleic = np.abs((data_cam_df['mz'] - 299.3091)) / 299.3091 
new_mz_oleic = new_mz_oleic * 1000000
sunf_ppm_oleic = data_cam_df[new_mz_oleic<10]

sunf_ppm_oleic.index

#Staeric

new_mz_staeric = np.abs((data_cam_df['mz'] - 301.3246)) / 301.3246 
new_mz_staeric = new_mz_staeric * 1000000
sunf_ppm_staeric = data_cam_df[new_mz_staeric<10]

sunf_ppm_staeric.index

RO=pd.Series([i.split('_')[3][1:] if 'R' in i else i for i in data_cam_df.columns],index=data_cam_df.columns)

for i in RO.index:
    try:
        RO[i]=int(RO[i])
    except:
        RO[i]=np.nan
        
RO=RO.dropna()

#Oleic

plt.scatter(RO,data_cam_df.loc[351,RO.index], color = 'rebeccapurple')

#Palm
plt.scatter(RO,data_cam_df.loc[259,RO.index], color = 'rebeccapurple')

plt.scatter(RO,data_cam_df.loc[362,RO.index], color = 'rebeccapurple')

#filtration by rt

rt_in_min = data_cam_df['rt']/60
data_cam_df['rt'] = rt_in_min
data_cam_filt_rt_df = data_cam_df[(data_cam_df['rt'] > 0.6) & (data_cam_df['rt'] < 18)]
data_cam_filt_rt_df.head(4)

plt.gcf().set_size_inches(25,20)
plt.scatter(data_cam_filt_rt_df['rt'], data_cam_filt_rt_df['mz'], color = 'rebeccapurple', s = 50)

#filtering by isotopes from CAMERA¶

del_isotopes = data_cam_filt_rt_df['isotopes'].str.match(r'\[\d+\]\[M\+\d+\]\+').fillna(False)
#[m][M+n]+ где n от 1, m from 1
data_cam_filt_rt_iso_df = data_cam_filt_rt_df[~del_isotopes]

data_col_val = data_cam_filt_rt_iso_df.columns

col_with_blank = data_col_val.str.contains('blank')
blank_cal = data_col_val[col_with_blank]

name_blank_col = list(blank_cal)

peaks=data_cam_filt_rt_df.index[(data_cam_filt_rt_df['rt']<6.2)&(data_cam_filt_rt_df['rt']>5.8)]


tmp=data_cam_filt_rt_iso_df.copy().loc[:, 'X130418_AC_neg_blank':'X130418_AC_neg_R99_1.400']
tmp=np.log2(tmp)
tmp=tmp.T
plt.scatter(tmp.loc['X130418_AC_neg_blank'],tmp.mean(),marker='.')
plt.plot([8,30],[8,30],color='black')
blank_table = data_cam_filt_rt_iso_df[name_blank_col]

#choosing 1 batch

service_cols = ["index", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "0", "isotopes", "adduct", "pcgroup"]
samples_cols = list(set(data_cam_filt_rt_iso_df.columns) - set(service_cols))
data_df = data_cam_filt_rt_iso_df[samples_cols]

mean_data = data_df.transpose().mean()

plt.gcf().set_size_inches(15,10)
plt.scatter(np.log(blank_table['X130418_AC_neg_blank']), np.log(mean_data))
plt.scatter(np.log(blank_table['X130418_AC_neg_blank'][peaks]), np.log(mean_data[peaks]))
plt.plot([5,18],[5,18])
plt.xlabel('rt')
plt.ylabel('mz')

#PCA

service_cols = ["index", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "0", "isotopes", "adduct", "pcgroup"]
samples_cols = list(set(data_cam_filt_rt_iso_df.columns) - set(service_cols))
samples_filled_nan = data_cam_filt_rt_iso_df[samples_cols].fillna(0).transpose().as_matrix()

from sklearn.decomposition import PCA

pca = PCA(2)

transformed_samples = pca.fit_transform(samples_filled_nan)

plt.scatter(transformed_samples[:,0], transformed_samples[:,1])

#Annotation

annot_data = pd.read_csv('xcms_centroid_rape_400_minfrac_0.5_camera.tsv.ann.txt', sep = ',', index_col=0)
#loading data

lmfa_index = annot_data['lm_id'].str.contains('LMFA0101').fillna(False)
#строки в annot_data которые содержать LMFA0101 + что-то
lmfa_annot_data = annot_data[lmfa_index]
lmfa_annot_data = annot_data[lmfa_index]
#беру из annot_data строки которые соответствуют строкам в которых есть LMFA0101 + что-то


lmfa_index3 = annot_data['lm_id'].str.contains('LMFA0103').fillna(False)
#строки в annot_data которые содержать LMFA0101 + что-то
lmfa_annot_data3 = annot_data[lmfa_index3]
lmfa_annot_data3 = annot_data[lmfa_index3]
#беру из annot_data строки которые соответствуют строкам в которых есть LMFA0101 + что-то

index_lmfa_annot_data = lmfa_annot_data.index.unique()
#в index_lmfa_annot_data ищу уникальные индексы
index_lmfa_annot_data   

index_lmfa_annot_data3 = lmfa_annot_data3.index.unique()
#в index_lmfa_annot_data ищу уникальные индексы
index_lmfa_annot_data3

index_data_cam_filt_rt_iso_df = data_cam_filt_rt_iso_df.index
index_data_cam_filt_rt_iso_df = pd.DataFrame(index_data_cam_filt_rt_iso_df)

#беру индексы из data_cam_filt_rt_iso_df

lmfa_1_data = data_cam_filt_rt_iso_df.loc[index_lmfa_annot_data]
#из data_cam_filt_rt_iso_df беру индексы которые соответствуют index_lmfa_annot_data
not_lmfa_1_data = data_cam_filt_rt_iso_df[~data_cam_filt_rt_iso_df.index.isin(index_lmfa_annot_data)]

lmfa_3_data = data_cam_filt_rt_iso_df.loc[index_lmfa_annot_data3]
#из data_cam_filt_rt_iso_df беру индексы которые соответствуют index_lmfa_annot_data
not_lmfa_3_data = data_cam_filt_rt_iso_df[~data_cam_filt_rt_iso_df.index.isin(index_lmfa_annot_data3)]

plt.gcf().set_size_inches(25,20)
plt.scatter(lmfa_1_data['rt'], lmfa_1_data['mz'], c = 'blue')
plt.scatter(not_lmfa_1_data['rt'], not_lmfa_1_data['mz'], c = 'red')

plt.gcf().set_size_inches(25,20)
plt.scatter(lmfa_3_data['rt'], lmfa_3_data['mz'], c = 'blue')
plt.scatter(lmfa_1_data['rt'], lmfa_1_data['mz'], c = 'red')

fatty_acids = []

for n in range(-6,20):
    for m in range(0,6):
        fatty_acids.append((n+16,m, 255.2314435527 + n*14.01565 - m*2.01565))
fatty_acids = pd.DataFrame(fatty_acids, columns=["n","m","mass"])    

fatty_acids_annotation = {
    "index_fatty_acids" : [],
    "index_mz": [],
    "fatty_acid_mass": [],
    "mz": [],
    "m": [],
    "n": [],
}

for row in fatty_acids.itertuples():
    mz = data_cam_filt_rt_iso_df["mz"]
    m = row.mass
    condition = np.abs(mz-m)/m*1e6<10
    mz = mz[condition]
    
    k = len(mz)
    
    fatty_acids_annotation["index_fatty_acids"].extend([row.Index]*k)
    fatty_acids_annotation["index_mz"].extend(mz.index)
    fatty_acids_annotation["fatty_acid_mass"].extend([m]*k)
    fatty_acids_annotation["mz"].extend(mz)
    fatty_acids_annotation["m"].extend([row.m]*k)
    fatty_acids_annotation["n"].extend([row.n]*k)
    
fatty_acids_annotation_df = pd.DataFrame(fatty_acids_annotation)



fatty_acids_annotation_df["rt"] = data_cam_filt_rt_iso_df.loc[fatty_acids_annotation_df.index_mz, "rt"].values
#plt.gcf().set_size_inches(10,15)
plt.scatter(fatty_acids_annotation_df["rt"], fatty_acids_annotation_df["mz"])

plt.gcf().set_size_inches(15,15)

rt = fatty_acids_annotation_df["rt"] 
mz = fatty_acids_annotation_df["mz"]
n = fatty_acids_annotation_df["n"]
m = fatty_acids_annotation_df["m"]

plt.scatter(rt, mz, color = 'rebeccapurple', s = 50)

for rt_i, mz_i, n_i, m_i in zip(rt, mz, n, m):
    plt.annotate("{0}:{1}".format(n_i, m_i), (rt_i, mz_i))

fitered_fatty_acids_annot_df = fatty_acids_annotation_df[fatty_acids_annotation_df.rt<=7.81]
fitered_fatty_acids_annot_df = fitered_fatty_acids_annot_df.drop([40])

plt.gcf().set_size_inches(15,15)

rt = fitered_fatty_acids_annot_df["rt"] 
mz = fitered_fatty_acids_annot_df["mz"]
n = fitered_fatty_acids_annot_df["n"]
m = fitered_fatty_acids_annot_df["m"]

plt.scatter(rt, mz,  color = 'rebeccapurple', s = 50)

for rt_i, mz_i, n_i, m_i in zip(rt, mz, n, m):
    plt.annotate("{0}:{1}".format(n_i, m_i), (rt_i, mz_i))
    
fatty_acids_results = data_cam_filt_rt_iso_df.loc[fitered_fatty_acids_annot_df['index_mz']]

fatty_acids_results.to_csv('fatty_acids_results_rape_400.csv')

fatty_acids_annotation_df

columns_no_blanks = [col for col in samples_cols if 'blank' not in col]
mean_data_no_blanks = data_cam_filt_rt_iso_df[columns_no_blanks].transpose().mean()

fatty_acids_annotation_df["n:m"] = fatty_acids_annotation_df.n.astype(str) + ":" + fatty_acids_annotation_df.m.astype(str)

fatty_acids_even = fatty_acids_annotation_df[fatty_acids_annotation_df.n%2==0]
fatty_acids_odd = fatty_acids_annotation_df[fatty_acids_annotation_df.n%2==1]
fatty_acids_zero_m = fatty_acids_annotation_df[fatty_acids_annotation_df.m==0]


plt.gcf().set_size_inches(20,10)
plt.xticks(fatty_acids_zero_m["fatty_acid_mass"], fatty_acids_zero_m["n:m"])
plt.scatter(fatty_acids_even.fatty_acid_mass, np.log(mean_data_no_blanks.loc[fatty_acids_even.index_mz]), s =100, c='blueviolet', label='Even N')
plt.scatter(fatty_acids_odd.fatty_acid_mass, np.log(mean_data_no_blanks.loc[fatty_acids_odd.index_mz]),s = 100, c='darkred', label='Odd N')
plt.legend()

for fa_mass, mean_intense, nm in zip(fatty_acids_annotation_df["fatty_acid_mass"],
                                     np.log(mean_data_no_blanks.loc[fatty_acids_annotation_df.index_mz]),
                                     fatty_acids_annotation_df["n:m"]):
    plt.annotate(nm, (fa_mass, mean_intense))


#Normalization
service_cols = ["index", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "0", "isotopes", "adduct", "pcgroup"]
samples_cols = list(set(data_cam_filt_rt_iso_df.columns) - set(service_cols))

fatty_acids_results_no_blank = fatty_acids_results[[col for  col in samples_cols if 'blank' not in col]]

fatty_acids_results_no_blank_normal = np.log(fatty_acids_results_no_blank) - np.log(fatty_acids_results_no_blank).quantile(0.75)

fatty_acids_results_no_blank_normal_plus = fatty_acids_results_no_blank_normal + np.mean(np.log(fatty_acids_results_no_blank.quantile(0.75)))

fatty_acids_results_no_blank_normal_plus.to_csv('fatty_acids_results_no_blank_normal_plus_rape_400.csv')
fatty_acids_results.to_csv('fatty_acids_results_rape_400.csv')

#MDS

corr_fatty_acids = fatty_acids_results_no_blank_normal_plus.iloc[:, 3:].corr()
corr_fatty_acids.shape

from sklearn.manifold import MDS

pos = MDS(dissimilarity='precomputed', max_iter=3000, random_state=13).fit_transform(1-corr_fatty_acids)
pos = pos.T
pos = pd.DataFrame(pos, columns = corr_fatty_acids.columns)


qc_cols = [col for col in pos.columns if 'QC' in col]
non_qc_cols = [col for col in pos.columns if 'QC' not in col]

plt.gcf().set_size_inches(10,5)
plt.scatter(pos[non_qc_cols].iloc[0], pos[non_qc_cols].iloc[1], label='Non QC', s=100, facecolors='none', edgecolors='rebeccapurple')
plt.scatter(pos[qc_cols].iloc[0], pos[qc_cols].iloc[1], label='QC', s=100, facecolors='none', edgecolors='black')
plt.legend()

#MDS not normalization

fatty_acids_results_no_blank = np.log(fatty_acids_results_no_blank)
corr_fatty_acids = fatty_acids_results_no_blank.iloc[:, 3:].corr()
corr_fatty_acids.shape

from sklearn.manifold import MDS

pos = MDS(dissimilarity='precomputed', max_iter=3000, random_state=13).fit_transform(1-corr_fatty_acids)
pos = pos.T
pos = pd.DataFrame(pos, columns = corr_fatty_acids.columns)

qc_cols = [col for col in pos.columns if 'QC' in col]
non_qc_cols = [col for col in pos.columns if 'QC' not in col]

plt.gcf().set_size_inches(10, 5)
plt.scatter(pos[non_qc_cols].iloc[0], pos[non_qc_cols].iloc[1], label='Non QC', s=100, facecolors='none', edgecolors='blueviolet')
plt.scatter(pos[qc_cols].iloc[0], pos[qc_cols].iloc[1], label='QC', s=100, facecolors='none', edgecolors='black')
plt.legend()
