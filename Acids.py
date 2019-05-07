# С16 H32 O2 пальмитиновая C16 H30 O2 пальмитолеиновая CH3(CH2)16 COOH стеариновая C18 H34 O2 олеиновая (282.467) C18 H32 O2 линолевая
# C18 H30 O2 линоленовая C20 H40 O2 арахиновая C20 H38 О2 экозеновая (310.522) C22 H44 O2 бегеновая
# C22 H42 O2 эруковая (338.575) C24 H42 O2 лигноцериновая (368.646) C24 H46 O2 селахолевая (366.63)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_cam_df = pd.read_csv('xs_annotated_test_sunflower.csv', sep = ',')

data_cam_df[data_cam_df==0]=np.nan

new_mz_palm = np.abs((data_cam_df['mz'] - 529.3993)) / 529.3993
new_mz_palm = new_mz_palm * 1000000
sunf_ppm_palm = data_cam_df[new_mz_palm<10]

#Palm

new_mz_palm = np.abs((data_cam_df['mz'] - 271.2866)) / 271.2866
new_mz_palm = new_mz_palm * 1000000
sunf_ppm_palm = data_cam_df[new_mz_palm<10]

sunf_ppm_palm.head()
sunf_ppm_palm['index']

#Oleic compound

new_mz_oleic = np.abs((data_cam_df['mz'] - 299.3091)) / 299.3091
new_mz_oleic = new_mz_oleic * 1000000
sunf_ppm_oleic = data_cam_df[new_mz_oleic<10]

sunf_ppm_oleic['index']

#Staeric
new_mz_staeric = np.abs((data_cam_df['mz'] - 301.3246)) / 301.3246
new_mz_staeric = new_mz_staeric * 1000000
sunf_ppm_staeric = data_cam_df[new_mz_staeric<10]

sunf_ppm_staeric['index']

RO=pd.Series([i.split('_')[3][1:] if 'S' in i else i for i in data_cam_df.columns],index=data_cam_df.columns)

for i in RO.index:
    try:
        RO[i]=int(RO[i])
    except:
        RO[i]=np.nan
RO=RO.dropna()

plt.scatter(RO,data_cam_df.loc[265,RO.index])
plt.scatter(RO,data_cam_df.loc[259,RO.index])
plt.scatter(RO,data_cam_df.loc[106,RO.index])

#filtration by rt
#rt(sec) -> rt(min)

a = data_cam_df['rt']/60
data_cam_df['rt'] = a
data_cam_filt_rt_df = data_cam_df[(data_cam_df['rt'] > 0.6) & (data_cam_df['rt'] < 18)]
data_cam_filt_rt_df.head(4)

#filtering by isotopes from CAMERA

del_isotopes = data_cam_filt_rt_df['isotopes'].str.match(r'\[\d+\]\[M\+\d+\]\+').fillna(False)
#[m][M+n]+ где n от 1, m from 1
data_cam_filt_rt_iso_df = data_cam_filt_rt_df[-del_isotopes]

peaks=data_cam_filt_rt_df.index[(data_cam_filt_rt_df['rt']<6.2)&(data_cam_filt_rt_df['rt']>5.8)]
tmp=data_cam_filt_rt_df.copy().loc[:,'X260318_AC_neg_ACN':'X260318_AC_neg_Sblank2_1.5']
tmp=np.log2(tmp)
tmp=tmp.T

plt.hist(tmp.mean(),bins=100);
peaks=tmp.columns[tmp.mean()>14]

plt.gcf().set_size_inches(25,20)
plt.scatter(data_cam_filt_rt_df['rt'], data_cam_filt_rt_df['mz'])

plt.gcf().set_size_inches(25,20)
plt.scatter(data_cam_filt_rt_df['rt'], data_cam_filt_rt_df['mz'],c=tmp.mean(),cmap='binary')
plt.scatter(data_cam_filt_rt_df.loc[peaks,'rt'], data_cam_filt_rt_df.loc[peaks,'mz'])
plt.xlabel('rt')
plt.ylabel('mz')

data_col_val = data_cam_filt_rt_iso_df.columns

if isinstance(data_col_val, np.ndarray):
    data_col_val = pd.DataFrame(a)
type(a)


data_cam_filt_rt_iso_df.columns

col_with_blank = data_col_val.str.contains('blank')
#col_with_blank = data_col_val[0].str.match(r'blank')
blank_cal = data_col_val[col_with_blank]

#data_col_val
name_blank_col = list(blank_cal)

name_blank_col

import numpy as np
tmp=data_cam_filt_rt_iso_df.copy().loc[:,'X260318_AC_neg_ACN':'X260318_AC_neg_Sblank2_1.5']
tmp=np.log2(tmp)
tmp=tmp.T
plt.scatter(tmp.loc['X260318_AC_neg_Sblank2_1.5'],tmp.mean(),marker='.')
plt.plot([8,30],[8,30],color='black')

blank_table = data_cam_filt_rt_iso_df[name_blank_col]
blank_table.head(10)
blank_table['X260318_AC_neg_Sblank2_1.5'].mean() #средняя интенсивность последнего блэнка

#choosing 7 batch

batch7 = pd.read_csv('batch7_samples.txt')
batch7 = list(batch7['sample'])
#batch7 = ['X'+s[:-6] for s in batch7]

res = []
for s in batch7:
    s = 'X'+s[:-6]
    s = s.replace("1-5", "1.5")
    res.append(s)
batch7 = res

data_cam_filt_rt_iso_df.head(5)

batch7_columns = list(set(data_cam_filt_rt_iso_df.columns) & set(batch7))
data_batch7 = data_cam_filt_rt_iso_df[batch7_columns]


trans_data_batch7 = data_batch7.transpose()
mean_data_batch7 = trans_data_batch7.mean()

blank_table['X260318_AC_neg_Sblank2_1.5'] #средняя интенсивность последнего блэнка

plt.gcf().set_size_inches(15,10)
plt.scatter(np.log(blank_table['X260318_AC_neg_Sblank2_1.5']), np.log(mean_data_batch7))
plt.scatter(np.log(blank_table['X260318_AC_neg_Sblank2_1.5'][peaks]), np.log(mean_data_batch7[peaks]))
plt.plot([5,18],[5,18])
plt.xlabel('rt')
plt.ylabel('mz')

service_cols = ["index", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "0", "isotopes", "adduct", "pcgroup"]
samples_cols = list(set(data_cam_filt_rt_iso_df.columns) - set(service_cols))
samples_filled_nan = data_cam_filt_rt_iso_df[samples_cols].fillna(0).transpose().as_matrix()

from sklearn.decomposition import PCA

pca = PCA(2)

transformed_samples = pca.fit_transform(samples_filled_nan)

plt.scatter(transformed_samples[:,0], transformed_samples[:,1])

transformed_samples_df = pd.DataFrame(transformed_samples.T, columns = samples_cols)
transformed_samples_df

non_batch7_columns = list(set(samples_cols) - set(batch7_columns))

plt.scatter(transformed_samples[:,0], transformed_samples[:,1], c='blue')
plt.scatter(transformed_samples_df[batch7_columns].iloc[0], transformed_samples_df[batch7_columns].iloc[1], c='red')

plt.scatter(transformed_samples_df[non_batch7_columns].iloc[0], transformed_samples_df[non_batch7_columns].iloc[1],
            c='blue')
plt.scatter(transformed_samples_df[batch7_columns].iloc[0], transformed_samples_df[batch7_columns].iloc[1], c='red')

# annotation
annot_data = pd.read_csv('xs_annotated_test_sunflower.ann.csv', sep=',')
# loading data

annot_data.set_index('index').loc[[35, 52, 68, 99, 187, 189, 248, 265, 297, 319, 346,
                                   360, 361, 395, 440, 455, 484, 536, 552, 553, 584, 704,
                                   981, 1136, 1377, 194]]

lmfa_index = annot_data['lm_id'].str.contains('LMFA0101').fillna(False)
# строки в annot_data которые содержать LMFA0101 + что-то
lmfa_annot_data = annot_data[lmfa_index]
lmfa_annot_data = annot_data[lmfa_index]
# беру из annot_data строки которые соответствуют строкам в которых есть LMFA0101 + что-то

lmfa_index3 = annot_data['lm_id'].str.contains('LMFA0103').fillna(False)
# строки в annot_data которые содержать LMFA0101 + что-то
lmfa_annot_data3 = annot_data[lmfa_index3]
lmfa_annot_data3 = annot_data[lmfa_index3]
# беру из annot_data строки которые соответствуют строкам в которых есть LMFA0101 + что-то

index_lmfa_annot_data = lmfa_annot_data['index'].unique()
# в index_lmfa_annot_data ищу уникальные индексы
index_lmfa_annot_data

index_lmfa_annot_data3 = lmfa_annot_data3['index'].unique()
# в index_lmfa_annot_data ищу уникальные индексы
index_lmfa_annot_data3

index_data_cam_filt_rt_iso_df = data_cam_filt_rt_iso_df['index']
index_data_cam_filt_rt_iso_df = pd.DataFrame(index_data_cam_filt_rt_iso_df)

# беру индексы из data_cam_filt_rt_iso_df

index_data_cam_filt_rt_iso_df
index_lmfa_annot_data
index_lmfa_annot_data3

data_cam_filt_rt_iso_df.set_index('index', inplace=True)

b = data_cam_filt_rt_iso_df.loc[index_lmfa_annot_data]
# из data_cam_filt_rt_iso_df беру индексы которые соответствуют index_lmfa_annot_data
a = data_cam_filt_rt_iso_df[~data_cam_filt_rt_iso_df.index.isin(index_lmfa_annot_data)]

c = data_cam_filt_rt_iso_df.loc[index_lmfa_annot_data3]
# из data_cam_filt_rt_iso_df беру индексы которые соответствуют index_lmfa_annot_data
d = data_cam_filt_rt_iso_df[~data_cam_filt_rt_iso_df.index.isin(index_lmfa_annot_data3)]

plt.gcf().set_size_inches(25, 20)
plt.scatter(b['rt'], b['mz'], c='blue')
plt.scatter(a['rt'], a['mz'], c='red')

plt.gcf().set_size_inches(25, 20)
plt.scatter(c['rt'], c['mz'], c='blue')
plt.scatter(b['rt'], b['mz'], c='red')

annot_data[annot_data['lm_id'] == 'LMFA01010001']

data_cam_filt_rt_iso_df.loc[[68, 265, 319]]

# Fatty acids
fatty_acids = []

for n in range(-6, 20):
    for m in range(0, 6):
        fatty_acids.append((n + 16, m, 255.2314435527 + n * 14.01565 - m * 2.01565))
fatty_acids = pd.DataFrame(fatty_acids, columns=["n", "m", "mass"])

fatty_acids_annotation = {
    "index_fatty_acids": [],
    "index_mz": [],
    "fatty_acid_mass": [],
    "mz": [],
    "m": [],
    "n": [],
}

for row in fatty_acids.itertuples():
    mz = data_cam_filt_rt_iso_df["mz"]
    m = row.mass
    condition = np.abs(mz - m) / m * 1e6 < 10
    mz = mz[condition]

    k = len(mz)

    fatty_acids_annotation["index_fatty_acids"].extend([row.Index] * k)
    fatty_acids_annotation["index_mz"].extend(mz.index)
    fatty_acids_annotation["fatty_acid_mass"].extend([m] * k)
    fatty_acids_annotation["mz"].extend(mz)
    fatty_acids_annotation["m"].extend([row.m] * k)
    fatty_acids_annotation["n"].extend([row.n] * k)

fatty_acids_annotation_df = pd.DataFrame(fatty_acids_annotation)

fatty_acids_annotation_df["rt"] = data_cam_filt_rt_iso_df.loc[fatty_acids_annotation_df.index_mz, "rt"].values

# plt.gcf().set_size_inches(10,15)
plt.scatter(fatty_acids_annotation_df["rt"], fatty_acids_annotation_df["mz"])

plt.gcf().set_size_inches(15, 15)

rt = fatty_acids_annotation_df["rt"]
mz = fatty_acids_annotation_df["mz"]
n = fatty_acids_annotation_df["n"]
m = fatty_acids_annotation_df["m"]

plt.scatter(rt, mz, color='rebeccapurple', s=50)

for rt_i, mz_i, n_i, m_i in zip(rt, mz, n, m):
    plt.annotate("{0}:{1}".format(n_i, m_i), (rt_i, mz_i))

fatty_acids_annotation_df[(fatty_acids_annotation_df['n'] == 18) & (fatty_acids_annotation_df['m'] == 1)].sort_values(
    'rt')

fatty_acids_annotation_df[(fatty_acids_annotation_df['n'] == 27) & (fatty_acids_annotation_df['m'] == 5)].sort_values(
    'rt')

fatty_acids_annotation_df[(fatty_acids_annotation_df['n'] == 18) & (fatty_acids_annotation_df['m'] == 3)].sort_values(
    'rt')

fatty_acids_annotation_df[(fatty_acids_annotation_df['n'] == 28) & (fatty_acids_annotation_df['m'] == 0)].sort_values(
    'rt')

fatty_acids_annotation_df[(fatty_acids_annotation_df['n'] == 28) & (fatty_acids_annotation_df['m'] == 0)].sort_values(
    'rt')

through_out_inds = [20, 19, 23, 22, 15, 18, 17, 14, 21, 13, 12, 16, 43, 25, 27]

fitered_fatty_acids_annot_df = fatty_acids_annotation_df[~fatty_acids_annotation_df.index.isin(through_out_inds)]

plt.gcf().set_size_inches(15, 15)

rt = fitered_fatty_acids_annot_df["rt"]
mz = fitered_fatty_acids_annot_df["mz"]
n = fitered_fatty_acids_annot_df["n"]
m = fitered_fatty_acids_annot_df["m"]

plt.scatter(rt, mz, color='rebeccapurple', s=50)

for rt_i, mz_i, n_i, m_i in zip(rt, mz, n, m):
    plt.annotate("{0}:{1}".format(n_i, m_i), (rt_i, mz_i))

fitered_fatty_acids_annot_df.shape

fatty_acids_results = data_cam_filt_rt_iso_df.loc[fitered_fatty_acids_annot_df['index_mz']]
fatty_acids_results.to_csv('fatty_acids_results.csv')
# Draw blanks
batch7_columns_no_blanks = [col for col in batch7_columns if 'blank' not in col]
fatty_acids_results_batch7 = fatty_acids_results[batch7_columns_no_blanks]
fatty_acids_peaks = fatty_acids_results.index[(fatty_acids_results['rt'] < 0.2) & (fatty_acids_results['rt'] > 5.8)]

fatty_acids_peaks

tmp = fatty_acids_results.copy().loc[:, 'X260318_AC_neg_ACN':'X260318_AC_neg_Sblank2_1.5']
tmp = np.log2(tmp)
tmp = tmp.T

plt.hist(tmp.mean(), bins=100);
fatty_acids_peaks = tmp.columns[tmp.mean() > 14]

fatty_acids_peaks

fatty_acids_results.index

mean_fatty_acids_batch7 = fatty_acids_results_batch7.transpose().mean()

plt.gcf().set_size_inches(15, 10)
plt.scatter(np.log(fatty_acids_results['X260318_AC_neg_Sblank2_1.5']), np.log(mean_fatty_acids_batch7))
# plt.scatter(np.log(fatty_acids_results['X260318_AC_neg_Sblank2_1.5'][fatty_acids_results.index]), np.log(mean_fatty_acids_batch7.loc[fatty_acids_results.index]))
plt.plot([5, 18], [5, 18])
plt.xlabel('rt')
plt.ylabel('mz')

# Add other peaks
fitered_fatty_acids_annot_df = fatty_acids_annotation_df[~fatty_acids_annotation_df.index.isin(through_out_inds)]

# Add other peaks from batch 7
data_batch7 = data_cam_filt_rt_iso_df[batch7_columns_no_blanks]

data_batch7_last = data_batch7[~data_batch7.index.isin(fatty_acids_annotation_df['index_mz'])]
data_batch7_blank_last = data_cam_filt_rt_iso_df['X260318_AC_neg_Sblank2_1.5'][
    ~data_batch7.index.isin(fatty_acids_annotation_df['index_mz'])]

mean_data_batch7_last = data_batch7_last.transpose().mean()

plt.gcf().set_size_inches(15, 10)
plt.scatter(np.log(fatty_acids_results['X260318_AC_neg_Sblank2_1.5']), np.log(mean_fatty_acids_batch7))
# plt.scatter(np.log(fatty_acids_results['X260318_AC_neg_Sblank2_1.5'][fatty_acids_results.index]), np.log(mean_fatty_acids_batch7.loc[fatty_acids_results.index]))
plt.scatter(np.log(data_batch7_blank_last), np.log(mean_data_batch7_last))
plt.plot([5, 18], [5, 18])
plt.xlabel('rt')
plt.ylabel('mz')

columns_no_blanks = [col for col in samples_cols if 'blank' not in col]
mean_data_no_blanks = data_cam_filt_rt_iso_df[columns_no_blanks].transpose().mean()

fatty_acids_annotation_df["n:m"] = fatty_acids_annotation_df.n.astype(str) + ":" + fatty_acids_annotation_df.m.astype(
    str)

fatty_acids_even = fatty_acids_annotation_df[fatty_acids_annotation_df.n % 2 == 0]
fatty_acids_odd = fatty_acids_annotation_df[fatty_acids_annotation_df.n % 2 == 1]
fatty_acids_zero_m = fatty_acids_annotation_df[fatty_acids_annotation_df.m == 0]

plt.gcf().set_size_inches(20, 10)
plt.xticks(fatty_acids_zero_m["fatty_acid_mass"], fatty_acids_zero_m["n:m"])
plt.scatter(fatty_acids_even.fatty_acid_mass, mean_data_no_blanks.loc[fatty_acids_even.index_mz], s=100, c='blueviolet',
            label='Even N')
plt.scatter(fatty_acids_odd.fatty_acid_mass, mean_data_no_blanks.loc[fatty_acids_odd.index_mz], s=100, c='darkred',
            label='Odd N')
plt.legend()

for fa_mass, mean_intense, nm in zip(fatty_acids_annotation_df["fatty_acid_mass"],
                                     mean_data_no_blanks.loc[fatty_acids_annotation_df.index_mz],
                                     fatty_acids_annotation_df["n:m"]):
    plt.annotate(nm, (fa_mass, mean_intense))

# Normalisation

service_cols = ["index", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "0", "isotopes", "adduct", "pcgroup"]
samples_cols = list(set(data_cam_filt_rt_iso_df.columns) - set(service_cols))

fatty_acids_results_no_blank = fatty_acids_results[[col for col in samples_cols if 'blank' not in col]]
fatty_acids_results_no_blank_normal = np.log(fatty_acids_results_no_blank) - np.log(
    fatty_acids_results_no_blank).quantile(0.75)
fatty_acids_results_no_blank_normal_plus = fatty_acids_results_no_blank_normal + np.mean(
    np.log(fatty_acids_results_no_blank.quantile(0.75)))
fatty_acids_results_no_blank_normal_plus.to_csv('fatty_acids_results_no_blank_normal_plus.csv')
fatty_acids_results.to_csv('fatty_acids_results.csv')

# MDS

import re

batch_number_pattern = re.compile(r"_S(\d+)")


def get_extraction_number(name):
    search_result = batch_number_pattern.search(name)
    if search_result is None:
        return None

    num = int(search_result.group(1))
    return num


def get_batch_number(name):
    if "QC" in name or "ACN" in name:
        return 8

    num = get_extraction_number(name)
    assert num is not None

    limits = [97, 193, 276, 346, 441, 536]

    for i, lim in enumerate(limits, 1):
        if num <= lim:
            return i

    return len(limits) + 1  # 7


corr_fatty_acids = fatty_acids_results_no_blank_normal_plus.iloc[:, 3:].corr()

corr_fatty_acids.shape

from sklearn.manifold import MDS

pos = MDS(dissimilarity='precomputed', max_iter=3000, random_state=13).fit_transform(1 - corr_fatty_acids)
pos = pos.T
pos = pd.DataFrame(pos, columns=corr_fatty_acids.columns)

qc_cols = [col for col in pos.columns if 'QC' in col]
non_qc_cols = [col for col in pos.columns if 'QC' not in col]

plt.gcf().set_size_inches(15, 10)
plt.scatter(pos[non_qc_cols].iloc[0], pos[non_qc_cols].iloc[1], label='Non QC', s=100, facecolors='none',
            edgecolors='rebeccapurple')
plt.scatter(pos[qc_cols].iloc[0], pos[qc_cols].iloc[1], label='QC', s=100, facecolors='none', edgecolors='black')
plt.legend()

cols_by_batch = {}

for col in pos.columns:
    batch_num = get_batch_number(col)
    batch_cols = cols_by_batch.setdefault(batch_num, [])
    batch_cols.append(col)

plt.gcf().set_size_inches(15, 10)

# colors = ['b', 'g', 'r', 'c', 'm', 'y', 'purple', 'black']
colors = ["#FF0000", "#FFFF00", "#00FF00", "#D6BCC0", "#00FFFF",
          "#0000FF", "blueviolet", 'black']

for batch_num, batch_cols in sorted(cols_by_batch.items()):
    label = "Batch {0}".format(batch_num)
    plt.scatter(pos[batch_cols].iloc[0], pos[batch_cols].iloc[1], label=label, s=100, facecolors='none',
                edgecolors=colors[batch_num - 1])

plt.legend()

# MDS not normalization
fatty_acids_results_no_blank = np.log(fatty_acids_results_no_blank)
corr_fatty_acids = fatty_acids_results_no_blank.iloc[:, 3:].corr()
corr_fatty_acids.shape

from sklearn.manifold import MDS

pos = MDS(dissimilarity='precomputed', max_iter=3000, random_state=13).fit_transform(1 - corr_fatty_acids)
pos = pos.T
pos = pd.DataFrame(pos, columns=corr_fatty_acids.columns)

# other samples
qc_cols = [col for col in pos.columns if 'QC' in col]
non_qc_cols = [col for col in pos.columns if 'QC' not in col]

plt.gcf().set_size_inches(10, 5)
plt.scatter(pos[non_qc_cols].iloc[0], pos[non_qc_cols].iloc[1], label='Non QC', s=100, facecolors='none',
            edgecolors='blueviolet')
plt.scatter(pos[qc_cols].iloc[0], pos[qc_cols].iloc[1], label='QC', s=100, facecolors='none', edgecolors='black')
plt.legend()

need_samples = pd.read_csv('need_sunf_samples.csv', sep=',')

needed_columns = []

for col in fatty_acids_results_no_blank_normal_plus.columns:
    col_extraction_number = get_extraction_number(col)

    if (col_extraction_number is not None) and (col_extraction_number in need_samples.extraction_number):
        needed_columns.append(col)

fatty_acids_results_no_blank_normal_plus[needed_columns]
