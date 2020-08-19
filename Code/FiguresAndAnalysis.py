#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 08:22:24 2020

@author: prabh

If you run outside the main project folder change the CommonDf.tsv address in line 27
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import normaltest
pd.set_option('display.max_columns', 15)
pd.set_option('display.width', 200)
# sns.set(context='paper')
sns.set(context='notebook')
pd.set_option('display.max_columns', 15)
pd.set_option('display.width', 120)

CommonDf = pd.read_csv('../Tables/CommonDf.tsv', sep='\t') #13163
#drop gene families with no overlap
CommonDf.dropna(inplace=True, subset=['overlap']) #13160

AaLs = [x for x in list(CommonDf.columns) if x.endswith('_aa')]
CommonDf['Min'] = CommonDf[AaLs].apply(lambda row: min(row), axis=1)
CommonDf['Max'] = CommonDf[AaLs].apply(lambda row: max(row), axis=1)


################# figure 2
CommonDf['MaxDif'] = CommonDf[AaLs].apply(lambda row: (max(row) - min(row)), axis=1)
CommonDf['RelMaxDif'] = CommonDf['MaxDif']*100/CommonDf['Min']

CommonDf['RelOverlap'] = CommonDf['overlap']*100/CommonDf['Min']
CommonDf['%MinId'] = CommonDf['AbsId']*100/CommonDf['Min']

plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(nrows=1, ncols=4,figsize=(10,4))
#figure 2a scatterplot for min and max distribution of each cluster
sns.regplot(x="Min", y="Max", data=CommonDf, y_jitter=0.05, line_kws={'color':'red'}, ax=ax1)
ax1.set(xticks=range(0,8001,2000), xlim=(0,8000), ylim=(0,8000), xlabel="Shortest Ortholog\n(aa residues)", ylabel="Longest Ortholog (aa residues)")
plt.setp(ax1.get_xmajorticklabels(), rotation=30)

#figure 2b %maxdiff histogram
sns.distplot(CommonDf[CommonDf['RelMaxDif'] > 0]['RelMaxDif'], bins=range(0,15,1), kde=False, ax=ax2)
#counts, bins, patches = ax3.hist(CommonDf[CommonDf['%MaxDif'] > 0]['%MaxDif'], bins=10)
ax2.set(xticks=range(1,15,1), xlim=(0,9), xlabel="Length difference (%)", ylabel="Ortholog Families")
plt.setp(ax2.get_xmajorticklabels(), rotation=30)

# figure 2c Rel-overlap
sns.distplot(CommonDf['RelOverlap'], bins=range(90,101,1), kde=False, ax=ax3)
ax3.set(xticks=range(90,101,2), xlim=(90,100), xlabel="Aligment Saturation (%)", ylabel="Ortholog Families")
plt.setp(ax3.get_xmajorticklabels(), rotation=30)

# figure 2d %Abs identiy
sns.distplot(CommonDf['%AbsId'], bins=range(90,101,1), kde=False, ax=ax4)
ax4.set(xticks=range(90,101,2), xlim=(90,100), xlabel="Idenity (%)\n(within aligned region)", ylabel="Ortholog Families")
plt.setp(ax4.get_xmajorticklabels(), rotation=30)

plt.tight_layout(pad=1.08)
plt.savefig('Fig2_OrthologFamilies.pdf')


#Min in sp
for Index, Row in CommonDf.iterrows():
    m = 0
    for Sp in AaLs:
        if Row[Sp] == Row['Min']:
            CommonDf.at[Index, 'Min_sp'] = Sp
            m += 1
    if m > 1: CommonDf.at[Index, 'Min_sp'] = m
CommonDf.Min_sp.value_counts()

for Index, Row in CommonDf.iterrows():
    m = 0
    for Sp in AaLs:
        if Row[Sp] == Row['Max']:
            CommonDf.at[Index, 'Max_sp'] = Sp
            m += 1
    if m > 1: CommonDf.at[Index, 'Max_sp'] = m
CommonDf.Max_sp.value_counts()



################# figure 3 Multi sp substitution
CommonDf['MultiSp_Subs'] = (CommonDf['Subs'] - CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs']].sum(axis=1))
CommonDf['%MultiSp_Subs'] = (CommonDf['Subs'] - CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs']].sum(axis=1))*100/CommonDf['overlap']

len(CommonDf[CommonDf['MultiSp_Subs'] > 0])#5669
plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
fig, ((ax2,)) = plt.subplots(nrows=1, ncols=1,figsize=(3,5))
# figure 3g sum of all substitutions
CommonDf[['#1#_Subs', 'Convergent_Subs', 'OneInOutId', 'OnlyInGpId', 'OnlyOutGpId',
       'NoId',]].sum().plot.bar(ax=ax2)
ax2.set(ylabel="Total Number of sites", xticklabels=["#1#_Subs", 'Two identities - inconsistent', 'One identity - inconsistent', 'Only Ingroup identical', 'Only Outgroup identical',
       'No identity',])
plt.setp(ax2.get_xmajorticklabels(), rotation=30, ha='right')

plt.tight_layout(pad=1.08)
plt.savefig('Fig4_MultiSpSubs.pdf')


################# Supplemental figure Subsampling for mean susbtitution frequencies
Sub_mean = pd.Series(dtype=np.float)
n=0
while n < 5000:
    Sub_mean.at[n] = CommonDf['%Subs'].sample(n=6580).mean()
    n += 1
HS_mean = pd.Series(dtype=np.float)
n=0
while n < 5000:
    HS_mean.at[n] = CommonDf['%HS_Subs'].sample(n=6580).mean()
    n += 1
PT_mean = pd.Series(dtype=np.float)
n=0
while n < 5000:
    PT_mean.at[n] = CommonDf['%PT_Subs'].sample(n=6580).mean()
    n += 1
GG_mean = pd.Series(dtype=np.float)
n=0
while n < 5000:
    GG_mean.at[n] = CommonDf['%GG_Subs'].sample(n=6580).mean()
    n += 1
NL_mean = pd.Series(dtype=np.float)
n=0
while n < 5000:
    NL_mean.at[n] = CommonDf['%NL_Subs'].sample(n=6580).mean()
    n += 1
P1_mean = pd.Series(dtype=np.float)
n=0
while n < 5000:
    P1_mean.at[n] = CommonDf['%#1#_Subs'].sample(n=6580).mean()
    n += 1

MeanSubDf = pd.DataFrame(columns=['mean', 'sd', 'normality_pvalue'])

plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3,figsize=(10,8))

counts1, bins1, patches1 = ax1.hist(Sub_mean, bins=100)
ax1.set(xticks=np.around(np.linspace(bins1.min(), bins1.max(), 5), decimals =3),
        xlabel="Mean % per site subs - overall",
        ylabel="Subsamples",)
plt.setp(ax1.get_xmajorticklabels(), rotation=30)
zscore, pvalue = normaltest(Sub_mean)
MeanSubDf.loc['overall', ['mean', 'sd', 'normality_pvalue']] = Sub_mean.mean(), Sub_mean.std(), pvalue
counts1, bins1, patches1 = ax2.hist(HS_mean, bins=100)
ax2.set(xticks=np.around(np.linspace(bins1.min(), bins1.max(), 5), decimals =3),
        xlabel="Mean % per site subs - human",
        ylabel="Subsamples",)
plt.setp(ax2.get_xmajorticklabels(), rotation=30)
zscore, pvalue = normaltest(HS_mean)
MeanSubDf.loc['human', ['mean', 'sd', 'normality_pvalue']] = HS_mean.mean(), HS_mean.std(), pvalue
counts1, bins1, patches1 = ax3.hist(PT_mean, bins=100)
ax3.set(xticks=np.around(np.linspace(bins1.min(), bins1.max(), 5), decimals =3),
        xlabel="Mean % per site subs - chimpanzee",
        ylabel="Subsamples",)
plt.setp(ax3.get_xmajorticklabels(), rotation=30)
zscore, pvalue = normaltest(PT_mean)
MeanSubDf.loc['chimpanzee', ['mean', 'sd', 'normality_pvalue']] = PT_mean.mean(), PT_mean.std(), pvalue
counts1, bins1, patches1 = ax4.hist(GG_mean, bins=100)
ax4.set(xticks=np.around(np.linspace(bins1.min(), bins1.max(), 5), decimals =3),
        xlabel="Mean % per site subs - gorilla",
        ylabel="Subsamples",)
plt.setp(ax4.get_xmajorticklabels(), rotation=30)
zscore, pvalue = normaltest(GG_mean)
MeanSubDf.loc['gorilla', ['mean', 'sd', 'normality_pvalue']] = GG_mean.mean(), GG_mean.std(), pvalue
counts1, bins1, patches1 = ax5.hist(NL_mean, bins=100)
ax5.set(xticks=np.around(np.linspace(bins1.min(), bins1.max(), 5), decimals =3),
        xlabel="Mean % per site subs - gibbon",
        ylabel="Subsamples",)
plt.setp(ax5.get_xmajorticklabels(), rotation=30)
zscore, pvalue = normaltest(NL_mean)
MeanSubDf.loc['gibbon', ['mean', 'sd', 'normality_pvalue']] = NL_mean.mean(), NL_mean.std(), pvalue
counts1, bins1, patches1 = ax6.hist(P1_mean, bins=100)
ax6.set(xticks=np.around(np.linspace(bins1.min(), bins1.max(), 5), decimals =3),
        xlabel="Mean % per site subs - #1#",
        ylabel="Subsamples",)
plt.setp(ax6.get_xmajorticklabels(), rotation=30)
zscore, pvalue = normaltest(P1_mean)
MeanSubDf.loc['#1#', ['mean', 'sd', 'normality_pvalue']] = P1_mean.mean(), P1_mean.std(), pvalue
plt.tight_layout(pad=1.08)
plt.savefig('Supplemental_Figure_1_MeanSubstituionFreq.pdf')

MeanSubDf.to_csv('../Tables/MeanSubDf.tsv', sep='\t', index=True, header=True)

################# Branch length normalization

NormDf = CommonDf[(CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs']] > 0).any(axis=1)].copy()

NormDf['HS_norm'] = NormDf['HS_Subs']/np.sum(NormDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs']], axis=1)
NormDf['PT_norm'] = NormDf['PT_Subs']/np.sum(NormDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs']], axis=1)
NormDf['GG_norm'] = NormDf['GG_Subs']/np.sum(NormDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs']], axis=1)
NormDf['NL_norm'] = NormDf['NL_Subs']/np.sum(NormDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs']], axis=1)
NormDf['#1#_norm'] = NormDf['#1#_Subs']/np.sum(NormDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs']], axis=1)

NormLs = [x for x in list(NormDf.columns) if x.endswith('_norm')]
SuffixLs = ['HS', 'PT', 'GG', 'NL', '#1#']

NormNoOneDf = NormDf.copy(deep=True)
for suf in SuffixLs:
    # NormDf = NormDf[~((NormDf[f'{suf}_Subs'] < 3) & (NormDf[f'{suf}_norm'] == 1))]
    NormNoOneDf = NormNoOneDf[~((NormNoOneDf[f'{suf}_norm'] == 1))]

NormNoOneDf.to_csv('../Tables/NormNoOneDf.tsv', sep='\t', header=True, index=False)

NameSer = pd.Series(data=['Human', 'Chimpanzee', 'Gorilla','Gibbon', '#1#',], index=NormLs)

################# figure 4 qqplot
from statsmodels.graphics.gofplots import qqplot
plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3,figsize=(15,10))
qqplot(CommonDf['%Subs'], line='s', ax=ax1)
ax1.set_title('% per site subs',fontweight="bold",size= 16)
ax1.set_xlabel('Theoretical Quantiles', Fontsize = 12)
ax1.set_xlabel('Sample Quantiles', Fontsize = 12)
ax1.set(ylim=(-3,35))
qqplot(NormDf.HS_norm, line='s', ax=ax2)
ax2.set_title('Normalised human subs',fontweight="bold",size= 16)
ax2.set_xlabel('Theoretical Quantiles', Fontsize = 12)
ax2.set_xlabel('Sample Quantiles', Fontsize = 12)
ax2.set(ylim=(-0.1,1.1))
qqplot(NormDf.PT_norm, line='s', ax=ax3)
ax3.set_title('Normalised chimpanzee subs',fontweight="bold",size= 16)
ax3.set_xlabel('Theoretical Quantiles', Fontsize = 12)
ax3.set_xlabel('Sample Quantiles', Fontsize = 12)
ax3.set(ylim=(-0.1,1.1))
qqplot(NormDf.GG_norm, line='s', ax=ax4)
ax4.set_title('Normalised gorilla subs',fontweight="bold",size= 16)
ax4.set_xlabel('Theoretical Quantiles', Fontsize = 12)
ax4.set_xlabel('Sample Quantiles', Fontsize = 12)
ax4.set(ylim=(-0.1,1.1))
qqplot(NormDf.NL_norm, line='s', ax=ax5)
ax5.set_title('Normalised gibbon subs',fontweight="bold",size= 16)
ax5.set_xlabel('Theoretical Quantiles', Fontsize = 12)
ax5.set_xlabel('Sample Quantiles', Fontsize = 12)
ax5.set(ylim=(-0.1,1.1))
qqplot(NormDf['#1#_norm'], line='s', ax=ax6)
ax6.set_title('Normalised #1# subs',fontweight="bold",size= 16)
ax6.set_xlabel('Theoretical Quantiles', Fontsize = 12)
ax6.set_xlabel('Sample Quantiles', Fontsize = 12)
ax6.set(ylim=(-0.1,1.1))
plt.tight_layout(pad=1.08)
plt.savefig('Fig4_qqplt.pdf')


################# figure 5 Molecular clock and Normalized branch lengths

#confidence interval
MolClockSer = pd.Series({'HS':0.27, 'PT':0.27, 'GG':0.37, 'NL':1.41, '#1#':0.09,})
MolClockNormDf = pd.DataFrame()
for Indx, value in MolClockSer.iteritems():
    print(Indx, value)
    MolClockNormDf.at[Indx, 'length'] = value/MolClockSer.sum()
MolClockNormDf['Min'] = MolClockNormDf['length'] - MolClockNormDf.at['#1#', 'length']
MolClockNormDf['Max'] = MolClockNormDf['length'] + MolClockNormDf.at['#1#', 'length']
MolClockNormDf.to_csv('./MolClockNormDf.tsv', sep='\t', header=True, index=False)

SuffixLs = ['HS', 'PT', 'GG', 'NL', '#1#']
ConformMolColDf = NormDf.copy(deep=True)
for branch in SuffixLs:#within confidence range on each branch
    ConformMolColDf = ConformMolColDf[ConformMolColDf[f'{branch}_norm'] > MolClockNormDf.at[branch, 'Min']]
    ConformMolColDf = ConformMolColDf[ConformMolColDf[f'{branch}_norm'] < MolClockNormDf.at[branch, 'Max']]
ConformMolColDf[['HS_norm', 'PT_norm', 'GG_norm', 'NL_norm', '#1#_norm']].describe()

#sub sampling
SubSamDf = pd.DataFrame(columns=['HS_norm', 'PT_norm', 'GG_norm', 'NL_norm', '#1#_norm'], dtype=np.float)
n=0
while n < 5000:
    SubSamMean = NormDf[['HS_norm', 'PT_norm', 'GG_norm', 'NL_norm', '#1#_norm']].sample(n=117).mean()
    SubSamDf.loc[n, SubSamMean.index] = SubSamMean.values
    n += 1

#subsamples within range
ConformSubSamDf = SubSamDf 
for branch in SuffixLs:
    ConformSubSamDf = ConformSubSamDf[ConformSubSamDf[f'{branch}_norm'] > MolClockNormDf.at[branch, 'Min']]
    ConformSubSamDf = ConformSubSamDf[ConformSubSamDf[f'{branch}_norm'] < MolClockNormDf.at[branch, 'Max']]

len(ConformSubSamDf) 
#comparision of conformed families with subsamles
CombDf = NormNoOneDf[NormLs].copy(deep=True)
CombDf['type'] = 'All non-unitary'
ConformMolColDf['type']  = 'Within range'
SubSamDf['type']  = 'Subsample Mean'
CombDf = CombDf.append(ConformMolColDf)
CombDf = CombDf.append(SubSamDf)
CombDf = CombDf.melt(id_vars='type', value_vars=NormLs, var_name='branch', value_name='length')
plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
plt.figure(figsize=(20,10))

sns.boxplot(x='branch', y='length', data=CombDf, hue='type')
sns.stripplot(x='branch',y='length',data=CombDf, jitter=True, dodge=True, color='.3', hue='type')

plt.xticks(np.arange(5), ['Human', 'Chimpanzee', 'Gorilla','Gibbon', '#1#',])
plt.tick_params(axis='both', which='major', labelsize=22)
plt.ylabel('Normalised substituions', size=25)
plt.xlabel('Branch', size=25)
plt.legend(fontsize='large')
plt.tight_layout(pad=1.08)
plt.savefig('Fig5_ConformedFamilies_Normalized_Distribution.pdf')



NoBranchOverDf = pd.DataFrame(columns=NormDf.columns, dtype=np.int)
OneBranchOverDf = pd.DataFrame(columns=NormDf.columns, dtype=np.int)
TwoBranchOverDf = pd.DataFrame(columns=NormDf.columns, dtype=np.int)
ThreeBranchBranchOverDf = pd.DataFrame(columns=NormDf.columns, dtype=np.int)
FourBranchBranchOverDf = pd.DataFrame(columns=NormDf.columns, dtype=np.int)

for Indx, row in NormDf.iterrows():
    n=0
    BranchSer = pd.Series()
    for branch in SuffixLs:
        if row[f'{branch}_norm'] > MolClockNormDf.at[branch, 'Max']:
            BranchSer.at[n] = branch
            n += 1
    if n == 0: NoBranchOverDf.loc[Indx] = row.values
    elif n == 1:
        OneBranchOverDf.loc[Indx, row.index] = row.values
        OneBranchOverDf.loc[Indx, 'OverSp'] = BranchSer.at[0]
    elif n == 2:
        TwoBranchOverDf.loc[Indx, row.index] = row.values
        TwoBranchOverDf.loc[Indx, 'OverSp'] = f"{BranchSer.at[0]}, {BranchSer.at[1]}"
    elif n == 3:
        ThreeBranchBranchOverDf.loc[Indx, row.index] = row.values
        ThreeBranchBranchOverDf.loc[Indx, 'OverSp'] = f"{BranchSer.at[0]}, {BranchSer.at[1]}, {BranchSer.at[2]}"
    elif n == 4:
        FourBranchBranchOverDf.loc[Indx, row.index] = row.values
        FourBranchBranchOverDf.loc[Indx, 'OverSp'] = f"{BranchSer.at[0]}, {BranchSer.at[1]}, {BranchSer.at[2]}, {BranchSer.at[3]}"
    

print(len(NoBranchOverDf))
print(len(OneBranchOverDf))
print(len(TwoBranchOverDf))
print(len(ThreeBranchBranchOverDf))
print(len(FourBranchBranchOverDf))

NormDf.loc[NoBranchOverDf.index, 'NumOver'] = 0
NormDf.loc[OneBranchOverDf.index, 'NumOver'] = 1
NormDf.loc[TwoBranchOverDf.index, 'NumOver'] = 2
NormDf.loc[ThreeBranchBranchOverDf.index, 'NumOver'] = 3
NormDf.loc[FourBranchBranchOverDf.index, 'NumOver'] = 4

NormDf.loc[NoBranchOverDf.index, 'SpOver'] = '-'
NormDf.loc[OneBranchOverDf.index, 'SpOver'] = OneBranchOverDf.OverSp
NormDf.loc[TwoBranchOverDf.index, 'SpOver'] = TwoBranchOverDf.OverSp
NormDf.loc[ThreeBranchBranchOverDf.index, 'SpOver'] = ThreeBranchBranchOverDf.OverSp
NormDf.loc[FourBranchBranchOverDf.index, 'SpOver'] = FourBranchBranchOverDf.OverSp


#outliers
HSoutliersDf = OneBranchOverDf[(OneBranchOverDf.OverSp == 'HS') & (OneBranchOverDf.HS_Subs > 3)]
HSoutliersDf = HSoutliersDf.sort_values(['%HS_Subs', 'HS_norm'], ascending=[False, False]).copy(deep=True)
HSoutliersDf[['HS', 'Description', '%HS_Subs', 'HS_norm', 'overlap', 'RelOverlap']].head(n=20)

PToutliersDf = OneBranchOverDf[(OneBranchOverDf.OverSp == 'PT') & (OneBranchOverDf.PT_Subs > 3)]
PToutliersDf = PToutliersDf.sort_values(['%PT_Subs', 'PT_norm'], ascending=[False, False]).copy(deep=True)
PToutliersDf[['HS', 'PT', 'Description', '%PT_Subs', 'PT_norm', 'overlap', 'RelOverlap', 'PT_Subs']].head(n=20)

P1outliersDf = OneBranchOverDf[(OneBranchOverDf.OverSp == '#1#') & (OneBranchOverDf['#1#_Subs'] > 3)]
P1outliersDf = P1outliersDf.sort_values(['%#1#_Subs', '#1#_norm'], ascending=[False, False]).copy(deep=True)
P1outliersDf[['HS', 'Description', '%#1#_Subs', '#1#_norm', 'overlap', 'RelOverlap']].head(n=20)

GGoutliersDf = OneBranchOverDf[(OneBranchOverDf.OverSp == 'GG') & (OneBranchOverDf.GG_Subs > 3)]
GGoutliersDf = GGoutliersDf.sort_values(['%GG_Subs', 'GG_norm'], ascending=[False, False]).copy(deep=True)
GGoutliersDf[['HS', 'GG', 'Description', '%GG_Subs', 'GG_norm', 'overlap', 'RelOverlap', 'GG_Subs']].head(n=20)

NLoutliersDf = OneBranchOverDf[(OneBranchOverDf.OverSp == 'NL') & (OneBranchOverDf.NL_Subs > 3)]
NLoutliersDf = NLoutliersDf.sort_values(['%NL_Subs', 'NL_norm'], ascending=[False, False]).copy(deep=True)
NLoutliersDf[['HS', 'NL', 'Description', '%NL_Subs', 'NL_norm', 'overlap', 'RelOverlap', 'NL_Subs']].head(n=20)

#conserved in Gibbon over the thresholod on all other branches
NormDf[(NormDf.NL_Subs == 0) & (NormDf.NumOver == 4)][['HS', 'NL', 'Description', 'HS_Subs',
                                                       'PT_Subs', 'GG_Subs', '#1#_Subs', 'Overlap', 'RelOverlap']]

#### comments
NormDf.set_index('HS', drop=False, inplace=True)
NormDf.loc['ENST00000450565','Comment'] = 'Most divergent gene on the human branch'
NormDf.loc['ENST00000369951','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000637878','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000259845','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000454136','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000391916','Comment'] = 'Many subsituions are within a short block'

NormDf.loc['ENST00000542996','Comment'] = 'Most divergent gene on the chimpanzee branch'
NormDf.loc['ENST00000296280','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
NormDf.loc['ENST00000380041','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
NormDf.loc['ENST00000371633','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
NormDf.loc['ENST00000343470','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
NormDf.loc['ENST00000316401','Comment'] = 'True human-chimpanzee ortholgs were lost'
NormDf.loc['ENST00000640237','Comment'] = 'First half of the alighment is unrelaible even after filtering'
NormDf.loc['ENST00000228468','Comment'] = 'First half of the alighment is unrelaible even after filtering'

NormDf.loc['ENST00000274520','Comment'] = 'Most divergent gene on the #1# branch'
NormDf.loc['ENST00000299191','Comment'] = 'Among the top 5 divergent genes on the #1# branch'
NormDf.loc['ENST00000433976','Comment'] = 'Among the top 5 divergent genes on the #1# branch'
NormDf.loc['ENST00000295622','Comment'] = 'Among the top 5 divergent genes on the #1# branch'
NormDf.loc['ENST00000392027','Comment'] = 'Among the top 5 divergent genes on the #1# branch'
NormDf.loc['ENST00000405093','Comment'] = 'Tandem duplication'
NormDf.loc['ENST00000359741','Comment'] = 'Many subsituions are within a short block'
NormDf.loc['ENST00000397893','Comment'] = 'Many subsituions are within a short block'


NormDf.loc['ENST00000255499','Comment'] = 'Most divergent gene on the gorilla branch'
NormDf.loc['ENST00000432056','Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
NormDf.loc['ENST00000340083','Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
NormDf.loc['ENST00000449199','Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
NormDf.loc['ENST00000331952','Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
NormDf.loc['ENST00000370177','Comment'] = 'Many subsituions are within a short block'
NormDf.loc['ENST00000375098','Comment'] = 'Many subsituions are within a short block'
NormDf.loc['ENST00000359878','Comment'] = 'Many subsituions are within a short block'

NormDf.loc['ENST00000513010','Comment'] = 'Most divergent gene on the human branch'
NormDf.loc['ENST00000345088','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000398462','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000393330','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000397301','Comment'] = 'Among the top 5 divergent genes on the human branch'
NormDf.loc['ENST00000641401','Comment'] = 'Many subsituions are within a short block'
NormDf.loc['ENST00000641401','Comment'] = 'Tandem duplication'
NormDf.loc['ENST00000527139','Comment'] = 'Many subsituions are within a short block'

#save the norm data frame
NormDf.to_csv('../Tables/NormDf.tsv', sep='\t', header=True, index=False)



