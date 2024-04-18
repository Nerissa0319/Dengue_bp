import pandas as pd
import os
from datetime import datetime, timedelta
import numpy as np
from psmatching.utilities import *

pd.options.display.float_format = "{:.2f}".format
data_dir = 'data'

raw_lab = pd.read_csv(f'{data_dir}/Laboratory - Normal - Hourly.csv', header=0)
labdata = raw_lab.drop('Time', axis=1)
daily_lab = labdata.groupby(['Study_Number', 'Date']).apply(lambda x: x.apply(np.nanmean))
daily_lab.to_csv(f'{data_dir}/daily_lab.csv')

raw_signs = pd.read_csv(f'{data_dir}/Signs - Hourly.csv', header=0)
raw_signs.columns = ['Study_Number', 'Date', 'Time', 'Temperature', 'SBP', 'DBP', 'Pulse_rate', 'Respiratory_rate',
                     'Oxygen_saturation', 'Others', 'Narrow_Pulse_Pressure']
raw_signs['baselineSBP'] = None
raw_signs['SBP90'] = raw_signs['SBP'] < 90  # SBP < 90

SBP_decrease = []  # SBP decrease from the first measurement(baseline) > 40
baselineSBP_ls = []
for sn in raw_signs['Study_Number'].unique():
    temp_df = raw_signs[raw_signs['Study_Number'] == sn]
    temp_len = len(temp_df)
    baseline_ind = list(pd.isna(list(temp_df['SBP']))).index(False)
    baselinesbp = list(temp_df['SBP'])[baseline_ind]
    baselineSBP_ls.extend([baselinesbp] * temp_len)
    temp_ls = (baselinesbp - temp_df['SBP']) > 40
    SBP_decrease.extend(temp_ls)
raw_signs['baselineSBP'] = baselineSBP_ls
raw_signs['SBP_decrease40'] = SBP_decrease

# if SBP is NA, SBP90, SBPdecrease are NA
for ind, row in raw_signs.iterrows():
    if pd.isna(row['SBP']):
        raw_signs.loc[ind, 'SBP90'] = None
        raw_signs.loc[ind, 'SBP_decrease40'] = None
raw_signs.replace(True, int(1), inplace=True)
raw_signs.replace(False, int(0), inplace=True)
raw_signs.replace('No', int(0), inplace=True)
raw_signs.replace('Yes', int(1), inplace=True)

raw_signs.to_csv(f'{data_dir}/bp.csv')
bp = raw_signs.loc[:, ['Study_Number', 'Date', 'Temperature', 'SBP', 'DBP', 'Pulse_rate', 'Respiratory_rate',
                       'Oxygen_saturation', 'baselineSBP']]
daily_signs = bp.groupby(['Study_Number', 'Date']).apply(lambda x: x.apply(np.nanmean))
daily_signs['SBP90'] = None
daily_signs['SBP_decrease40'] = None
daily_signs['meanSBP90'] = None
daily_signs['meanSBP_decrease40'] = None
daily_signs['Narrow_Pulse_Pressure_bp'] = None
for ind, row in daily_signs.iterrows():
    sn = ind[0]
    dte = ind[1]
    base = list(raw_signs[raw_signs['Study_Number'] == sn]['baselineSBP'])[0]
    temp_df = raw_signs[(raw_signs['Study_Number'] == sn) & (raw_signs['Date'] == dte)]
    if (temp_df['SBP90'] == 1).any():
        daily_signs.loc[ind, 'SBP90'] = 1
    else:
        daily_signs.loc[ind, 'SBP90'] = 0
    if (temp_df['SBP_decrease40'] == 1).any():
        daily_signs.loc[ind, 'SBP_decrease40'] = 1
    else:
        daily_signs.loc[ind, 'SBP_decrease40'] = 0
    if (temp_df['Narrow_Pulse_Pressure'] == 1).any():
        daily_signs.loc[ind, 'Narrow_Pulse_Pressure_bp'] = 1
    else:
        daily_signs.loc[ind, 'Narrow_Pulse_Pressure_bp'] = 0
    if row['SBP'] < 90:
        daily_signs.loc[ind, 'meanSBP90'] = 1
    else:
        daily_signs.loc[ind, 'meanSBP90'] = 0
    if (base - row['SBP']) > 40:
        daily_signs.loc[ind, 'meanSBP_decrease40'] = 1
    else:
        daily_signs.loc[ind, 'meanSBP_decrease40'] = 0
daily_signs.to_csv(f'{data_dir}/daily_signs.csv')

raw_postural = pd.read_csv(f'{data_dir}/Signs - Postural BP.csv', header=0)
raw_postural.columns = ['Study_Number', 'Date', 'Time', 'SBP_lying', 'DBP_lying', 'SBP_standing', 'DBP_standing']
raw_postural['SBP_postural_decrease'] = raw_postural['SBP_lying'] - raw_postural['SBP_standing']
raw_postural['DBP_postural_decrease'] = raw_postural['DBP_lying'] - raw_postural['DBP_standing']
raw_postural['postural_decrease'] = [None] * len(raw_postural)
for ind, row in raw_postural.iterrows():
    if (pd.isna(row['SBP_postural_decrease'])) and (pd.isna(row['DBP_postural_decrease'])):
        raw_postural.loc[ind, 'postural_decrease'] = None
    elif (row['SBP_postural_decrease'] > 20) or (row['DBP_postural_decrease'] > 10):
        raw_postural.loc[ind, 'postural_decrease'] = 1
    else:
        raw_postural.loc[ind, 'postural_decrease'] = 0

raw_postural.to_csv(f'{data_dir}/postural.csv')

posturaldata = raw_postural.loc[:, ['Study_Number', 'Date', 'SBP_lying', 'SBP_standing', 'DBP_lying', 'DBP_standing']]
daily_postural = posturaldata.groupby(['Study_Number', 'Date']).apply(lambda x: x.apply(np.nanmean))
daily_postural['postural_decrease'] = None
for ind, row in daily_postural.iterrows():
    sn = ind[0]
    dte = ind[1]
    temp_df = raw_postural[(raw_postural['Study_Number'] == sn) & (raw_signs['Date'] == dte)]
    if (temp_df['postural_decrease'] == 1).any():
        daily_postural.loc[ind, 'postural_decrease'] = 1
    else:
        daily_postural.loc[ind, 'postural_decrease'] = 0

daily_postural.to_csv(f'{data_dir}/daily_postural.csv')

merged_bp = pd.merge(daily_signs, daily_lab, on=['Study_Number', 'Date'], how='outer')
merged_bp = pd.merge(merged_bp, daily_postural, on=['Study_Number', 'Date'], how='outer')
merged_bp.reset_index(inplace=True)
hosp_data = pd.read_csv(f'{data_dir}/daily_data.csv', header=0)
dateformat = '%d/%m/%Y'
hosp_data['Date'] = pd.to_datetime(hosp_data['Date'])
merged_bp['Date'] = pd.to_datetime(merged_bp['Date'], format=dateformat)

outermerged_bp = pd.merge(hosp_data, merged_bp, on=['Study_Number', 'Date'], how='outer')
leftmerged_bp = pd.merge(hosp_data, merged_bp, on=['Study_Number', 'Date'], how='left')


def replacebp(bp):
    bp.replace('Female', 0, inplace=True)
    bp.replace('Male', 1, inplace=True)
    bp.drop(['Narrow_Pulse_Pressure'], axis=1, inplace=True)
    return bp


outermerged_bp = replacebp(outermerged_bp)
leftmerged_bp = replacebp(leftmerged_bp)


def ccmi_convert(df):
    ccmi_class = []
    for ind, row in df.iterrows():
        ccmi = row['CCMI']
        if ccmi == 0:
            ccmi_class.append('0')
        elif ccmi <= 2:
            ccmi_class.append('[1,2]')
        elif ccmi <= 4:
            ccmi_class.append('[3,4]')
        elif ccmi >= 5:
            ccmi_class.append('[5,)')
        else:
            ccmi_class.append(None)
    df['CCMI'] = ccmi_class
    return df


outermerged_bp = ccmi_convert(outermerged_bp)
leftmerged_bp = ccmi_convert(leftmerged_bp)

outermerged_bp = pd.read_csv(f'{data_dir}/outer_merged_bp.csv',header=0,index_col=0)
leftmerged_bp = pd.read_csv(f'{data_dir}/merged_bp.csv',header=0,index_col=0)
def ws_bf_day(bp):
    warning_signs = ['Abdominal_Pain_Tenderness', 'Clinical_fluid_accumulation',
                     'Mucosal_Bleed', 'Lethargy', 'Hepatomegaly', 'HCT', 'Hypotension_for_age',
                     'Persistent_vomiting']
    for ws in warning_signs:
        bp[f'{ws}_bf'] = None
    bp['ws_bf'] = None
    warning_signs_bf = []
    for ws in warning_signs:
        warning_signs_bf.append(f'{ws}_bf')
    warning_signs_bf.append('ws_bf')

    for sn in list(bp['Study_Number'].unique()):
        temp_df = bp[bp['Study_Number'] == sn]
        for d in temp_df['Day']:
            if d is not None:
                if d == 1:
                    for wsbf in warning_signs_bf:

                        bp.loc[(bp['Study_Number'] == sn) & (bp['Day'] == d), wsbf] = 0
                elif d>1:

                    temp_df2 = temp_df[temp_df['Day'] <= d]
                    for ws in warning_signs:
                        if temp_df2[ws].sum() > 0:
                            bp.loc[(bp['Study_Number'] == sn) & (bp['Day'] == d), f'{ws}_bf'] = 1
                        elif temp_df2[ws].sum()==0:
                            bp.loc[(bp['Study_Number'] == sn) & (bp['Day'] == d), f'{ws}_bf'] = 0
                temp_df = bp[bp['Study_Number'] == sn]
                if d >= 1:
                    wscheck = temp_df[warning_signs_bf].iloc[int(d - 1)]
                    if wscheck.sum() > 0:
                        bp.loc[(bp['Study_Number'] == sn) & (bp['Day'] == d), 'ws_bf'] = 1
                    elif wscheck.sum()==0:
                        bp.loc[(bp['Study_Number'] == sn) & (bp['Day'] == d), 'ws_bf'] = 0

    return bp

leftmerged_bp = ws_bf_day(leftmerged_bp)
leftmerged_bp.to_csv(f'{data_dir}/merged_bp.csv')
outermerged_bp = ws_bf_day(outermerged_bp)
outermerged_bp.to_csv(f'{data_dir}/outer_merged_bp.csv')

dengue_df = pd.read_csv(f'{data_dir}/dengue_fluid.csv',header=0,index_col=0)
bp = pd.read_csv(f'{data_dir}/merged_bp.csv',header=0)
dengue_df['SBP90'] = None
dengue_df['SBP_decrease40'] = None
dengue_df['Narrow_Pulse_Pressure_bp'] = None
dengue_df['shock_ind'] = None
dengue_df['diastolic_shock_ind'] = None
dengue_df['modified_shock_ind'] = None


for sn in dengue_df.index:
    temp_df = bp[bp['Study_Number']==sn]
    if dengue_df.loc[sn,'outcome'] == 0:
        if temp_df['SBP90'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'SBP90'] = None
        elif temp_df['SBP90'].sum()>0:
            dengue_df.loc[sn,'SBP90'] = 1
        else:
            dengue_df.loc[sn,'SBP90'] = 0
        if temp_df['SBP_decrease40'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'SBP_decrease40'] = None
        elif temp_df['SBP_decrease40'].sum()>0:
            dengue_df.loc[sn,'SBP_decrease40'] = 1
        else:
            dengue_df.loc[sn,'SBP_decrease40'] = 0
        if temp_df['Narrow_Pulse_Pressure_bp'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'Narrow_Pulse_Pressure_bp'] = None
        elif temp_df['Narrow_Pulse_Pressure_bp'].sum()>0:
            dengue_df.loc[sn,'Narrow_Pulse_Pressure_bp'] = 1
        else:
            dengue_df.loc[sn,'Narrow_Pulse_Pressure_bp'] = 0
        if temp_df['shock_ind'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'shock_ind'] = None
        else:
            dengue_df.loc[sn,'shock_ind'] = temp_df.shock_ind.mean()
        if temp_df['modified_shock_ind'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'modified_shock_ind'] = None
        else:
            dengue_df.loc[sn,'modified_shock_ind'] = temp_df.modified_shock_ind.mean()
        if temp_df['diastolic_shock_ind'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'diastolic_shock_ind'] = None
        else:
            dengue_df.loc[sn,'diastolic_shock_ind'] = temp_df.diastolic_shock_ind.mean()
    else:
        sd_day = dengue_df.loc[sn,'days_to_SD'] + 1
        temp_df = temp_df[temp_df['Day'] < sd_day]
        if temp_df['SBP90'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'SBP90'] = None
        elif temp_df['SBP90'].sum()>0:
            dengue_df.loc[sn,'SBP90'] = 1
        else:
            dengue_df.loc[sn,'SBP90'] = 0
        if temp_df['SBP_decrease40'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'SBP_decrease40'] = None
        elif temp_df['SBP_decrease40'].sum()>0:
            dengue_df.loc[sn,'SBP_decrease40'] = 1
        else:
            dengue_df.loc[sn,'SBP_decrease40'] = 0
        if temp_df['Narrow_Pulse_Pressure_bp'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'Narrow_Pulse_Pressure_bp'] = None
        elif temp_df['Narrow_Pulse_Pressure_bp'].sum()>0:
            dengue_df.loc[sn,'Narrow_Pulse_Pressure_bp'] = 1
        else:
            dengue_df.loc[sn,'Narrow_Pulse_Pressure_bp'] = 0
        if temp_df['shock_ind'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'shock_ind'] = None
        else:
            dengue_df.loc[sn,'shock_ind'] = temp_df.shock_ind.mean()
        if temp_df['modified_shock_ind'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'modified_shock_ind'] = None
        else:
            dengue_df.loc[sn,'modified_shock_ind'] = temp_df.modified_shock_ind.mean()
        if temp_df['diastolic_shock_ind'].isna().sum() == len(temp_df):
            dengue_df.loc[sn,'diastolic_shock_ind'] = None
        else:
            dengue_df.loc[sn,'diastolic_shock_ind'] = temp_df.diastolic_shock_ind.mean()



dengue_df = dengue_df[['outcome','stay','days_to_SD',
                       'age','gender','CCMI','SBP90','SBP_decrease40','Narrow_Pulse_Pressure_bp','shock_ind','diastolic_shock_ind','modified_shock_ind']]
dengue_sbp90 = dengue_df.copy()
dengue_sbp40 = dengue_df.copy()
dengue_np = dengue_df.copy()
dengue_shock = dengue_df.copy()
dengue_diastolic = dengue_df.copy()
dengue_modified = dengue_df.copy()
dfs = [dengue_sbp90,dengue_sbp40,dengue_np,dengue_modified,dengue_diastolic,dengue_shock]
for df in dfs:
    for cat_col in ['Abdominal_Pain_Tenderness', 'Clinical_fluid_accumulation', 'Mucosal_Bleed', 'Lethargy',
                'Hepatomegaly', 'HCT', 'Hypotension_for_age', 'Persistent_vomiting','Temperature','SBP','DBP','Pulse_rate',
                    'Respiratory_rate','Oxygen_saturation','white.cell.count','Neutrophils.polymorphs','lymphocytes',
                    'monocytes','basophils','eosinophils','atypical.reactive.lymphocytes','Haematocrit',
                    'Haemoglobin','Platelet']:
        df[cat_col] = None
for sn in dengue_df.index:
    temp_df = bp[bp['Study_Number'] == sn]
    condition = ['SBP90','SBP_decrease40','Narrow_Pulse_Pressure_bp']
    for i in range(3):
        temp_df = bp[(bp['Study_Number']==sn)]
        if dengue_df.loc[sn,condition[i]]==1:
            exposure_day = list(temp_df[temp_df[condition[i]] == 1]['Day'])[0]
            temp_df = temp_df[temp_df['Day']<=exposure_day]

        for cat_col in ['Abdominal_Pain_Tenderness','Clinical_fluid_accumulation','Mucosal_Bleed','Lethargy',
                        'Hepatomegaly','HCT','Hypotension_for_age','Persistent_vomiting']:
            dfs[i].loc[sn,cat_col] = temp_df[cat_col].sum()
            if dfs[i].loc[sn,cat_col]>1:
                dfs[i].loc[sn, cat_col]=1
        for con_col in ['Temperature','SBP','DBP','Pulse_rate',
                    'Respiratory_rate','Oxygen_saturation','white.cell.count','Neutrophils.polymorphs','lymphocytes',
                    'monocytes','basophils','eosinophils','atypical.reactive.lymphocytes','Haematocrit',
                    'Haemoglobin','Platelet']:
            dfs[i].loc[sn, con_col] = temp_df[con_col].mean()
    for df in [dengue_modified,dengue_diastolic,dengue_shock]:
        temp_df = bp[(bp['Study_Number'] == sn)]
        for cat_col in ['Abdominal_Pain_Tenderness', 'Clinical_fluid_accumulation', 'Mucosal_Bleed', 'Lethargy',
                        'Hepatomegaly', 'HCT', 'Hypotension_for_age', 'Persistent_vomiting']:
            df.loc[sn, cat_col] = temp_df[cat_col].sum()
            if df.loc[sn, cat_col] > 1:
                df.loc[sn, cat_col] = 1
        for con_col in ['Temperature', 'SBP', 'DBP', 'Pulse_rate',
                        'Respiratory_rate', 'Oxygen_saturation', 'white.cell.count', 'Neutrophils.polymorphs',
                        'lymphocytes',
                        'monocytes', 'basophils', 'eosinophils', 'atypical.reactive.lymphocytes', 'Haematocrit',
                        'Haemoglobin', 'Platelet']:
            df.loc[sn, con_col] = temp_df[con_col].mean()

dengue_df_summary = dengue_df.copy()
for cat_col in ['Abdominal_Pain_Tenderness', 'Clinical_fluid_accumulation', 'Mucosal_Bleed', 'Lethargy',
                'Hepatomegaly', 'HCT', 'Hypotension_for_age', 'Persistent_vomiting', 'Temperature', 'SBP', 'DBP',
                'Pulse_rate',
                'Respiratory_rate', 'Oxygen_saturation', 'white.cell.count', 'Neutrophils.polymorphs', 'lymphocytes',
                'monocytes', 'basophils', 'eosinophils', 'atypical.reactive.lymphocytes', 'Haematocrit',
                'Haemoglobin', 'Platelet']:
    dengue_df_summary[cat_col] = None
for sn in dengue_df_summary.index:
    temp_df = bp[bp['Study_Number'] == sn]
    if dengue_df_summary.loc[sn,'outcome']==1:
        sd_day = list(temp_df[temp_df['Severe.Dengue'] == 1]['Day'])[0]
        temp_df = temp_df[temp_df['Day']<sd_day]

    for cat_col in ['Abdominal_Pain_Tenderness','Clinical_fluid_accumulation','Mucosal_Bleed','Lethargy',
                    'Hepatomegaly','HCT','Hypotension_for_age','Persistent_vomiting']:
        dengue_df_summary.loc[sn,cat_col] = temp_df[cat_col].sum()
        if dengue_df_summary.loc[sn,cat_col]>1:
            dengue_df_summary.loc[sn, cat_col]=1
    for con_col in ['Temperature','SBP','DBP','Pulse_rate',
                'Respiratory_rate','Oxygen_saturation','white.cell.count','Neutrophils.polymorphs','lymphocytes',
                'monocytes','basophils','eosinophils','atypical.reactive.lymphocytes','Haematocrit',
                'Haemoglobin','Platelet']:
        dengue_df_summary.loc[sn, con_col] = temp_df[con_col].mean()

dengue_sbp40.to_csv(f'{data_dir}/dengue_sbp40.csv')
dengue_sbp90.to_csv(f'{data_dir}/dengue_sbp90.csv')
dengue_np.to_csv(f'{data_dir}/dengue_narrow.csv')
dengue_shock.to_csv(f'{data_dir}/dengue_shock.csv')
dengue_diastolic.to_csv(f'{data_dir}/dengue_diastolic.csv')
dengue_modified.to_csv(f'{data_dir}/dengue_modified.csv')
dengue_df_summary.to_csv(f'{data_dir}/dengue_sd.csv')

