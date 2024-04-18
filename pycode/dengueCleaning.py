import pandas as pd
import os
from datetime import datetime, timedelta
import numpy as np
from psmatching.utilities import *

pd.options.display.float_format = "{:.2f}".format
data_dir = 'data'
# data_dir = f'{os.path.dirname(os.getcwd())}/data'
pwd = 'DenIV2023fluid'
warning_sign = ['Abdominal_Pain_Tenderness', 'Clinical_fluid_accumulation',
                'Mucosal_Bleed', 'Lethargy', 'Hepatomegaly', 'HCT', 'Hypotension_for_age','Persistent_vomiting']
# clean files
at_pre_raw = pd.read_csv(f'{data_dir}/Final Extracted 2005 - 2008 data - At Presentation+.csv', header=None,
                         index_col=None)
at_pre_raw = at_pre_raw.replace('-', None)
hosp_raw = pd.read_csv(f'{data_dir}/Final Extracted 2005 - 2008 data - Hospitalised Daily+.csv', header=None,
                       index_col=None)
hosp_raw = hosp_raw.replace('-', None)
fluid_raw = pd.read_csv(f'{data_dir}/Final Extracted Treatment - IV Fluid.csv', header=0)
fluid_raw = fluid_raw.replace('-', None)
fluid_raw = fluid_raw.replace('0.45 NS', '0.45NS')
fluid_raw = fluid_raw.replace('nonrmal Saline', 'Normal Saline')
fluid_raw = fluid_raw.replace('half', 'Half')
fluid_raw = fluid_raw.replace('H2O', 'Water (H2O)')
fluid_raw = fluid_raw.replace('dextrose', 'Dextrose')
ftype = []

for ind, row in fluid_raw.iterrows():
    if row['Fluid Name'] in ['Gelafundin']:
        ftype.append('Colloid')
    elif row['Fluid Name'].lower() in ['Dextrose s'.lower(), 'Dextrose Saline'.lower(), 'Dextrose se'.lower(),
                                       'Hartmann'.lower(), "Hartmann's solution".lower(), "Hartmen's solution".lower(),
                                       'Normal Saline'.lower(), 'Ringer Lactate'.lower()]:
        ftype.append('Crystalloid1')
    elif row['Fluid Name'] in ['Dextrose 5%']:
        ftype.append('Crystalloid3')
    elif row['Fluid Name'].lower() in ['0.45NS'.lower(), 'Half strength D/S'.lower(),
                                       'Half strength Dextrose Saline'.lower(),
                                       'Half strength Normal Saline'.lower()]:
        ftype.append('Crystalloid2')
    elif row['Fluid Name'].lower() in ['Water (H2O)'.lower(), 'Mist KCl'.lower()]:
        ftype.append('Oral')
    else:
        ftype.append('error')
fluid_raw['Fluid_Type'] = ftype
fluid_raw.to_csv(f'{data_dir}/fluids_raw.csv')
blood_raw = pd.read_csv(f'{data_dir}/Final Extracted Treatment - Blood.csv', header=0)
blood_raw = blood_raw.replace('-', None)
dateformat = '%d/%m/%Y'
at_pre = at_pre_raw.iloc[2:, ]
at_pre.columns = at_pre_raw.iloc[1]
at_pre.set_index('Study Number', inplace=True)
condition = (at_pre['DHF 1997'] == 'Yes') | (at_pre['DSS 1997'] == 'Yes') | (at_pre['Severe Dengue 2009'] == 'Yes')
at_pre = at_pre[~condition]
HCT_at_pre = []
for ind, row in at_pre.iterrows():
    if (row['HCT Change >=20% and PLT <50'] is None) & (row['HCT Change >=20%'] is None):
        HCT_at_pre.append(np.nan)
    elif (row['HCT Change >=20% and PLT <50'] == 'Yes') | (row['HCT Change >=20%'] == 'Yes'):
        HCT_at_pre.append('Yes')
    else:
        HCT_at_pre.append('No')
at_pre['HCT'] = HCT_at_pre
at_pre_warning_sign = ['Abdominal Pain/Tenderness', 'Clinical fluid accumulation',
                       'Mucosal Bleed', 'Lethargy', 'Hepatomegaly', 'HCT', 'Hypotension for age','Persistent Vomiting']
for i in range(len(warning_sign)):
    at_pre.rename(columns={at_pre_warning_sign[i]: warning_sign[i]}, inplace=True)

studyno_filtered = list(at_pre.index)
hosp = hosp_raw.iloc[2:, ]
hosp.columns = hosp_raw.iloc[1,]
hosp = hosp[hosp['Study_Number'].isin(studyno_filtered)]
intday = hosp['Day'].astype(int)
hosp['Day'] = intday
hosp.set_index(['Study_Number', 'Day'], inplace=True)
hosp['Admission Date'] = pd.to_datetime(hosp['Admission Date'], format=dateformat)
blood_raw['Start_Date'] = pd.to_datetime(blood_raw['Start_Date'], format=dateformat)
fluid_raw['Start_Date'] = pd.to_datetime(fluid_raw['Start_Date'], format=dateformat)
fluid_raw['End_Date'] = pd.to_datetime(fluid_raw['End_Date'], format=dateformat)
date = []
for ind, row in hosp.iterrows():
    tempdate = row['Admission Date'] + timedelta(days=ind[1] - 1)
    date.append(tempdate)
hosp.insert(2, 'Date', date)
fluid = fluid_raw[fluid_raw['Study_Number'].isin(studyno_filtered)]
fluid.to_csv(f'{data_dir}/fluids.csv')
treatment = []
studyno = studyno_filtered.copy()
exclusion1_count = 0  # criteria 1: severe dengue at presentation (remove)
for sno in studyno_filtered:
    df = hosp.xs(sno, level='Study_Number')
    condition = (df['DHF 1997'] == 'Yes') | (df['DSS 1997'] == 'Yes') | (
            df['Severe Dengue 2009'] == 'Yes')
    given_df = df[(df['Given'] == 'Yes')]
    dengue_df = df[condition]
    if len(dengue_df) > 0:
        if df[condition].index[0] == 1:
            print(sno)
            hosp.drop(((sno, ind) for ind in df.index), inplace=True)
            studyno.remove(sno)
            exclusion1_count += 1

print('exclusion criteria1: ', exclusion1_count)

if os.path.exists(f'{data_dir}/daily_data.csv'):
    cleaned_df = pd.read_csv(f'{data_dir}/daily_data.csv',header=0,index_col=[0,1])
else:
    # cleaned_df is the cleaned data of hospitalization data
    cleaned_df = pd.DataFrame(index=hosp.index)
    gender = []
    age = []
    CCMI = []
    dengue = []
    blood = []
    HCT = []
    heartfailure = []
    given_fluid = []
    for ind, row in hosp.iterrows():  # ind is (study_number, day)
        # outcome: DHF 1997/DSS 1997/Severe Dengue == 1
        if (row['DHF 1997'] is None) & (row['DSS 1997'] is None) & (
                row['Severe Dengue 2009'] is None):
            dengue.append(np.nan)
        elif (row['DHF 1997'] == 'Yes') | (row['DSS 1997'] == 'Yes') | (
                row['Severe Dengue 2009'] == 'Yes'):
            dengue.append('Yes')
        else:
            dengue.append('No')
        # warning sign - HCT
        if (row['HCT Change >=20% and PLT <50'] is None) & (row['HCT Change >=20%'] is None):
            HCT.append(np.nan)
        elif (row['HCT Change >=20% and PLT <50'] == 'Yes') | (row['HCT Change >=20%'] == 'Yes'):
            HCT.append('Yes')
        else:
            HCT.append('No')
        sn = ind[0]  # study number
        age.append(at_pre.loc[sn, 'Age'])  # age
        gender.append(at_pre.loc[sn, 'Gender'])  # gender
        CCMI.append(at_pre.loc[sn, 'Charlson Score'])  # CCMI
        heartfailure.append(at_pre.loc[sn,'Heart_Failure']) # heart failure
        # blood_product is given or not
        # if blood_product is given, total volume on the specific day is recorded
        if sn in list(blood_raw['Study_Number']):
            temp = blood_raw[blood_raw['Study_Number'] == sn]
            if row['Date'] in list(temp['Start_Date']):
                temp_date = row['Date']
                blood_date = temp[temp['Start_Date'] == temp_date]
                volume = sum(blood_date['Unit/Volume'])
                blood.append(volume)
            else:
                blood.append(0)
        else:
            blood.append(0)
        if sn in list(fluid['Study_Number']):
            temp = fluid[fluid['Study_Number'] == sn]
            if row['Date'] in list(temp['Start_Date']):
                temp_date = row['Date']
                fluid_date = temp[temp['Start_Date'] == temp_date]
                given_fluid.append('Yes')
            else:
                given_fluid.append('No')
        else:
            given_fluid.append('No')

    # add columns to cleaned_df
    hosp['HCT'] = HCT
    cleaned_df['Date'] = hosp['Date']
    cleaned_df['Severe Dengue'] = dengue
    cleaned_df['given_fluid'] = given_fluid
    cleaned_df['Age'] = age
    cleaned_df['Gender'] = gender
    cleaned_df['CCMI'] = CCMI
    cleaned_df['Heart_Failure'] = heartfailure
    cleaned_df['blood_treatment'] = blood

    # add warning signs before or after treatment
    cleaned_df[
        ['Abdominal_Pain_Tenderness', 'Clinical_fluid_accumulation', 'Mucosal_Bleed', 'Lethargy', 'Hepatomegaly',
         'HCT',
         'Hypotension_for_age', 'Narrow_Pulse_Pressure','Persistent_vomiting']] = hosp[
        ['Abdominal Pain/Tenderness', 'Clinical fluid accumulation (Daily)', 'Mucosal Bleed', 'Lethargy', 'Hepatomegaly',
         'HCT',
         'Hypotension for age (Daily)', 'Narrow Pulse Pressure (Daily)','Persistent Vomiting']]
    cleaned_df.replace('No', 0, inplace=True)
    cleaned_df.replace('Yes', 1, inplace=True)
    cleaned_df.to_csv(f'{data_dir}/daily_data.csv')  # write to csv file

# df is the data for each patient over the stay period
df = pd.DataFrame(index=studyno,
                  columns=['outcome', 'stay','days_to_SD',
                           'age', 'gender', 'CCMI','Heart_Failure',
                           'Abdominal_Pain_Tenderness_bf', 'Abdominal_Pain_Tenderness_af',
                           'Clinical_fluid_accumulation_bf',
                           'Clinical_fluid_accumulation_af',
                           'Mucosal_Bleed_bf', 'Mucosal_Bleed_af',
                           'Lethargy_bf', 'Lethargy_af', 'Hepatomegaly_bf', 'Hepatomegaly_af',
                           'HCT_bf', 'HCT_af', 'Hypotension_for_age_bf', 'Hypotension_for_age_af', 'Persistent_vomiting_bf',
                           'Persistent_vomiting_af',
                           ])

for sno in studyno:
    temp_df = cleaned_df.xs(sno, level='Study_Number')  # temporary df for the specific patient
    df.loc[sno,'stay'] = len(temp_df)
    df.loc[sno, 'age'] = list(temp_df['Age'])[0]
    if list(temp_df['Gender'])[0] == 'Male':
        df.loc[sno, 'gender'] = 1
    elif list(temp_df['Gender'])[0] == 'Female':
        df.loc[sno, 'gender'] = 0
    df.loc[sno, 'CCMI'] = list(temp_df['CCMI'])[0]
    if list(temp_df['Heart_Failure'])[0] == 1:
        df.loc[sno, 'Heart_Failure'] = 1
    elif list(temp_df['Heart_Failure'])[0] == 0:
        df.loc[sno, 'Heart_Failure'] = 0
    condition = (temp_df['Severe Dengue'] == 1)  # The patient has developed SD
    condition2 = (temp_df['given_fluid'] == 1)  # The patient has been given fluid
    # outcome
    if temp_df['Severe Dengue'].sum()>0:
        df.loc[sno, 'outcome'] = 1
    else:
        df.loc[sno, 'outcome'] = 0

    if temp_df['given_fluid'].sum()>0:
        df.loc[sno, 'treatment'] = 1  # given fluid is true
        # if fluid is given on the first day, any warning sign at presentation is recorded
        if temp_df[condition2].index[0] == 1:
            for col in warning_sign:
                if at_pre.loc[sno, col] == 1:
                    df.loc[sno, f'{col}_bf'] = 1
                    df.loc[sno, f'{col}_af'] = 1
                else:
                    df.loc[sno, f'{col}_bf'] = 0
                    if temp_df[col].sum()>0:
                        df.loc[sno, f'{col}_af'] = 1
                    else:
                        df.loc[sno, f'{col}_af'] = 0
        else:  # if fluid is given after admission day
            for col in warning_sign:
                if not (temp_df[col].sum()>0):  # no warning sign over the period
                    df.loc[sno, f'{col}_bf'] = 0
                    df.loc[sno, f'{col}_af'] = 0
                # warning sign before treatment
                elif temp_df[(temp_df[col] == 1)].index[0] < temp_df[condition2].index[0]:
                    df.loc[sno, f'{col}_bf'] = 1
                    if temp_df[(temp_df[col] == 1)].index[-1] < temp_df[condition2].index[0]:
                        df.loc[sno, f'{col}_af'] = 0
                    else:
                        df.loc[sno, f'{col}_af'] = 1
                elif temp_df[(temp_df[col] == 1)].index[0] >= temp_df[condition2].index[0]:
                    df.loc[sno, f'{col}_bf'] = 0
                    df.loc[sno, f'{col}_af'] = 1

    else:  # no fluid treatment
        df.loc[sno, 'treatment'] = 0
        for col in warning_sign:
            df.loc[sno, f'{col}_af'] = 0
            if not (temp_df[col].sum()>0):
                df.loc[sno, f'{col}_bf'] = 0
            else:
                df.loc[sno, f'{col}_bf'] = 1

for col in ['age', 'CCMI']:
    df[col] = pd.to_numeric(df[col])

# find the mean days of treatment before severe dengue for SD group
sd_fdays = {}
for sno in studyno:
    temp_df = cleaned_df.xs(sno, level='Study_Number')
    condition = (temp_df['Severe Dengue'] == 1)
    if temp_df['Severe Dengue'].sum()>0:  # if the patient has developed SD
        date_sd = list(temp_df[condition]['Date'])[0]
        if sno in list(fluid['Study_Number']):  # if the patient has been given fluid
            temp_fluid = fluid[fluid['Study_Number'] == sno]
            fluid_bSD = temp_fluid[temp_fluid['Start_Date'] < date_sd]
            if len(fluid_bSD) > 0:
                fluid_bSD_dates = set(list(fluid_bSD['Start_Date']))
                sd_fdays[sno] = len(fluid_bSD_dates)
mean_treatment_days = np.mean(list(sd_fdays.values()))

# find the mean days of treatment before severe dengue for SD group
sd_bpdays = {}
sd_bFbpdays = {}
sd_aFbpdays = {}
for sno in studyno:
    temp_df = cleaned_df.xs(sno, level='Study_Number')
    condition = (temp_df['Severe Dengue'] == 1)
    if temp_df['Severe Dengue'].sum()>0:  # if the patient has developed SD
        date_sd = list(temp_df[condition]['Date'])[0]
        if sno in list(blood_raw['Study_Number']):  # if the patient has been given fluid
            temp_bp = blood_raw[blood_raw['Study_Number'] == sno]
            condition2 = (temp_df['given_fluid'] == 1)
            bp_bSD = temp_bp[temp_bp['Start_Date'] < date_sd]
            if len(bp_bSD) > 0:
                bp_bSD_dates = set(list(bp_bSD['Start_Date']))
                sd_bpdays[sno] = len(bp_bSD_dates)
            if temp_df['given_fluid'].sum()>0:
                date_fluid = list(temp_df[condition2]['Date'])[0]
                bF_bp = bp_bSD[bp_bSD['Start_Date'] < date_fluid]
                aF_bp = bp_bSD[bp_bSD['Start_Date'] >= date_fluid]
                if len(bF_bp) > 0:
                    bf_dates = set(list(bF_bp['Start_Date']))
                    sd_bFbpdays[sno] = len(bf_dates)
                if len(aF_bp) > 0:
                    af_dates = set(list(aF_bp['Start_Date']))
                    sd_aFbpdays[sno] = len(af_dates)
mean_bp_days = np.mean(list(sd_bpdays.values()))
mean_bFbp_days = np.mean(list(sd_bFbpdays.values()))
mean_aFbp_days = np.mean(list(sd_aFbpdays.values()))


for ind, row in df.iterrows():
    if row['outcome'] == 1:
        temp_df = cleaned_df.loc[ind, :]
        days_to_SD = temp_df[temp_df['Severe Dengue'] == 1].index[0] - 1
        df.loc[ind, 'days_to_SD'] = days_to_SD

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

df.to_csv(f'{data_dir}/dengue_fluid.csv')


