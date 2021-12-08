# -*- coding: utf-8 -*-
"""
Created on Wed May 27 08:40:28 2020

@author: manz184215
"""
import numpy as np
import pandas as pd

def convert_data_phenotips_raw(file, hpo_cols=True):
    """
    Extract HPO terms from a raw phenotips JSON file. Returns a dataframe.
    
    Parameters
    ----------
    file : str
        Path to the JSON file to be processed.
    hpo_cols: bool
        Whether to hot encode the HPO terms and create a seperate column for each one
        
    Returns
    -------
        Pandas DataFrame containing the extracted HPO terms
    """
    import json
    import pandas as pd
    import numpy as np
        
    with open(file, "r") as f:
        data = f.read()
           
    data = data.replace("'s", "s")
    data = data.replace("'marbled'", "marbled")
    
    with open(file, 'w') as myfile:
      myfile.write(data)
        
    with open(file) as train_file:
        dict_train = json.load(train_file)
    # converting json dataset from dictionary to dataframe
    df_phenotips = pd.DataFrame(dict_train)
    
    df_phenotips = df_phenotips.reset_index(drop=True)    
     
    Count = 0
    
    df_phenotips['hpo_all'], df_phenotips['hpo_all_name'] = '' ,'' 
    
    for value_dict in df_phenotips['features']:
        name_string, string_id = [], []
        for value in value_dict:
            temp_str = str(value)
            json_value = json.loads(temp_str.replace("'",'"'))
            if json_value['observed'] == 'yes':
                string_id.append(json_value['id'])
                name_string.append(json_value['label'])
        df_phenotips.at[Count,'hpo_all'] = string_id
        df_phenotips.at[Count,'hpo_all_name'] = name_string
        Count = Count + 1
        
    Count = 0
    
        
    name_string, string_id = [], []
    
    for value_dict in df_phenotips['global_age_of_onset']:
        if (type(value_dict) == list):
            name_string, string_id = [], []
            for value in value_dict:
                temp_str = str(value)
                json_value = json.loads(temp_str.replace("'",'"'))
                string_id.append(json_value['id'])
                name_string.append(json_value['label'])
            cur_id, cur_name = df_phenotips.iloc[Count,-2],  df_phenotips.iloc[Count,-1]
            cur_id.extend(string_id)
            cur_name.extend(name_string)
            df_phenotips.at[Count,'hpo_all'] = cur_id
            df_phenotips.at[Count,'hpo_all_name'] = cur_name
            Count = Count + 1
    
    Count = 0
    
    for value_dict in df_phenotips['global_mode_of_inheritance']:
        if (type(value_dict) == list):
            name_string, string_id = [], []
            for value in value_dict:
                temp_str = str(value)
                json_value = json.loads(temp_str.replace("'",'"'))
                string_id.append(json_value['id'])
                name_string.append(json_value['label'])
            cur_id, cur_name = df_phenotips.iloc[Count,-2],  df_phenotips.iloc[Count,-1]
            cur_id.extend(string_id)
            cur_name.extend(name_string)
            df_phenotips.at[Count,'hpo_all'] = cur_id
            df_phenotips.at[Count,'hpo_all_name'] = cur_name
            Count = Count + 1     
        
    df_phenotips = df_phenotips.reset_index(drop=True)    

    for i in range(len(df_phenotips)):
        df_phenotips.at[i, 'hpo_all'] = list(set(df_phenotips.loc[i,'hpo_all']))
        df_phenotips.at[i, 'hpo_all_name'] = list(set(df_phenotips.loc[i,'hpo_all_name']))

    if hpo_cols == True:
        for i in range(len(df_phenotips)):
            for h in df_phenotips.loc[i,'hpo_all']:
                if 'HP:' in h:
                    df_phenotips[h] = df_phenotips.loc[:,'hpo_all'].astype(str).str.contains(h).astype(int)
                    
    df_phenotips = df_phenotips.reset_index(drop=True)
        
    return df_phenotips

df_ankrd = convert_data_phenotips_raw("phenotips-export_2021-11-16.json")
label_list = pd.read_excel("ANKRD11_individualsHPO-forLex.xlsx")
df_ankrd = df_ankrd.merge(label_list.loc[:, ['P-number', 'Variant type']], right_on='P-number', left_on='report_id')
df_ankrd.loc[:, "HP:0004322":].to_csv("hpo_cols_with_labels.csv", index=False)

    
