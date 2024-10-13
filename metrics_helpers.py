import pandas as pd
import numpy as np
import math
from scipy.stats import rankdata

''' This module contains functions to calcualte tumor metrics for each bs cycle. 
    It also normalizes TTN and TTB to KT mice and calculate metrics relative to inert tumors. 
'''

def Cal_Combined_Gene_Effect_v2(x,trait_of_interest):
    d = {}
    temp_weight_list = x['TTN_normalized_relative']
    for temp_trait in trait_of_interest:
        d[temp_trait] = sum(x[temp_trait]*temp_weight_list/sum(temp_weight_list)) # weighted sum by normalized TTN relative to inert
    return pd.Series(d, index=list(d.keys())) 

def Calculate_Relative_Normalized_Metrics(input_df1,input_df2,percentile_list,input_control_gRNA_list):
    # Cas9 mouse
    # print('haha'+str(input_df1.shape[0]))
    temp_df = input_df1.reset_index().groupby(['gRNA'],as_index = False).apply(Cal_Tumor_Size_simple,(percentile_list))
    # print(input_df1[['gRNA','Cell_number']])
    # print(temp_df)
    # Control mouse
    temp_df2 = input_df2.reset_index().groupby(['gRNA'],as_index = False).apply(Cal_Tumor_Size_simple,(percentile_list))
    
    # normalize TTN and TTB to KT
    temp_out = Generate_Normalized_Metrics(temp_df,temp_df2,['TTN','TTB'])
    temp_df = temp_df.merge(temp_out,on ='gRNA') # merge

    # calculate relative expression to sgInert
    Add_Corhort_Specific_Relative_Metrics(temp_df,input_control_gRNA_list)

    # annotate sample type
    temp_df['Type'] = temp_df.apply(lambda x: 'Inert' if (x['gRNA'] in input_control_gRNA_list) else 'Experiment',axis=1)
    temp_df = temp_df.merge(input_df1[['gRNA','Targeted_gene_name',
     'Identity', 'Numbered_gene_name']].drop_duplicates(),how = 'inner',on = 'gRNA')
    return(temp_df)

def LN_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    temp_var = log_vector.var()
    if len(log_vector)==1:
        temp_var = 0 # if only one clonal
    return (math.exp(temp_mean + 0.5*temp_var))

# calculate the Geometric mean from a vector of number
def Geometric_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    return (math.exp(temp_mean))

def fdr(p_vals):
    p = np.asfarray(p_vals) # make input as float array
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    p = p[by_descend] # sort pvalue from small to large
    ranked_p_values = rankdata(p,method ='max') # this max is very important, when identical, use largest
    fdr = p * len(p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    return fdr[by_orig]

def Cal_Tumor_Size_simple(x,input_percentile):
    # return 1D labeled array for tumor size. similar to a column
    d = {}
    temp_vect = x['Cell_number']
    if type (temp_vect) == 'int':
        temp_vect = [temp_vect]
    # measure size
    d['LN_mean'] = LN_Mean(temp_vect)
    d['Geo_mean'] = Geometric_Mean(temp_vect)
    Percentile_list = list(np.percentile(temp_vect,input_percentile))
    for c,y in enumerate(input_percentile):
        temp_name = str(y)+'_percentile'
        d[temp_name] = Percentile_list[c]
    d['TTN'] = len(temp_vect) # this is total tumor number
    d['TTB'] = sum(temp_vect)
    return pd.Series(d, index=list(d.keys()))  

def Generate_Normalized_Metrics(input_df1,input_df2,trait_list):
    # this functional use input_df2 to normalized input_df1 using metrics defined by trait_list 
    # input_df1 is the experimental group (KTHC), input_df2 is the control group (KT)
    temp1 = input_df1.set_index('gRNA')
    temp2 = input_df2.set_index('gRNA').loc[temp1.index] # only looking at gRNA sampled in KTHC. temp1&2 gRNAs are in the same order
    temp_output_df = pd.DataFrame({'gRNA':temp1.index.values})
    for temp_cname in trait_list:
        temp_cname_new = temp_cname + '_normalized'
        temp_output_df[temp_cname_new] = np.array(temp1[temp_cname].to_list())/np.array(temp2[temp_cname].to_list())
    return(temp_output_df)

def Add_Corhort_Specific_Relative_Metrics(input_df,input_control_list):
    # Add relative metrics for LN_mean, GEO_mean etc using the median of inert
    temp_sub = input_df[input_df['gRNA'].isin(input_control_list)]
    #for temp_cname in input_df.drop(columns=['gRNA'],inplace = False).columns: # use this for more traits of interest
    trait_of_interst = input_df.drop(columns=['gRNA'],inplace = False).columns
    #trait_of_interst = ['LN_mean', 'Geo_mean', '50_percentile', '95_percentile', 'TTB_normalized', 'TTN_normalized', 'TTN']
    for temp_cname in trait_of_interst:
        # Check if the column is numeric before proceeding
        if pd.api.types.is_numeric_dtype(input_df[temp_cname]):
            temp_name = temp_cname+'_relative'
            # Calculate the relative metric using the median of the inert controls
            input_df[temp_name] = input_df[temp_cname]/temp_sub[temp_cname].median()
        else:
            print(f"Skipping {temp_cname} as it is not numeric.")
        
def add_cohort_specific_relative_metrics_scaled_to_untreated(treatment_specific_df, untreated_df):
    # this version is for balanced gRNA representation. every gRNA is present in each bs cycle
    # Rt = sgTS/sgInert in treatment group. Rc = sgTS/sgInert in untreated group
    # scaled metric = Rt/Rc
    trait_list = ['LN_mean_relative', 'Geo_mean_relative', '50_percentile_relative', '95_percentile_relative', 'TTB_normalized_relative',
                'TTN_normalized_relative', 'TTN']
    temp_output_df = treatment_specific_df[['gRNA','Targeted_gene_name', 'Numbered_gene_name','TTN', 'TTN_normalized_relative']].copy(deep=True)
    temp_treatment_df = treatment_specific_df.set_index('gRNA')
    temp_untreated_df = untreated_df.set_index('gRNA').reindex(temp_treatment_df.index) # Align rows in treatement and untreated df
    for trait in trait_list:
        scaled_value = (temp_treatment_df[trait] / temp_untreated_df[trait]).tolist()
        temp_output_df[f'{trait}_scaled'] = scaled_value
        # log_scaled_value = np.log2(temp_treatment_df[trait] / temp_untreated_df[trait]).tolist()
        # temp_output_df[f'{trait}_scaled_log'] = log_scaled_value
    return temp_output_df

def calculate_ScoreRTN(treatment_df, untreated_df, scaled_df, input_control_list):
    # Calculate TTN_inert for treatment and untreated
    TTN_inert_treatment = treatment_df[treatment_df['gRNA'].isin(input_control_list)]['TTN'].sum()
    TTN_inert_untreated = untreated_df[untreated_df['gRNA'].isin(input_control_list)]['TTN'].sum()
    
    # Align the indices of untreated_df with treatment_df
    untreated_df = untreated_df.set_index('gRNA').reindex(treatment_df['gRNA']).reset_index()

    # Calculate RTN for treatment and untreated
    treatment_df['RTN'] = treatment_df['TTN'] / TTN_inert_treatment
    untreated_df['RTN'] = untreated_df['TTN'] / TTN_inert_untreated

    # Ensure the order of gRNA in scaled_df
    scaled_df = scaled_df.set_index('gRNA').reindex(treatment_df['gRNA']).reset_index()
    
    # Calculate ScoreRTN
    scaled_df['ScoreRTN'] = np.log2(treatment_df['RTN'] / untreated_df['RTN'])

    return treatment_df, untreated_df, scaled_df

    