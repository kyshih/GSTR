import pandas as pd
import numpy as np
from scipy.stats import rankdata
import copy
from metrics_helpers import *

"""
    This module contains helpers for summarziing bs results
"""

def Generate_Final_Summary_Dataframe(input_df,trait_of_interest): # gRNA level
    temp_summary = input_df[input_df['Bootstrap_id']!='Real'].groupby('gRNA',as_index = False).apply(Cal_Bootstrapping_Summary,(trait_of_interest))
    temp_output_df = copy.deepcopy(input_df[input_df['Bootstrap_id'] =='Real'])
    temp_output_df = temp_output_df.merge(temp_summary, on = 'gRNA')
    for temp_trait in trait_of_interest:
        if temp_trait == 'ScoreRTN' or temp_trait == 'ScoreRGM':
            temp_name0 = temp_trait + '_fraction_greater_than_zero' # t_test pvalue column name
        else:
            temp_name0 = temp_trait + '_fraction_greater_than_one' # t_test pvalue column name
        temp_name1 = temp_trait + '_pvalue'
        temp_name2 = temp_name1 + '_FDR'
        temp_name3 = temp_name1 + '_twoside'
        temp_name4 = temp_name1 + '_twoside_FDR'
        temp_output_df[temp_name1] = temp_output_df.apply(lambda x: min(x[temp_name0],1-x[temp_name0]), axis=1)
        temp_output_df[temp_name2] = fdr(temp_output_df[temp_name1])
        temp_output_df[temp_name3] = temp_output_df[temp_name1]*2
        temp_output_df[temp_name4] = fdr(temp_output_df[temp_name3])
    return(temp_output_df)

def Generate_Gene_Level_Summary_Dataframe(input_df, trait_of_interest):
    # calculate gene effect from bootstrapping to estiamte CI, pval, and FDR
    temp_df = input_df[input_df['Bootstrap_id']!='Real'].groupby([
        'Targeted_gene_name','Bootstrap_id'],as_index = False).apply(Cal_Combined_Gene_Effect_v2,(trait_of_interest))
    temp_summary = temp_df.groupby('Targeted_gene_name',as_index = False).apply(Cal_Bootstrapping_Summary,(trait_of_interest))

    # calculate gene effect for the observed data. sample estiamte for the population parameters
    temp_output_df = input_df[input_df['Bootstrap_id'] =='Real'].groupby([
        'Targeted_gene_name'],as_index = False).apply(Cal_Combined_Gene_Effect_v2,(trait_of_interest))
    temp_output_df = temp_output_df.merge(temp_summary, on = 'Targeted_gene_name')
    for temp_trait in trait_of_interest:
        if temp_trait == 'ScoreRTN' or temp_trait == 'ScoreRGM':
            temp_name0 = temp_trait + '_fraction_greater_than_zero' # t_test pvalue column name
        else:
            temp_name0 = temp_trait + '_fraction_greater_than_one' # t_test pvalue column name
        temp_name1 = temp_trait + '_pvalue'
        temp_name2 = temp_name1 + '_FDR'
        temp_name3 = temp_name1 + '_twoside'
        temp_name4 = temp_name1 + '_twoside_FDR'
        temp_output_df[temp_name1] = temp_output_df.apply(lambda x: min(x[temp_name0],1-x[temp_name0]), axis=1) # 
        temp_output_df[temp_name2] = fdr(temp_output_df[temp_name1])
        temp_output_df[temp_name3] = temp_output_df[temp_name1]*2
        temp_output_df[temp_name4] = fdr(temp_output_df[temp_name3])
    return(temp_output_df)

def Cal_Bootstrapping_Summary(x,trait_of_interest):
    d = {}
    for temp_trait in trait_of_interest:
        temp0 = temp_trait + '_95P'
        temp1 = temp_trait + '_5P'
        temp2 = temp_trait +'_fraction_greater_than_one' # t_test pvalue column name
        temp3 = temp_trait +'_bootstrap_median'
        temp4 = temp_trait +'_bootstrap_mean'
        temp5 = temp_trait + '_97.5P'
        temp6 = temp_trait + '_2.5P'
        d[temp0] = x[temp_trait].quantile(0.95)
        d[temp1] = x[temp_trait].quantile(0.05)
        if temp_trait == 'ScoreRTN' or temp_trait == 'ScoreRGM':
            temp2 = temp_trait +'_fraction_greater_than_zero'
            d[temp2] = sum(x[temp_trait]>0)/len(x[temp_trait])
        else:
            d[temp2] = sum(x[temp_trait]>1)/len(x[temp_trait])
        d[temp3] = x[temp_trait].mean()
        d[temp4] = x[temp_trait].median()
        d[temp5] = x[temp_trait].quantile(0.975)
        d[temp6] = x[temp_trait].quantile(0.025)
    #creates a pd Series obj containing the computed statistics, with the index being the column names ie. temp2, temp3, etc
    return pd.Series(d, index=list(d.keys())) 

def Generate_Index_Dictionary(input_df):
    # This function generate a dictionary to speed up the boostrap process
    temp_dic = {}
    temp_group = input_df.groupby(['Sample_ID'])
    # iterate through each group
    for key in temp_group.groups.keys(): # Get the keys (unique 'Sample_ID')
        # For each 'Sample_ID', get the array of indices and assign it to the dictionary
        temp_dic[key] = temp_group.get_group(key).index.values
        #temp_dic[key] = temp_group.get_group((key,)).index.values  # Wrap the key in a tuple
    return(temp_dic)

def Generate_ref_input_df(input_df,input_sample_list,input_cell_cutoff):
    return(input_df[(input_df['Cell_number']>input_cell_cutoff)&(input_df['Sample_ID'].isin(input_sample_list))])

def Cal_Combined_Gene_Effect_v2(x,trait_of_interest): 
    # weighted effect of each gRNA based on TTN to see the combined gene effect
    d = {}
    temp_weight_list = x['RTN_treated'] # for GSTR
    #temp_weight_list = x['TTN_normalized_relative'] # using TTN as the weights associated with each gRNA
    #temp_weight_list = x['TTN_normalized']
    for temp_trait in trait_of_interest: # loop through each tumor metric, e.g LN_mean, Geo_mean, etc
        # normalizes the weighted effects by dividing each term by the sum of weights
        # then calculates the sum of the normalized, weighted effects
        d[temp_trait] = sum(x[temp_trait]*temp_weight_list/sum(temp_weight_list))
    return pd.Series(d, index=list(d.keys())) 

def Find_Controls(input_gRNA_df, input_pattern):
    # this function will find the gRNA associated with control based on the key word
    # input_pattern is a regex expression 
    return(input_gRNA_df.loc[
        input_gRNA_df['Targeted_gene_name'].str.contains(input_pattern, na=False, regex=True),'gRNA'].unique())
    
def Nested_Boostrap_Index_single(input_dic):
    # input_dic has {SampleID : [row_number that corresponds to a gRNA and read counts, etc]}
    temp_sample_list = list(input_dic.keys()) # list of SampleIDs
    # I first sample mouse
    temp_list = np.random.choice(temp_sample_list,len(temp_sample_list),replace = True) # sample the SampleID with replacement
    temp_coho = []
    for y in temp_list: # within each mouse
        temp_array = input_dic.get(y) # get index of gRNA read associated of that mouse. array of tuple, each is a (gRNA, clonal_barcode)
        temp_resampled = np.random.choice(temp_array,len(temp_array),replace = True) # resample gRNA
        temp_coho = np.concatenate([temp_coho,temp_resampled])
    return(temp_coho)  

def Nested_Boostrap_Index_Special_single(input_dic,input_df,input_total_gRNA_number): # for the control mice
    temp_sample_list = list(input_dic.keys())
    # I first sample mouse
    temp_coho = []
    while len(set(input_df.loc[temp_coho].gRNA)) < input_total_gRNA_number: # stop until we reamples all the gRNA in the control mice. usually KT
        temp_list = np.random.choice(temp_sample_list,len(temp_sample_list),replace = True)
        temp_coho = []
        for y in temp_list: # within each mouse
            if y not in input_dic:
                print(f"Mouse '{y}' not found in the mouse index dictionary.")
                continue # skip to the next mouse if one mouse is not found in the dictionary
            temp_array = input_dic.get(y) # array of tuple, each is a (gRNA, clonal_barcode)
            temp_resampled = np.random.choice(temp_array,len(temp_array),replace = True)
            temp_coho = np.concatenate([temp_coho,temp_resampled]) 
    return(temp_coho)  