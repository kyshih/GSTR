"""
    This contains GSTR helpers
"""
import pandas as pd
import numpy as np
import math
from scipy.stats.mstats import gmean
from find_S_helpers import find_S_given_base_cutoff, find_S

def compute_top_N_per_sample(df, ratio_dict, input_control_gRNA_list):
    '''
    Generates a nested dictionary mapping each Sample_ID to a dictionary of gRNAs and their scaled TTN counts.
    
    Params:
    - df: df containing the data with columns ['Sample_ID', 'gRNA', 'Clonal_barcode']
    - ratio_dict (dict): Dictionary mapping gRNA to median TTN_to_inert_ratio, e.g., {gRNA1: ratio1, gRNA2: ratio2, ...}
    - input_control_gRNA_list (list): List of gRNAs considered as inert controls
    
    Returns:
    - dict: Nested dictionary structured as {Sample_ID: {gRNA: scaled_TTN, ...}, ...}
    '''
    top_N_tumors_dict = {}
    
    # Precompute the number of inert tumors per Sample_ID
    inert_counts = df[df['gRNA'].isin(input_control_gRNA_list)].groupby('Sample_ID').size()
    #print(f"{inert_counts}")
    # Iterate over each unique Sample_ID
    for sample in df['Sample_ID'].unique():
        # Get the number of inert tumors for the current sample
        TTN_inert_per_mouse = inert_counts.get(sample)
        if not TTN_inert_per_mouse:
            print(f"{sample} is missing inert tumors") # inert tumors in this mouse are below cutoff
            continue
        # Scale the number of inert tumors by the ratio to inert for each gRNA
        scaled_TTN = {gRNA: ratio * TTN_inert_per_mouse for gRNA, ratio in ratio_dict.items()}
        # Assign the scaled TTN dictionary to the current sample
        top_N_tumors_dict[sample] = scaled_TTN
        
    return top_N_tumors_dict

def find_ratio_to_inert(untreated_df, input_control_gRNA_list, count_cutoff=2):
    ''' 
    the goal is to correct for mouse-to-mouse variability in calcualting RGM
    by finding the ratio to inert per mouse then finding the median of the ratio for each gRNA
    return a ratio dict {gRNAi: Ri}
    '''
    untreated_df = untreated_df[(untreated_df['Count'] > count_cutoff)].copy(deep=True)
    # Find the total number of inert gRNAs per mouse (TTN_inert)
    TTN_inert_per_mouse = untreated_df[untreated_df['gRNA'].isin(input_control_gRNA_list)].groupby(['Sample_ID']).gRNA.count().reset_index(name='TTN_inert')
    # Find the total number of each gRNAs per mouse (TTN), both inerts and non-inerts
    TTN_df = untreated_df.groupby(['gRNA', 'Sample_ID']).Clonal_barcode.count().reset_index(name='TTN')
    # Merge TTN with the inert gRNA counts per mouse
    TTN_df = TTN_df.merge(TTN_inert_per_mouse, on='Sample_ID', how='left')
    # Calculate the ratio of TTN to TTN_inert
    TTN_df['TTN_to_inert_ratio'] = TTN_df['TTN'] / TTN_df['TTN_inert']
    # Group by gRNA and find the median of the TTN to inert ratios
    ratio_dict = TTN_df.groupby(['gRNA']).apply(lambda x: np.median(x['TTN_to_inert_ratio'])).to_dict()
    
    return ratio_dict

def subset_top_N_per_gRNA_across_samples(df, top_N_dict):
    """
    Args:
        df: post-cutoff df with gRNA, Clonal_barcode, Cell_number, etc
        top_N_dict: desire top N tumor counts of gRNAi in mouse k

    Returns:
        df with top N largets tumors of gRNAi in all mice
    """
    result = []
    
    # Iterate through each Sample_ID
    for sample_id, gRNA_dict in top_N_dict.items():
        # Filter for the current Sample_ID
        df_sample = df[df['Sample_ID'] == sample_id]
        
        # Iterate through each gRNA and the corresponding top N value
        for gRNA, top_n in gRNA_dict.items():
            # Filter for the current gRNA
            df_gRNA = df_sample[df_sample['gRNA'] == gRNA]
            
            # Sort the entries by Cell_number in descending order
            df_gRNA_sorted = df_gRNA.sort_values(by='Cell_number', ascending=False)
            
            # if this gRNA is not found in this mouse just skip the rest of the current iterationa and move onto the next gRNA
            if df_gRNA_sorted.empty:
                continue
            # Sort the entries by Cell_number in descending order and take the top N without duplicating the smallest tumors
            #top_entries = df_gRNA.nlargest(math.ceil(top_n), 'Cell_number')
            
            # Check if the number of tumors is less than the desired top N
            if len(df_gRNA_sorted) < top_n:
                # Get the smallest tumor in the gRNA
                smallest_tumor = df_gRNA_sorted.iloc[-1]
                # Calculate how many times to duplicate the smallest tumor
                num_to_add = math.ceil(top_n) - len(df_gRNA_sorted)
                # Duplicate the smallest tumor to reach top N
                duplicates = pd.DataFrame([smallest_tumor] * num_to_add)
                # Append the duplicates to the original DataFrame
                df_gRNA_sorted = pd.concat([df_gRNA_sorted, duplicates], ignore_index=True)
            top_entries = df_gRNA_sorted.head(math.ceil(top_n))
            # Append the selected top entries to the result
            result.append(top_entries)
    
    # Concatenate the results into a final DataFrame
    return pd.concat(result)

def calculate_ScoreRTN(treated_df, untreated_df, input_control_gRNA_list):
    """
	Compare the # of tumors above the adjusted cutoff in treated to untreated.
	The adaptive nature of the cutoff takes into account the global size shift associated
	with treatment. Per usual, adjust gRNA i to inerts, and then treated to untreated
 
    Params:
    treated_df, untreated_df: post cutoff
    
    append RTN and ScoreRTN columns
	"""
    TTN_inert_treatment = treated_df[treated_df['gRNA'].isin(input_control_gRNA_list)]['TTN'].sum()
    TTN_inert_untreated = untreated_df[untreated_df['gRNA'].isin(input_control_gRNA_list)]['TTN'].sum()
    # Align the indices of untreated_df with treatment_df
    # dont need this cuz I merge the two df later
    #untreated_df = untreated_df.set_index('gRNA').reindex(treated_df['gRNA']).reset_index()
    
    # Calculate RTN for treatment and untreated
    treated_df['TTN_inert_treated'] = TTN_inert_treatment
    treated_df['RTN_treated'] = treated_df['TTN'] / TTN_inert_treatment
    
    untreated_df['TTN_inert_untreated'] = TTN_inert_untreated
    untreated_df['RTN_untreated'] = untreated_df['TTN'] / TTN_inert_untreated
    
    treated_df = treated_df.merge(untreated_df, on='gRNA', suffixes=('_treated', '_untreated'))
    #treated_df['ScoreRTN'] = np.log2(treated_df['RTN'] / untreated_df['RTN']) # need to align df if doing this
    treated_df['ScoreRTN'] = np.log2(treated_df['RTN_treated'] / treated_df['RTN_untreated'])
    return treated_df

def compute_geo_mean(df, input_control_gRNA_list, group_cols, col_name):
    # Calculate geometric mean for inert gRNAs
    gm_inert = gmean(df[df['gRNA'].isin(input_control_gRNA_list)]['Cell_number'])
    
    # Group by the necessary columns and calculate the geometric mean
    gm_df = df.groupby(group_cols).apply(lambda x: gmean(x['Cell_number'])).reset_index(name=col_name)
    
    # Add the geometric mean for the inert controls
    gm_df[f'{col_name}_inert'] = gm_inert
    
    return gm_df

def prepare_RGM_data(treated_df_cut, untreated_df_cut, ratio_dict, input_control_gRNA_list):
    # Prepare treated DataFrame
    top_N_treated_dict = compute_top_N_per_sample(treated_df_cut, ratio_dict, input_control_gRNA_list)
    top_N_treated_df = subset_top_N_per_gRNA_across_samples(treated_df_cut, top_N_treated_dict)
    
    # Prepare untreated DataFrame
    top_N_untreated_dict = compute_top_N_per_sample(untreated_df_cut, ratio_dict, input_control_gRNA_list)
    top_N_untreated_df = subset_top_N_per_gRNA_across_samples(untreated_df_cut, top_N_untreated_dict)
    
    return top_N_treated_df, top_N_untreated_df

def calculate_ScoreRGM(treated_df_cut, untreated_df_cut, ratio_dict, input_control_gRNA_list):
    # Step 1: Prepare Data
    top_N_treated_df, top_N_untreated_df = prepare_RGM_data(treated_df_cut, untreated_df_cut, ratio_dict, input_control_gRNA_list)
    
    # Step 2: Compute Geometric Means
    gm_treated_df = compute_geo_mean(top_N_treated_df, input_control_gRNA_list, 
                                           ['gRNA', 'Targeted_gene_name', 'Numbered_gene_name'], 
                                           'Geo_mean_treated')
    gm_untreated_df = compute_geo_mean(top_N_untreated_df, input_control_gRNA_list, 
                                             ['gRNA'], 
                                             'Geo_mean_untreated')
    
    # Step 3: Merge the treated and untreated geometric mean DataFrames
    gm_treated_df = gm_treated_df.merge(gm_untreated_df, on='gRNA')
    
    # Step 4: Calculate RGM
    gm_treated_df['RGM_treated'] = gm_treated_df['Geo_mean_treated'] / gm_treated_df['Geo_mean_treated_inert']
    gm_treated_df['RGM_untreated'] = gm_treated_df['Geo_mean_untreated'] / gm_treated_df['Geo_mean_untreated_inert']

    # Step 4: Calculate ScoreRGM
    gm_treated_df['ScoreRGM'] = np.log2(gm_treated_df['RGM_treated'] / gm_treated_df['RGM_untreated'])
    
    return gm_treated_df

def df_to_dict(df, key_col, value_col):
    """
    Convert specified DataFrame columns into a dictionary.

    Parameters:
    - df: to convert.
    - key_col (str): The column to use as the dictionary's keys.
    - value_col (str): The column to use as the dictionary's values.

    Returns:
    - dict: A dictionary mapping `key_col` to `value_col`.
    """
    return df.set_index(key_col)[value_col].to_dict()

def Nested_Bootstrap_Index_single_vary_mouse_count(input_dic, mouse_number):
    """resample mice based on their indices

    Args:
        input_dic: {SampleID : [row_number that corresponds to a gRNA and read counts, etc]}
        mouse_number: no. of mice you want to resample
    """
    temp_sample_list = list(input_dic.keys()) # list of SampleIDs
    # I first sample mouse
    temp_list = np.random.choice(temp_sample_list,mouse_number,replace = True) # sample the SampleID with replacement
    temp_coho = []
    for y in temp_list: # within each mouse
        temp_array = input_dic.get(y) # get index of gRNA read associated of that mouse. array of tuple, each is a (gRNA, clonal_barcode)
        temp_resampled = np.random.choice(temp_array,len(temp_array),replace = True) # resample gRNA
        temp_coho = np.concatenate([temp_coho,temp_resampled])
    return(temp_coho)  