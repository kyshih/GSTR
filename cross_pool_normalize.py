import pandas as pd
import numpy as np
import argparse
import copy
from scipy.stats import rankdata
from UltraSeq_Bootstrapping_EA_drug import Generate_Final_Summary_Dataframe,Generate_Gene_Level_Summary_Dataframe,Cal_Bootstrapping_Summary,Cal_Combined_Gene_Effect_v2

def add_cohort_specific_relative_metrics_cross_pool(input_df1, input_df2, input_control_list):
    # Make copies of the dataframes to avoid altering the original data
    input_df1 = input_df1.copy(deep=True)
    input_df2 = input_df2.copy(deep=True)
    
    # Filter the untreated DataFrame to get the inert guides. normalization factor is the median of group2-usually inert base or untreated
    temp_sub = input_df2[input_df2['gRNA'].isin(input_control_list)] 
    
    columns = ['LN_mean', 'Geo_mean', '50_percentile', '60_percentile', '70_percentile', '80_percentile', '90_percentile', '95_percentile',
                  '97_percentile', '99_percentile', 'TTN', 'TTB', 'TTN_normalized', 'TTB_normalized']
    # Normalize each column in input_df1 and input_df2
    # for temp_cname in input_df1.drop(columns=['gRNA', 'Type', 'Targeted_gene_name', 
    #                                           'Identity', 'Numbered_gene_name', 'Bootstrap_id'], inplace=False).columns:
        
    for temp_cname in columns:
        temp_name = temp_cname + '_relative'
        
        # Calculate the median of the current column in untreated inert guides for the current Bootstrap_id
        median_inert_value = temp_sub[temp_cname].median()
        print(f'{temp_cname} median inert value is {median_inert_value}')
        
        # Normalize both the treatment and untreated DataFrames using the median value
        input_df1[temp_name] = input_df1[temp_cname] / median_inert_value
        input_df2[temp_name] = input_df2[temp_cname] / median_inert_value
    
    return input_df1, input_df2

def calculate_cross_pool_normalized_metrics_per_bootstrap(input_df1, input_df2, input_control_list):
    """
    Normalize metrics across pools for each bootstrap cycle and return the combined results.

    @param input_df1: 
        df with the first set of data (e.g., treatment group) with a 'Bootstrap_id' column.
        
    @param input_df2: 
        df with the second set of data (e.g., untreated group) with a 'Bootstrap_id' column.
        
    @param input_control_list:
        A list of control gRNAs to be used for normalization purposes.

    @return: tuple (pd.DataFrame, pd.DataFrame)
        - final_normalized_input_df1: The normalized DataFrame corresponding to `input_df1`, with an additional 'Bootstrap_id' column.
        - final_normalized_input_df2: The normalized DataFrame corresponding to `input_df2`, with an additional 'Bootstrap_id' column.
        
    This function performs the following steps:
    1. Groups the input DataFrames (`input_df1` and `input_df2`) by `Bootstrap_id`.
    2. For each bootstrap cycle, it normalizes the metrics in both DataFrames using the control gRNAs specified in `input_control_list`.
    3. The normalized results for each bootstrap cycle are appended to a list.
    4. The lists are concatenated to form the final normalized DataFrames, which are returned as a tuple.
    """
    # Initialize lists to store the final results
    normalized_input_df1s = []
    normalized_input_df2s = []

    # Group by Bootstrap_id
    group1 = input_df1.groupby('Bootstrap_id')
    group2 = input_df2.groupby('Bootstrap_id')

    for cycle, sub_group1 in group1:
        sub_group2 = group2.get_group(cycle)
        
        #print(f'{cycle}')
        # Apply the normalization for the current bootstrap cycle. all are normalized to sub_group2
        normalized_input_df1, normalized_input_df2 = add_cohort_specific_relative_metrics_cross_pool(
            sub_group1, sub_group2, input_control_list)
        
        # Add bootstrap cycle identifier to the normalized DataFrames
        normalized_input_df1['Bootstrap_id'] = cycle
        normalized_input_df2['Bootstrap_id'] = cycle
        
        # Append to the final lists
        normalized_input_df1s.append(normalized_input_df1)
        normalized_input_df2s.append(normalized_input_df2)

    # Combine all results into single DataFrames
    final_normalized_input_df1 = pd.concat(normalized_input_df1s, ignore_index=True)
    final_normalized_input_df2 = pd.concat(normalized_input_df2s, ignore_index=True)
    
    return final_normalized_input_df1, final_normalized_input_df2

def Find_Controls(input_gRNA_df, input_pattern):
# this function will find the gRNA associated with control based on the key word
# input_pattern is a regex expression 
    return(input_gRNA_df.loc[
        input_gRNA_df['Targeted_gene_name'].str.contains(input_pattern, na=False, regex=True),'gRNA'].unique())

def main():
    parser = argparse.ArgumentParser(description='A function to calcualte fold change of tumor metrics')
    parser.add_argument("--input_treatment_interm_address", required=True, help="Address of treatment or gene ko base df for cross pool normalization")
    parser.add_argument("--input_untreatd_interm_address", required=True, help="Address of untreated or inert base df address used for normalization factor")
    parser.add_argument("--output_file_name1", required=True, help="output treatment file name prefix")
    parser.add_argument("--output_file_name2", required=False, help="output untreated file name prefix")
    parser.add_argument("--output_path", required=True, help="This the output address for summary data")
    
    args = parser.parse_args()
    untreatd_interm_address = args.input_untreatd_interm_address
    treatment_interm_address = args.input_treatment_interm_address
    output_path = args.output_path
    out_fname1 = args.output_file_name1
    
    if args.output_file_name2:
        out_fname2 = args.output_file_name2
    
    temp_q = [50, 60, 70, 80, 90, 95, 97, 99]
    
    untreated_interm_df = pd.read_csv(untreatd_interm_address)
    treatment_interm_df = pd.read_csv(treatment_interm_address)
    
    control_gRNA_list = Find_Controls(treatment_interm_df,'Safe|Neo|NT')
    print(f'calculating cross prool normalization')
    normalized_treatment_interm_df, normalized_untreated_interm_df = calculate_cross_pool_normalized_metrics_per_bootstrap(
                                                                treatment_interm_df, untreated_interm_df, control_gRNA_list)
    
    temp_trait_list = ['LN_mean_relative','Geo_mean_relative','TTB_normalized_relative','TTN_normalized_relative','95_percentile_relative'] + [str(x) + '_percentile_relative' for x in temp_q]
    temp_trait_list = list(set(temp_trait_list))
    print(f'generating gRNA level bs summary')
    final_treatment_bs_summary_df = Generate_Final_Summary_Dataframe(normalized_treatment_interm_df,temp_trait_list)
    #final_untreated_bs_summary_df = Generate_Final_Summary_Dataframe(normalized_untreated_interm_df, temp_trait_list)
    print(f'generating gene level bs summary')
    final_treatment_bs_gene_summary_df = Generate_Gene_Level_Summary_Dataframe(normalized_treatment_interm_df,temp_trait_list)
    #final_untreated_bs_gene_summary_df = Generate_Gene_Level_Summary_Dataframe(normalized_untreated_interm_df,temp_trait_list)
    
    normalized_treatment_interm_df.to_csv(f'{output_path}/{out_fname1}_intermediate',index=None)
    #normalized_untreated_interm_df.to_csv(f'{output_path}/{out_fname2}_intermediate',index=None)

    
    final_treatment_bs_summary_df.to_csv(f'{output_path}/{out_fname1}_result.csv',index=None)
    #final_untreated_bs_summary_df.to_csv(f'{output_path}/{out_fname2}_result.csv',index=None)
    
    final_treatment_bs_gene_summary_df.to_csv(f'{output_path}/{out_fname1}_result_gene_level.csv',index=None)
    #final_untreated_bs_gene_summary_df.to_csv(f'{output_path}/{out_fname2}_result_gene_level.csv',index=None)
    
if __name__ == '__main__':
    main()