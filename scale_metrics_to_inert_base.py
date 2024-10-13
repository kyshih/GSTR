"""
This script takes in treatment and untreated intermediate df and calculate ScoreRTN and scaled each metric on gRNA level.
I first run bootstrap on cell cutoff L for untreated and LxS for treated groups. Then, I used intermediate df to find RTN
and ScoreRTN.
S is constant here
"""
import pandas as pd
import numpy as np
import argparse
import copy
from scipy.stats import rankdata
from Bootstrapping_analysis.ADJ4_LORSHP2_050824.Python.metrics_helpers import Cal_Bootstrapping_Summary,fdr

def calcualte_scaled_metrics_per_bootstrap(treatment_df, untreated_df, input_control_list):
    treatment_df = treatment_df.copy(deep=True)
    untreated_df = untreated_df.copy(deep=True)
    # Initialize an empty DataFrame to store the final results
    all_scaled_dfs = []
    all_treatment_RTN_dfs = []
    all_untreated_RTN_dfs = []
    
    # Group by bootstrap_cycle
    grouped_treatment = treatment_df.groupby('Bootstrap_id')
    grouped_untreated = untreated_df.groupby('Bootstrap_id')

    for cycle, treatment_group in grouped_treatment:
        untreated_group = grouped_untreated.get_group(cycle)
        # Apply add_cohort_specific_relative_metrics_scaled_to_untreated
        scaled_df = add_cohort_specific_relative_metrics_scaled_to_untreated(treatment_group, untreated_group)
        
        # Apply calculate_ScoreRTN to get the RTN and ScoreRTN for this cycle. res=result
        treatment_res, untreated_res, scaled_res = calculate_ScoreRTN(treatment_group, untreated_group, scaled_df, input_control_list)

        # Add bootstrap cycle identifier to scaled_res
        scaled_res['Bootstrap_id'] = cycle
        treatment_res['Bootstrap_id'] = cycle
        untreated_res['Bootstrap_id'] = cycle
        
        # Append to the final list
        all_scaled_dfs.append(scaled_res)
        all_treatment_RTN_dfs.append(treatment_res[['Targeted_gene_name', 'Numbered_gene_name', 'gRNA', 'TTN_normalized_relative', 'RTN', 'Bootstrap_id']])
        all_untreated_RTN_dfs.append(untreated_res[['Targeted_gene_name', 'Numbered_gene_name', 'gRNA', 'TTN_normalized_relative', 'RTN', 'Bootstrap_id']])
        
    # Combine all results into a single DataFrame
    final_scaled_df = pd.concat(all_scaled_dfs, ignore_index=True)
    final_treatment_RTN_df = pd.concat(all_treatment_RTN_dfs, ignore_index=True)
    final_untreated_RTN_df = pd.concat(all_untreated_RTN_dfs, ignore_index=True)
    
    return final_treatment_RTN_df, final_untreated_RTN_df, final_scaled_df

def Add_Corhort_Specific_Relative_Metrics_Cross_Pool(input_df1,input_df2,input_control_list):
    # Add relative metrics for LN_mean, GEO_mean etc using the median of inert
    # sgTS LN_mean / sgInert LN_mean
    # input_df1 = treatment df
    # input_df2 = untreated df
    temp_sub = input_df2[input_df2['gRNA'].isin(input_control_list)] # untreated inert
    for temp_cname in input_df1.drop(columns=['gRNA'],inplace = False).columns:
        temp_name = temp_cname+'_relative'
        input_df1[temp_name] = input_df1[temp_cname]/temp_sub[temp_cname].median()
        input_df2[temp_name] = input_df2[temp_cname]/temp_sub[temp_cname].median()
        
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
        #log_scaled_value = np.log2(temp_treatment_df[trait] / temp_untreated_df[trait]).tolist()
        #temp_output_df[f'{trait}_scaled_log'] = log_scaled_value
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

    # Ensure the order of gRNA in scaled_df matches with treatment_df
    # reindex treatment_df so division happens as expected or pd give NaN
    scaled_df = scaled_df.set_index('gRNA').reindex(treatment_df['gRNA']).reset_index()
    treatment_df = treatment_df.set_index('gRNA').reindex(scaled_df['gRNA']).reset_index()
    
    # Calculate ScoreRTN
    scaled_df['ScoreRTN'] = np.log2(treatment_df['RTN'] / untreated_df['RTN'])

    return treatment_df, untreated_df, scaled_df
    
def Find_Controls(input_gRNA_df, input_pattern):
# this function will find the gRNA associated with control based on the key word
# input_pattern is a regex expression 
    return(input_gRNA_df.loc[
        input_gRNA_df['Targeted_gene_name'].str.contains(input_pattern, na=False, regex=True),'gRNA'].unique())
    
def Generate_Final_Summary_Dataframe(input_df,trait_of_interest):
    # gRNA level summary
    temp_summary = input_df[input_df['Bootstrap_id']!='Real'].groupby(['gRNA'],as_index = False).apply(Cal_Bootstrapping_Summary,(trait_of_interest))
    temp_output_df = copy.deepcopy(input_df[input_df['Bootstrap_id'] =='Real'])
    temp_output_df = temp_output_df.merge(temp_summary, on = 'gRNA')
    for temp_trait in trait_of_interest:
        temp_name0 = temp_trait + '_fraction_greater_than_one'
        temp_name1 = temp_trait + '_pvalue'
        temp_name2 = temp_name1 + '_FDR'
        temp_name3 = temp_name1 + '_twoside'
        temp_name4 = temp_name1 + '_twoside_FDR'
        if temp_trait == 'ScoreRTN' or 'log' in temp_trait:
            temp_name0 = temp_trait + '_fraction_greater_than_zero'
        temp_output_df[temp_name1] = temp_output_df.apply(lambda x: min(x[temp_name0],1-x[temp_name0]), axis=1) # twosided test
        temp_output_df[temp_name2] = fdr(temp_output_df[temp_name1])
        temp_output_df[temp_name3] = temp_output_df[temp_name1]*2
        temp_output_df[temp_name4] = fdr(temp_output_df[temp_name3])
    return(temp_output_df)


def main():
    parser = argparse.ArgumentParser(description='A function to calcualte the fold change of tumor metrics and ScoreRTN')
    parser.add_argument("--input_untreatd_interm_address")
    parser.add_argument("--input_treatment_interm_address")
    parser.add_argument("--output_file_name")
    #parser.add_argumnet("--output_file_name_gene_level", require=False)
    parser.add_argument("--output_path")
    
    args = parser.parse_args()
    untreatd_interm_address = args.input_untreatd_interm_address
    treatment_interm_address = args.input_treatment_interm_address
    output_path = args.output_path
    out_fname = args.output_file_name
    #out_fname_gene_level = args.output_file_name_gene_level
    
    untreated_interm_df = pd.read_csv(untreatd_interm_address)
    treatment_interm_df = pd.read_csv(treatment_interm_address)
    
    control_gRNA_list = Find_Controls(treatment_interm_df,'Safe|Neo|NT')
    # scaled_trait_list = ['LN_mean_relative_scaled','LN_mean_relative_scaled_log','Geo_mean_relative_scaled','Geo_mean_relative_scaled_log','50_percentile_relative_scaled',
    #                      '50_percentile_relative_scaled_log','95_percentile_relative_scaled','95_percentile_relative_scaled_log','TTB_normalized_relative_scaled',
    #                      'TTB_normalized_relative_scaled_log','TTN_normalized_relative_scaled','TTN_normalized_relative_scaled_log','TTN_scaled','TTN_scaled_log',
    #                      'ScoreRTN']
    scaled_trait_list = ['LN_mean_relative_scaled','Geo_mean_relative_scaled','50_percentile_relative_scaled',
                         '95_percentile_relative_scaled','TTB_normalized_relative_scaled', 'TTN_normalized_relative_scaled',
                         'TTN_scaled','ScoreRTN']
    print(f'calculating scaled metrics')
    
    treatment_RTN_interm_df, untreated_RTN_interm_df, scaled_interm_df = calcualte_scaled_metrics_per_bootstrap(treatment_interm_df, untreated_interm_df, control_gRNA_list)
    print(f'generating bs summary for the scaled df')
    scaled_bs_summary_df = Generate_Final_Summary_Dataframe(scaled_interm_df,scaled_trait_list) # gRNA level summary
    
    print(f'generating gene level intermediate for treatment and untreated groups')
    #untreated_RTN_interm_df_gene_level = untreated_RTN_interm_df.groupby(['Targeted_gene_name','Bootstrap_id'],as_index = False).apply(Cal_Combined_Gene_Effect_v2,['RTN'])
    #treatment_RTN_interm_df_gene_level = treatment_RTN_interm_df.groupby(['Targeted_gene_name','Bootstrap_id'],as_index = False).apply(Cal_Combined_Gene_Effect_v2,['RTN'])
    #treatment_RTN_interm_df_gene_level = treatment_RTN_interm_df_gene_level.rename(columns={'RTN': 'RTN_treatment'})
    #untreated_RTN_interm_df_gene_level = untreated_RTN_interm_df_gene_level.rename(columns={'RTN': 'RTN_untreated'})
    #untreated_RTN_interm_df_gene_level['ScoreRTN'] = treatment_RTN_interm_df_gene_level['RTN_treatment'] / treatment_RTN_interm_df_gene_level['RTN_untreated']
    #untreated_RTN_interm_df_gene_level.to_csv(f'{output_path}/{out_fname_gene_level}_intermediate', index=None)
    #TODO: gene level ScoreRTN summary
    
    # saving df
    scaled_interm_df.to_csv(f'{output_path}/{out_fname}_intermediate',index=None)
    scaled_bs_summary_df.to_csv(f'{output_path}/{out_fname}_result.csv',index=None)
    
if __name__ == '__main__':
    main()
