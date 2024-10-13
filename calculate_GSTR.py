"""
I adapted Emily's GSTR calculation to find S and find RTN and ScoreRTN with bootstrapping.
This is different from scale_metrics_to_inert_base in following ways:
1. S is re-estimated using binary search in each bs cycle
2. p-value is calculated by comparing the bootstrapped statistics to the null distribution
"""
from bootstrapping_helpers import *
from metrics_helpers import *
from GSTR_helpers import *

def bootstrap_GSTR_and_shrinkage(raw_df,input_sample_list1,input_sample_list2,cell_number_cutoff,input_control_gRNA_list,number_of_replicate,input_total_gRNA_number):
    # cell_number_cutoff is for treated (adjusted cutoff). This is given by me
    # find R for RGM calculation
    R_dict = find_ratio_to_inert(raw_df[raw_df['Sample_ID'].isin(input_sample_list2)], input_control_gRNA_list)
    # first find S and untreated (basal) cutoff
    raw_treated_df = Generate_ref_input_df(raw_df,input_sample_list1,0)
    raw_untreated_df = Generate_ref_input_df(raw_df,input_sample_list2,0)
    S, basal_cutoff = find_S(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    # return tumor size metric of each bootstrap cycle
    # applying the cell no. cutoff
    # treated mouse   
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,cell_number_cutoff)
    # untreated mouse
    temp_ref_df2 = Generate_ref_input_df(raw_df,input_sample_list2,basal_cutoff)
    
    # dfs passed are post cutoff
    temp_final_df_observed = calculate_GSTR_metrics(temp_ref_df1, temp_ref_df2, input_control_gRNA_list, R_dict)
    temp_final_df_observed['Shrinkage'] = S
    temp_final_df_observed['Bootstrap_id'] = ['Real'] # a new columns named Bootstrap_id. The values of the column are 'Real'
    
    # for concatenating dfs later
    temp_out_df = [temp_final_df_observed]
    
    if number_of_replicate!=0:
        # experimental
        Mouse_index_dic_1 = Generate_Index_Dictionary(raw_treated_df) # has tumor size info for each SampleID
        # control
        Mouse_index_dic_2 = Generate_Index_Dictionary(raw_untreated_df)
        for bootstrap_cycle in range(number_of_replicate):
            # resampling the experimental mice
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = raw_treated_df.loc[x]
            # resampleing the control mice
            y = Nested_Boostrap_Index_Special_single(Mouse_index_dic_2,temp_ref_df2,input_total_gRNA_number)
            temp_bootstrap_df_2 = raw_untreated_df.loc[y]
            # re-estimate S
            S, basal_cutoff = find_S(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            temp_bootstrap_ref_df1 = Generate_ref_input_df(temp_bootstrap_df_1, temp_bootstrap_df_1['Sample_ID'].unique(), cell_number_cutoff)
            temp_bootstrap_ref_df2 = Generate_ref_input_df(temp_bootstrap_df_2, temp_bootstrap_df_2['Sample_ID'].unique(), basal_cutoff)
            
            temp_metric_df = calculate_GSTR_metrics(temp_bootstrap_ref_df1, temp_bootstrap_ref_df2, input_control_gRNA_list)
            temp_metric_df['Shrinkage'] = S
            temp_metric_df['Bootstrap_id'] = 'B'+str(bootstrap_cycle)
            
            temp_out_df.append(temp_metric_df)
    
    temp_out_df = pd.concat(temp_out_df, ignore_index=True)
    
    return temp_out_df

def calculate_GSTR_metrics(treated_df,untreatd_df,input_control_gRNA_list, ratio_dict):
    # treated and untreated dfs are post cutoff
    treated_sum_df = treated_df.groupby(['gRNA']).Clonal_barcode.count().reset_index(name='TTN')
    untreated_sum_df = untreatd_df.groupby(['gRNA']).Clonal_barcode.count().reset_index(name='TTN')
    ScoreRTN_df = calculate_ScoreRTN(treated_sum_df, untreated_sum_df, input_control_gRNA_list)
    ScoreRGM_df = calculate_ScoreRGM(treated_df, untreatd_df, ratio_dict, input_control_gRNA_list)
    df_merged = ScoreRGM_df.merge(ScoreRTN_df, on='gRNA', how='outer')
    return df_merged