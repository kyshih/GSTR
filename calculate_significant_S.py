""" This module finds significant shrinkage (S) by bootstrapping tumors across cutoffs
"""

import pandas as pd
from bootstrapping_helpers import Generate_ref_input_df, Generate_Index_Dictionary, Nested_Boostrap_Index_single, Nested_Boostrap_Index_Special_single, Find_Controls
from GSTR_helpers import find_S
from find_S_helpers import find_S_gradient_descent_L1_loss, find_S_gradient_descent_se
import numpy as np
import copy

def find_significant_shrinkage(raw_df, input_sample_list1, input_sample_list2, cell_number_cutoff_list, input_control_gRNA_list, number_of_replicate, input_total_gRNA_number):
    # cell_number_cutoff is for treated (adjusted cutoff). This is given by me
    # first find S and untreated (basal) cutoff
    raw_treated_df = Generate_ref_input_df(raw_df,input_sample_list1,0) # this subset df based on cutoff
    raw_untreated_df = Generate_ref_input_df(raw_df,input_sample_list2,0)
    S_out_df = []
    for cutoff in cell_number_cutoff_list:
        shrinkage_result = bootstrap_shrinkage(raw_treated_df, raw_untreated_df, cutoff, input_control_gRNA_list, number_of_replicate, input_total_gRNA_number)
        S_out_df.append(pd.DataFrame(shrinkage_result))
    
    S_out_df_list = [pd.DataFrame(elem) for elem in S_out_df]
    final_S_df = pd.concat(S_out_df_list, ignore_index=False)
    final_S_summary_df = Generate_Final_Shrinkage_Summary_Dataframe(final_S_df, ['Shrinkage'], 'adjusted_cutoff')
    
    return final_S_df, final_S_summary_df

def bootstrap_shrinkage(raw_treated_df,raw_untreated_df,cell_number_cutoff,input_control_gRNA_list,number_of_replicate,input_total_gRNA_number):
    # Estimate Shrinkage (S) and find basal cutoff
    #S, basal_cutoff = find_S(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    #S, basal_cutoff = find_S_gradient_descent_L1_loss(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    S, basal_cutoff = find_S_gradient_descent_se(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    # Initialize list to store results
    results = [{'Shrinkage': S, 'Bootstrap_id': 'Real', 'basal_cutoff': basal_cutoff, 'adjusted_cutoff': cell_number_cutoff}]
    
    if number_of_replicate!=0:
        # treated
        Mouse_index_dic_1 = Generate_Index_Dictionary(raw_treated_df) # has tumor size info for each SampleID
        # untreated
        Mouse_index_dic_2 = Generate_Index_Dictionary(raw_untreated_df)
        for bootstrap_cycle in range(number_of_replicate):
            # resampling the treated mice
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = raw_treated_df.loc[x]
            # resampleing the untreated mice
            y = Nested_Boostrap_Index_Special_single(Mouse_index_dic_2,raw_untreated_df,input_total_gRNA_number)
            temp_bootstrap_df_2 = raw_untreated_df.loc[y]
            
            # re-estimate S
            #S, basal_cutoff = find_S(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            #S, basal_cutoff = find_S_gradient_descent_L1_loss(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            S, basal_cutoff = find_S_gradient_descent_se(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            # Append the result for this bootstrap cycle
            results.append({'Shrinkage': S, 'Bootstrap_id': f'B{bootstrap_cycle}',
                            'basal_cutoff': basal_cutoff, 'adjusted_cutoff': cell_number_cutoff})
    
    return results

def Generate_Final_Shrinkage_Summary_Dataframe(input_df,trait_of_interest,group_key):
    temp_summary = input_df[input_df['Bootstrap_id']!='Real'].groupby(group_key).apply(lambda x: pd.Series(np.percentile(x['Shrinkage'], [2.5, 50, 95, 97.5]), 
                                 index=['Shrinkage_2.5P', 'Shrinkage_50P', 'Shrinkage_95P', 'Shrinkage_97.5P'])).reset_index()
    temp_output_df = copy.deepcopy(input_df[input_df['Bootstrap_id'] =='Real'])
    temp_output_df = temp_output_df.merge(temp_summary, on = group_key)
    
    return(temp_output_df)

def main():
    parent_address = '/oak/stanford/scg/lab_mwinslow/Karen/Bootstrapping_analysis/ADJ4_LORSHP2_050824/Input_data'
    output_address = '/oak/stanford/scg/lab_mwinslow/Karen/Bootstrapping_analysis/ADJ4_LORSHP2_050824/Output_data/S'
    variant = 'v1'
    
    if variant == 'v1':
        #v1
        raw_df_address = parent_address + '/EA_drug_final_df_v1.csv'
        sample_to_exclude = ['ADJ4_188','ADJ4_194','ADJ4_200','ADJ4_206','ADJ4_210','ADJ4_213','ADJ4_218',
                        'ADJ4_219','ADJ4_227','ADJ4_228','ADJ4_230','ADJ4_237','ADJ4_240','ADJ4_256','ADJ4_262','ADJ4_279'] 
    elif variant == 'v3':
        # v3
        raw_df_address = parent_address + '/EA_drug_final_df_v3.csv'
        sample_to_exclude = ['ADJ4_188','ADJ4_206','ADJ4_213','ADJ4_218','ADJ4_219','ADJ4_227',
                        'ADJ4_228','ADJ4_230','ADJ4_237','ADJ4_240','ADJ4_256','ADJ4_262','ADJ4_194']

    
    raw_summary_df = pd.read_csv(raw_df_address)
    raw_summary_df= raw_summary_df[~raw_summary_df.Sample_ID.isin(sample_to_exclude)] # exclude the sample 
     
    temp_input = raw_summary_df[raw_summary_df['Identity']=='gRNA'] # consider only sgRNA but not spiekin
    sgRNA_number = len(temp_input[temp_input['Identity']=='gRNA']['gRNA'].unique())
     
    control_gRNA_list = Find_Controls(raw_summary_df,'Safe|Neo|NT')
    experiment_genotype = 'CE'
    control_genotype = 'CE'
    control_treatment = 'VEHICLE' # modify this
    exp_treatment = ['COMBO', 'LOR', 'T0', 'VEHICLE']
     
    interm_df_list = []
    sum_df_list = []
    for elem in exp_treatment:
        print(f'working on treatment: {elem}')
        cohort_1 = temp_input[(temp_input['Mouse_genotype'] == experiment_genotype)&(temp_input['Treatment'] == elem)]['Sample_ID'].unique()
        cohort_2 = temp_input[(temp_input['Mouse_genotype'] == control_genotype)&(temp_input['Treatment'] == control_treatment)]['Sample_ID'].unique()
        
        cutoff_list = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
        interm_df, sum_df = find_significant_shrinkage(temp_input, cohort_1, cohort_2, cutoff_list, control_gRNA_list, 100, sgRNA_number)
        interm_df['Treatment'] = elem
        sum_df['Treatment'] = elem
        interm_df_list.append(interm_df)
        sum_df_list.append(sum_df)
            
    final_interm_df = pd.concat(interm_df_list, ignore_index=True)
    final_sum_df = pd.concat(sum_df_list, ignore_index=True)
    
    final_interm_df.to_csv(f'{output_address}/S_intermediate_given_adjusted_cutoff_{variant}', index=None)
    final_sum_df.to_csv(f'{output_address}/S_summary_given_adjusted_cutoff_{variant}.csv', index=None)
    
if __name__ == "__main__":
    main() 