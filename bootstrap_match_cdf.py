import sys
import os
import copy
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from shrinkage_analysis.bootstrap_matching_dist import BootstrapAnalyzer
import numpy as np
import pandas as pd

def Generate_Final_Shrinkage_Summary_Dataframe(input_df,trait_of_interest,group_key):
    temp_summary = input_df[input_df['Bootstrap_id']!='Real'].groupby(group_key).apply(lambda x: pd.Series(np.nanpercentile(x['Shrinkage'], [2.5, 50, 95, 97.5]), 
                                 index=['Shrinkage_2.5P', 'Shrinkage_50P', 'Shrinkage_95P', 'Shrinkage_97.5P'])).reset_index()
    temp_output_df = copy.deepcopy(input_df[input_df['Bootstrap_id'] =='Real'])
    temp_output_df = temp_output_df.merge(temp_summary, on = group_key)
    
    return(temp_output_df)

parent_address = '/labs/mwinslow/Karen/Bootstrapping_analysis/ADJ4_LORSHP2_050824'
df_v1_raw = pd.read_csv(parent_address + '/Input/EA_drug_final_df_v1.csv')
df_v1_raw['Variant'] = 'v1'
df_v3_raw = pd.read_csv(parent_address + '/Input/EA_drug_final_df_v3.csv')
df_v3_raw['Variant'] = 'v3'

discarded_samples_v1_f = f'{parent_address}/Input/Discarded_sample_list_for_EA_drug_v1.txt'
discarded_samples_v3_f = f'{parent_address}/Input/Discarded_sample_list_for_EA_drug_v3.txt'

with open (discarded_samples_v1_f, 'r') as f:
    lines = f.readlines()
f.close()
discarded_samples_v1 = [l.strip() for l in lines]

with open (discarded_samples_v3_f, 'r') as f:
    discarded_samples_v3 = [line.strip() for line in f]
f.close()

# just look at CE mice
df_v1_filtered = df_v1_raw[df_v1_raw['Mouse_genotype']=='CE']
df_v3_filtered = df_v3_raw[df_v3_raw['Mouse_genotype']=='CE']

# remove the bad mice
df_v1_filtered = df_v1_filtered[~df_v1_filtered['Sample_ID'].isin(discarded_samples_v1)]
df_v3_filtered = df_v3_filtered[~df_v3_filtered['Sample_ID'].isin(discarded_samples_v3)]

# focus on inert tumors to quantify treatment response
df_v1_filtered_inert = df_v1_filtered[df_v1_filtered['Numbered_gene_name'].str.contains('NT|Neo|Safe')]
df_v3_filtered_inert = df_v3_filtered[df_v3_filtered['Numbered_gene_name'].str.contains('NT|Neo|Safe')]

# inert guides
control_guides_v1 = df_v1_filtered_inert['gRNA'].unique()
control_guides_v3 = df_v3_filtered_inert['gRNA'].unique()

df_v3_CE_combo = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'COMBO']
df_v3_CE_lor = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'LOR']
df_v3_CE_t0 = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'T0']
df_v3_CE_veh = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'VEHICLE']

df_v1_CE_combo = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'COMBO']
df_v1_CE_lor = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'LOR']
df_v1_CE_t0 = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'T0']
df_v1_CE_veh = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'VEHICLE']

def main():
    var_list = {'v3': [df_v3_CE_combo, df_v3_CE_lor, df_v3_CE_t0, df_v3_CE_veh],
               'v1': [df_v1_CE_combo, df_v1_CE_lor, df_v1_CE_t0, df_v1_CE_veh]}
    treatment_names = ['Combo', 'Lor', 'T0', 'Vehicle']
    inert_guides = {'v3': control_guides_v3, 'v1': control_guides_v1}
    veh_df = {'v3': df_v3_CE_veh, 'v1': df_v1_CE_veh}
    all_results = []
    for var in var_list.keys():
        control_guides = inert_guides[var]
        untreated_df = veh_df[var]
        df_list = var_list[var]
        
        for i, treatment_df in enumerate(df_list):
            treatment_name = treatment_names[i]  # Get corresponding treatment name
            
            #for cutoff in [500, 1000, 2000]:
            for cutoff in [200]:
                print(f'working on cutoff {cutoff}')
                analyzer = BootstrapAnalyzer(treatment_df, untreated_df, cutoff)
                results_df = pd.DataFrame(analyzer.generate_bootstrap_samples(1000, len(control_guides)))
                results_df['Treatment'] = treatment_name # Assign correct treatment name
                results_df['Variant'] = var  # Add the variant (v1 or v3)
                all_results.append(results_df)
                
    final_results_df = pd.concat(all_results, ignore_index=True)

    final_resutls_summary_df = Generate_Final_Shrinkage_Summary_Dataframe(final_results_df,'Shrinkage',['Variant', 'Treatment', 'cutoff_treated'])
    
    final_results_df.to_csv(f'{parent_address}/Output/S/match_size_dist/ks_cdf_diff/S_intermediate', index=None)
    final_resutls_summary_df.to_csv(f'{parent_address}/Output/S/match_size_dist/ks_cdf_diff/S.csv', index=None)

if __name__ == "__main__":
    main()