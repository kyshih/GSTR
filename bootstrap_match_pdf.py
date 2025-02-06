from shrinkage_analysis.bootstrap_matching_dist import BootstrapAnalyzer
import numpy as np
import pandas as pd

parent_address = '/labs/mwinslow/Karen/Bootstrapping_analysis/ADJ4_LORSHP2_050824'
df_v1_raw = pd.read_csv(parent_address + '/Input/EA_drug_final_df_v1.csv')
df_v3_raw = pd.read_csv(parent_address + '/Input/EA_drug_final_df_v3.csv')

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

df_v3_CE_combo = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'COMBO']
df_v3_CE_lor = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'LOR']
df_v3_CE_t0 = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'T0']
df_v3_CE_veh = df_v3_filtered_inert[df_v3_filtered_inert['Treatment'] == 'VEHICLE']

df_v1_CE_combo = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'COMBO']
df_v1_CE_lor = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'LOR']
df_v1_CE_t0 = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'T0']
df_v1_CE_veh = df_v1_filtered_inert[df_v1_filtered_inert['Treatment'] == 'VEHICLE']

def main():
    var_list = {['v3': df_v3_CE_combo, df_v3_CE_lor, df_v3_CE_t0, df_v3_CE_veh],
               'v1': [df_v1_CE_combo, df_v1_CE_lor, df_v1_CE_t0, df_v1_CE_veh]}
    all_results = []
    for var in var_list.keys()
        df_list = var_list[var]
        for treatment in df_list:
            df_temp = treatment
            for cutoff in [200, 500, 1000]:
                analyzer = BootstrapAnalyzer(df_temp, df_v1_CE_veh, cutoff)
                results_df = pd.DataFrame(analyzer.generate_bootstrap_samples(100, len(control_guides_v1)))
                all_results.append(results_df)
        