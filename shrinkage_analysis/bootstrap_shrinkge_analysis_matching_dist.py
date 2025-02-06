import pandas as pd
import argparse
import numpy as np
import sys
sys.path.append('.')
from shrinkage_analysis.data_loader import DataLoader
from shrinkage_analysis.bootstrap_matching_dist import BootstrapAnalyzer
from shrinkage_analysis.data_utils import filter_df
# TODO: this still gets cannot find module error
""" This module matches the distribution of treated and untreated. It processes v1 and v3 together since v1 and v3 are in the same mice
"""

def find_shrinkage_given_adjusted_cutoff(raw_df, sample_list1, sample_list2, cutoff_list, control_gRNA_list, n_replicates, total_gRNA_count,
                                                 objective=None, method=None, find_S_version=None):
    # subset the treated mice and inert tumors
    treated_df = filter_df(raw_df, 'Sample_ID', sample_list1, exclude=False)
    treated_df = filter_df(treated_df, 'gRNA', control_gRNA_list, exclude=False)
    
    untreated_df = filter_df(raw_df, 'Sample_ID', sample_list2, exclude=False)
    untreated_df = filter_df(untreated_df, 'gRNA', control_gRNA_list, exclude=False)
    
    all_results = []
    for cutoff in cutoff_list:
        print(f'Processing cutoff {cutoff}...')
        analyzer = BootstrapAnalyzer(treated_df, untreated_df, cutoff)
        results_df = pd.DataFrame(analyzer.generate_bootstrap_samples(n_replicates, total_gRNA_count))
        all_results.append(results_df)
        
    final_results_df = pd.concat(all_results, ignore_index=True)
    #final_summary_df = generate_final_shrinkage_summary(final_results_df, 'Shrinkage', 'cutoff_treated')
    final_summary_df = None
    return final_results_df, final_summary_df

def generate_final_shrinkage_summary(input_df, trait_key, group_key):
    summary = input_df[input_df['Bootstrap_id'] != 'Real'].groupby(group_key).apply(
        lambda x: pd.Series(
            np.nanpercentile(x[trait_key], [2.5, 50, 95, 97.5]), 
            index=[f'{trait_key}_2.5P', f'{trait_key}_50P', f'{trait_key}_95P', f'{trait_key}_97.5P']
        )
    ).reset_index()
    
    real_data = input_df[input_df['Bootstrap_id'] == 'Real']
    return real_data.merge(summary, on=group_key)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Bootstrap shrinkage analysis across cutoffs')
    parser.add_argument('-var', '--variant', choices=['v1', 'v3'], required=True, help='EA variant to run')
    parser.add_argument('-l', '--cutoffs', nargs='+', required=True, help='List of base cutoffs')
    parser.add_argument('-b', '--bs_rep', required=True, help='number of bootstrap replicates')
    parser.add_argument('-match', '--match_criteria', required=False, help='matching criteria fo find shrinkage')
    parser.add_argument('-obj', '--objective', required=False, help='objective function to minimize')
    parser.add_argument('-m', '--method', required=False, help='optimization method')
    parser.add_argument('-o', '--output', required=False, help='output address')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    return parser.parse_args()

def main():
    args = parse_arguments()
    cutoff_list = [int(cutoff) for cutoff in args.cutoffs]
    variant = args.variant
    no_rep = args.bs_rep
    # if args.match_criteria == "median_inert_TTN":
    #     imported_find_S = find_S_matching_median_TTN
    # elif args.match_criteria == "inert_size_dist":
    #     imported_find_S = find_S_matching_distribution
    # else:
    #     raise ValueError(f"Invalid match_criteria: {args.match_criteria}")
    
    if args.verbose:
        print(f"Verbose mode is ON")
        # print(f"Input file: {args.input_file}")
        # print(f"Output directory: {args.output_dir}")
        print(f"Given adjusted cutoff in treated group, shrink vehicle group to match treated group.")
        print(f"Variant: {variant}")
        print(f"Start processing tumors and calculate S...")
        # print(f"matching criteria: {args.match_criteria}")
        # print(f"objective is {args.objective}")
        # print(f"Optimization method: {args.method}")
        
    # Define file paths based on variant
    parent_address = '/oak/stanford/scg/lab_mwinslow/Karen/Bootstrapping_analysis/ADJ4_LORSHP2_050824/Input'
    output_address = '/oak/stanford/scg/lab_mwinslow/Karen/Bootstrapping_analysis/ADJ4_LORSHP2_050824/Output/S/match_size_dist/ks_cdf_diff'
    
    variant_files = {
        'v1': {
            'data': f'{parent_address}/EA_drug_final_df_v1.csv',
            #'discard': f'{parent_address}/Discarded_sample_list_for_EA_drug_v1.txt'
            'discard': f'{parent_address}/Discarded_sample_list_for_EA_drug.txt'
        },
        'v3': {
            'data': f'{parent_address}/EA_drug_final_df_v3.csv',
            #'discard': f'{parent_address}/Discarded_sample_list_for_EA_drug_v3.txt'
            'discard': f'{parent_address}/Discarded_sample_list_for_EA_drug.txt'
        }
    }
    
    # Load data using DataLoader
    # data_loader = DataLoader(
    #     raw_df_path=variant_files[variant]['data'],
    #     discard_samples_path=variant_files[variant]['discard']
    # )
    data_loader_v1 = DataLoader(
        raw_df_path=variant_files['v1']['data'],
        discard_samples_path=variant_files['v1']['discard']
    )
    data_loader_v3 = DataLoader(
        raw_df_path=variant_files['v3']['data'],
        discard_samples_path=variant_files['v3']['discard']
    )
    
    # raw_df_no_bad_samples = data_loader.load_and_exclude_samples()
    # filtered_raw_df = data_loader.get_gRNA_data(raw_df_no_bad_samples)
    # control_gRNA_list = data_loader.find_control_gRNAs(filtered_raw_df)
    # total_gRNAs = filtered_raw_df['gRNA'].nunique()
    
    raw_df_no_bad_samples_v1 = data_loader_v1.load_and_exclude_samples()
    filtered_raw_df_v1 = data_loader_v1.get_gRNA_data(raw_df_no_bad_samples_v1)
    filtered_raw_df_v1['Variant'] = 'v1'
    
    raw_df_no_bad_samples_v3 = data_loader_v3.load_and_exclude_samples()
    filtered_raw_df_v3 = data_loader_v3.get_gRNA_data(raw_df_no_bad_samples_v3)
    filtered_raw_df_v3['Variant'] = 'v3'
    
    control_gRNA_list = data_loader_v3.find_control_gRNAs(filtered_raw_df_v3)
    total_gRNAs = filtered_raw_df_v3['gRNA'].nunique()
    
    filtered_raw_df = pd.concat([filtered_raw_df_v1, filtered_raw_df_v3], ignore_index=True)
    
    # Filtered data based on genotype and treatment
    exp_genotype, control_genotype = 'CE', 'CE'
    control_treatment = 'VEHICLE'
    treatments = ['COMBO', 'LOR', 'T0', 'VEHICLE']
    
    all_intermediates = []
    all_summaries = []

    for treatment in treatments:
        print(f"Processing treatment: {treatment}")
        treated_samples = filtered_raw_df.query("Mouse_genotype == @exp_genotype & Treatment == @treatment")['Sample_ID'].unique()
        control_samples = filtered_raw_df.query("Mouse_genotype == @control_genotype & Treatment == @control_treatment")['Sample_ID'].unique()
        
        intermediate_df, summary_df = find_shrinkage_given_adjusted_cutoff(
            filtered_raw_df, treated_samples, control_samples, cutoff_list, control_gRNA_list, n_replicates=no_rep, total_gRNA_count=total_gRNAs)
            #args.objective, args.method, imported_find_S)
        intermediate_df['Treatment'] = treatment
        summary_df['Treatment'] = treatment
        all_intermediates.append(intermediate_df)
        all_summaries.append(summary_df)
    
    final_intermediate_df = pd.concat(all_intermediates, ignore_index=True)
    final_summary_df = pd.concat(all_summaries, ignore_index=True)
    
    # Save results
    final_intermediate_df.to_csv(f'{output_address}/S_intermediate_given_adjusted_cutoff_{variant}', index=False)
    #final_summary_df.to_csv(f'{output_address}/S_summary_given_adjusted_cutoff_{variant}.csv', index=False)
    
    if args.verbose:
        print("Processing complete. Results saved.")

if __name__ == "__main__":
    main()
    