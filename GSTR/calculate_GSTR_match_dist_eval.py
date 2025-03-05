"""
I adapted Emily's GSTR calculation to find S and find RTN and ScoreRTN with bootstrapping.
This is different from scale_metrics_to_inert_base in following ways:
1. S is re-estimated by using ks stats (min vertical distance between CDFs)
2. p-value is calculated by comparing the bootstrapped statistics to the null distribution
3. cell number cutoff is given to treated group
"""
import sys
import os
# Add the parent directory (Python/) to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils.bootstrapping_helpers import *
#from metrics_helpers import *
from GSTR.GSTR_helpers import *
import argparse
#from find_S_helpers import find_S_gradient_descent_se
from shrinkage_analysis.find_S_cdf import find_optimal_S_grid

def bootstrap_GSTR_and_shrinkage(raw_df,input_sample_list1,input_sample_list2,cell_number_cutoff,input_control_gRNA_list,number_of_replicate,input_total_gRNA_number):
    # cell_number_cutoff is for treated (adjusted cutoff). This is given by me
    # first find S and untreated (basal) cutoff
    raw_treated_df = Generate_ref_input_df(raw_df,input_sample_list1,0) # this subset df based on cutoff
    raw_untreated_df = Generate_ref_input_df(raw_df,input_sample_list2,0)
    
    # find R for RGM calculation
    R_dict = find_ratio_to_inert(raw_untreated_df, input_control_gRNA_list)
    print(f"observed ratio dict is {R_dict}")
    # estimate S
    #S, basal_cutoff = find_S(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    #S, basal_cutoff = find_S_gradient_descent_se(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    raw_untreated_df_inert = raw_untreated_df[raw_untreated_df['gRNA'].isin(input_control_gRNA_list)]
    raw_treated_df_inert = raw_treated_df[raw_treated_df['gRNA'].isin(input_control_gRNA_list)]
    # S, min_ks, ks_distances, S_values = find_optimal_S_grid(treated_df=raw_treated_df_inert, untreated_df=raw_untreated_df_inert,
    #                                               cutoff_tr=cell_number_cutoff)
    S = 0.7
    basal_cutoff = cell_number_cutoff/S
    # return tumor size metric of each bootstrap cycle
    # applying the cell no. cutoff
    # treated mouse   
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,cell_number_cutoff)
    # untreated mouse
    temp_ref_df2 = Generate_ref_input_df(raw_df,input_sample_list2,basal_cutoff)
    
    # dfs passed are post cutoff
    temp_final_df_observed = calculate_GSTR_metrics(temp_ref_df1, temp_ref_df2, input_control_gRNA_list, R_dict)
    temp_final_df_observed['Shrinkage'] = S
    temp_final_df_observed['Bootstrap_id'] = 'Real' # a new columns named Bootstrap_id. The values of the column are 'Real'
    
    # for concatenating dfs later
    temp_out_df = [temp_final_df_observed]
    
    if number_of_replicate!=0:
        # treated
        Mouse_index_dic_1 = Generate_Index_Dictionary(raw_treated_df) # has tumor size info for each SampleID
        # untreated
        Mouse_index_dic_2 = Generate_Index_Dictionary(raw_untreated_df)
        for bootstrap_cycle in range(number_of_replicate):
            # resampling the treated mice
            x = Nested_Bootstrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = raw_treated_df.loc[x]
            # resampleing the untreated mice
            y = Nested_Bootstrap_Index_Special_single(Mouse_index_dic_2,raw_untreated_df,input_total_gRNA_number)
            temp_bootstrap_df_2 = raw_untreated_df.loc[y]
            
            # re-estimate R
            R_dict = find_ratio_to_inert(temp_bootstrap_df_2, input_control_gRNA_list)
            print(f"bs {bootstrap_cycle} ratio dict is {R_dict}")
            # re-estimate S
            #S, basal_cutoff = find_S(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            #S, basal_cutoff = find_S_gradient_descent_se(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            temp_bootstrap_df_1_inert = temp_bootstrap_df_1[temp_bootstrap_df_1['gRNA'].isin(input_control_gRNA_list)]
            temp_bootstrap_df_2_inert = temp_bootstrap_df_2[temp_bootstrap_df_2['gRNA'].isin(input_control_gRNA_list)]
            # S, min_ks, ks_distances, S_values = find_optimal_S_grid(treated_df=temp_bootstrap_df_1_inert, untreated_df=temp_bootstrap_df_2_inert,
            #                                       cutoff_tr=cell_number_cutoff)
            S = 0.7
            basal_cutoff = cell_number_cutoff/S
            temp_bootstrap_ref_df1 = Generate_ref_input_df(temp_bootstrap_df_1, temp_bootstrap_df_1['Sample_ID'].unique(), cell_number_cutoff)
            temp_bootstrap_ref_df2 = Generate_ref_input_df(temp_bootstrap_df_2, temp_bootstrap_df_2['Sample_ID'].unique(), basal_cutoff)
            
            temp_metric_df = calculate_GSTR_metrics(temp_bootstrap_ref_df1, temp_bootstrap_ref_df2, input_control_gRNA_list, R_dict)
            temp_metric_df['Shrinkage'] = S
            temp_metric_df['Bootstrap_id'] = 'B'+str(bootstrap_cycle)
            
            temp_out_df.append(temp_metric_df)
    
    temp_out_df = pd.concat(temp_out_df, ignore_index=True)
    
    return temp_out_df

def calculate_GSTR_metrics(treated_df,untreated_df,input_control_gRNA_list, ratio_dict):
    # treated and untreated dfs are post cutoff
    treated_sum_df = treated_df.groupby(['gRNA']).Clonal_barcode.count().reset_index(name='TTN')
    untreated_sum_df = untreated_df.groupby(['gRNA']).Clonal_barcode.count().reset_index(name='TTN')
    ScoreRTN_df = calculate_ScoreRTN(treated_sum_df, untreated_sum_df, input_control_gRNA_list)
    ScoreRGM_df = calculate_ScoreRGM(treated_df, untreated_df, ratio_dict, input_control_gRNA_list)
    df_merged = ScoreRGM_df.merge(ScoreRTN_df, on='gRNA', how='outer')
    return df_merged

def calculate_GSTR_metrics_var(bootstrap_df):
    # Step 1: Compute variance (σ²) of ScoreRTN and ScoreRGM for each gene
    gene_variance = bootstrap_df.groupby(["Numbered_gene_name", "gRNA"], as_index=False)[["ScoreRTN", "ScoreRGM"]].var().rename(
    columns={"ScoreRTN": "Var_ScoreRTN", "ScoreRGM": "Var_ScoreRGM"})  
    print(gene_variance)
    return gene_variance

def calculate_G(bootstrap_summary_df):
    bootstrap_summary_df['ScoreGSTR'] = (
        (bootstrap_summary_df["ScoreRTN"] / bootstrap_summary_df["Var_ScoreRTN"]) +
        (bootstrap_summary_df["ScoreRGM"] / bootstrap_summary_df["Var_ScoreRGM"])
    ) / (
        (1 / bootstrap_summary_df["Var_ScoreRTN"]) + (1 / bootstrap_summary_df["Var_ScoreRGM"])
    )
    
    bootstrap_summary_df['G_hat'] = 2 ** bootstrap_summary_df["ScoreGSTR"] - 1

def calculate_GSTR_metrics_var_gene_level(bootstrap_df):
    temp_trait_list = ['ScoreRTN', 'ScoreRGM']
    temp_df = bootstrap_df[bootstrap_df['Bootstrap_id']!='Real'].groupby([
            'Targeted_gene_name','Bootstrap_id'],as_index = False).apply(Cal_Combined_Gene_Effect_v2,(temp_trait_list))
    gene_variance = bootstrap_df.groupby(["Targeted_gene_name"], as_index=False)[["ScoreRTN", "ScoreRGM"]].var().rename(
    columns={"ScoreRTN": "Var_ScoreRTN", "ScoreRGM": "Var_ScoreRGM"})  
    print(gene_variance)
    return gene_variance

def main():
    parser = argparse.ArgumentParser(description='A function to do resampling of mice')
    parser.add_argument("--a0", required=True, help="Address of processed data of Ultra-seq, can take multiple input")
    parser.add_argument("--a1", required=False, help="Sample to exclude list address")
    parser.add_argument("--a2", required=True, help="Basel cell number cutoff")
    parser.add_argument("--a3", required=True, help="Number of bootstrapping repeat")
    parser.add_argument("--a4", required=True, help="This is the experiment genotype")
    parser.add_argument("--a5", required=True, help="This is the control genotype")
    parser.add_argument("--a6", required=True, help='This is the experimental treatment')
    parser.add_argument("--a7", required=False, help='This is the control treatment')
    parser.add_argument("--o1", required=True, help="This the output address for summary data")
    parser.add_argument("--o2", required=False, help="This the output address for intermediate data") # df with tumor size metric of each bootstrap cycle
    parser.add_argument('--l2', nargs='+', required=False, help="A list of sgRNA sequence to exclude")
    
    # data input
    args = parser.parse_args()
    
    raw_df_input_address  = args.a0
    output_address = args.o1
    
    cell_number_cutoff = int(args.a2)
    print(f'cell number cutoff is {cell_number_cutoff}')

    experiment_genotype = args.a4
    print(f'The experiment mouse genotype is {experiment_genotype}')
    control_genotype = args.a5
    print(f'The control mouse genotype is {control_genotype}')
    exp_treatment = args.a6
    print(f'The experimental treatment in {experiment_genotype} mice is {exp_treatment}')

    number_of_bootstrap = int(args.a3)

    if args.l2 is None: # gRNA to exclude
        sgRNA_to_exclude = []
        print(f"No sgRNA is excluded from the analysis")
    else:
        sgRNA_to_exclude = args.l2
        print(f"sgRNAs excluded from the analysis:{sgRNA_to_exclude}")
        
    if args.a1 is None: # sample to exclude
        sample_to_exclude = []
        print(f"No sample is excluded from the analysis")
    else:
        sample_discarded_list_address = args.a1
        with open(sample_discarded_list_address, 'r') as f:
            sample_to_exclude = [line.rstrip('\n') for line in f]
        print(f"Samples excluded from the analysis:{sample_to_exclude}")
    
    if args.a7:
        control_treatment = args.a7
        print(f'The control treatment is {control_treatment}')
    else:
        control_treatment = None

    #raw_summary_df = pd.read_csv(raw_df_input_address) # read input data
    raw_summary_df = pd.read_csv(raw_df_input_address) # read input data
    control_gRNA_list = Find_Controls(raw_summary_df,'Safe|Neo|NT') # find the inert gRNA based on their targeted gene name
    
    # Generate bootstrapped df 
    raw_summary_df = raw_summary_df[~raw_summary_df['gRNA'].isin(sgRNA_to_exclude)] # exclude gRNA
    raw_summary_df= raw_summary_df[~raw_summary_df.Sample_ID.isin(sample_to_exclude)] # exclude the sample 
    temp_input = raw_summary_df[raw_summary_df['Identity']=='gRNA'] # consider only sgRNA but not spiekin
    
    # ----------------------------------------------------------------------------------------
    # only look at inert and pten for model eval
    gene_of_interest = ['Keap1', 'Kmt2c', 'Nf1', 'Pten', 'STAT3', 'Trp53', 'Safe33', 'Safe147', 'Safe152', 'Safe157', 'Neo1', 'Neo2', 'Neo3',
                        'NT1', 'NT2', 'NT3']
    temp_input = temp_input[temp_input['Targeted_gene_name'].isin(gene_of_interest)] 
    
    sgRNA_number = len(temp_input[temp_input['Identity']=='gRNA']['gRNA'].unique())
    # I want to generate two name list of mice, one for experimental group and another one for control group.
    # experimental mouse group
    cohort_1 = temp_input[(temp_input['Mouse_genotype'] == experiment_genotype)&(temp_input['Treatment'] == exp_treatment)]['Sample_ID'].unique()
    print(f"There are {len(cohort_1):d} experiment mice")
    # control mouse group
    if control_treatment:
        cohort_2 = temp_input[(temp_input['Mouse_genotype'] == control_genotype)&(temp_input['Treatment'] == control_treatment)]['Sample_ID'].unique()
    else:
        cohort_2 = temp_input[(temp_input['Mouse_genotype'] == control_genotype)]['Sample_ID'].unique()
    print(f"There are {len(cohort_2):d} control mice")
    
    # ----------------------------------------------------------------------------------------
    # repeat many times
    test_final_df = []
    for i in range(1):
        temp_test_final_df = bootstrap_GSTR_and_shrinkage(temp_input,cohort_1,cohort_2,cell_number_cutoff,control_gRNA_list,number_of_bootstrap,sgRNA_number)
        temp_test_final_df['Bootstrap_outer_id'] = 'B' + str(i)
        test_final_df.append(temp_test_final_df)
    test_final_df = pd.concat(test_final_df)
    print(f"Bootstrapping steps have finished")
    
    if args.o2:
        test_final_df.to_csv(args.o2,index = False)
    else:
        print(f"No intermediate file output")
        
    if number_of_bootstrap!=0:
        # generate summary statistics
        temp_trait_list = ['ScoreRTN', 'ScoreRGM', 'Shrinkage']
        temp_trait_list = list(set(temp_trait_list))
        
        Final_summary_df = []
        Final_gene_summary_df = []
        
        for i, group in test_final_df.groupby(['Bootstrap_outer_id']):
            temp_summary_df = Generate_Final_Summary_Dataframe(group,temp_trait_list) # gRNA level
            temp_GSTR_metrics_var_df = calculate_GSTR_metrics_var(group)
            # calculate_G(temp_GSTR_metrics_var_df)
            temp_summary_df = temp_summary_df.merge(temp_GSTR_metrics_var_df, on=['Numbered_gene_name', 'gRNA'])
            calculate_G(temp_summary_df)
            
            temp_gene_summary_df = Generate_Gene_Level_Summary_Dataframe(group,temp_trait_list)
            temp_GSTR_metrics_var_gene_df = calculate_GSTR_metrics_var_gene_level(group)
            temp_gene_summary_df = temp_gene_summary_df.merge(temp_GSTR_metrics_var_gene_df, on=['Targeted_gene_name'])
            calculate_G(temp_gene_summary_df)
            
            temp_summary_df['Bootstrap_outer_id'] = i[0]
            temp_gene_summary_df['Bootstrap_outer_id'] = i[0]
            Final_summary_df.append(temp_summary_df)
            Final_gene_summary_df.append(temp_gene_summary_df)
        
        Final_summary_df = pd.concat(Final_summary_df)
        Final_gene_summary_df = pd.concat(Final_gene_summary_df)
        
        Final_summary_df.to_csv(output_address+'.csv',index = False)
        Final_gene_summary_df.to_csv(output_address+'_gene_level.csv',index = False)
    else:
        test_final_df.to_csv(output_address+'.csv',index = False)
    print(f"All steps finished") 

if __name__ == "__main__":
    main() 