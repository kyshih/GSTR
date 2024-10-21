"""
I adapted Emily's GSTR calculation to find S and find RTN and ScoreRTN with bootstrapping.
This is different from scale_metrics_to_inert_base in following ways:
1. S is re-estimated using binary search in each bs cycle
2. p-value is calculated by comparing the bootstrapped statistics to the null distribution
"""
from bootstrapping_helpers import *
from metrics_helpers import *
from GSTR_helpers import *
import argparse
from find_S_helpers import find_S_gradient_descent_se_given_base_cutoff

def bootstrap_GSTR_and_shrinkage(raw_df,input_sample_list1,input_sample_list2,cell_number_cutoff,input_control_gRNA_list,number_of_replicate,input_total_gRNA_number):
    # cell_number_cutoff is for treated (adjusted cutoff). This is given by me
    # first find S and untreated (basal) cutoff
    raw_treated_df = Generate_ref_input_df(raw_df,input_sample_list1,0) # this subset df based on cutoff
    raw_untreated_df = Generate_ref_input_df(raw_df,input_sample_list2,0)
    
    # find R for RGM calculation
    R_dict = find_ratio_to_inert(raw_untreated_df, input_control_gRNA_list)
    print(f"observed ratio dict is {R_dict}")
    # estimate S
    S, adj_cutoff = find_S_gradient_descent_se_given_base_cutoff(raw_treated_df, raw_untreated_df, input_control_gRNA_list, cell_number_cutoff)
    # return tumor size metric of each bootstrap cycle
    # applying the cell no. cutoff
    # treated mouse   
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,adj_cutoff)
    # untreated mouse
    temp_ref_df2 = Generate_ref_input_df(raw_df,input_sample_list2,cell_number_cutoff)
    
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
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = raw_treated_df.loc[x]
            # resampleing the untreated mice
            y = Nested_Boostrap_Index_Special_single(Mouse_index_dic_2,raw_untreated_df,input_total_gRNA_number)
            temp_bootstrap_df_2 = raw_untreated_df.loc[y]
            
            # re-estimate R
            R_dict = find_ratio_to_inert(temp_bootstrap_df_2, input_control_gRNA_list)
            print(f"bs {bootstrap_cycle} ratio dict is {R_dict}")
            # re-estimate S
            S, adj_cutoff = find_S_gradient_descent_se_given_base_cutoff(temp_bootstrap_df_1, temp_bootstrap_df_2, input_control_gRNA_list, cell_number_cutoff)
            temp_bootstrap_ref_df1 = Generate_ref_input_df(temp_bootstrap_df_1, temp_bootstrap_df_1['Sample_ID'].unique(), adj_cutoff) # treated
            temp_bootstrap_ref_df2 = Generate_ref_input_df(temp_bootstrap_df_2, temp_bootstrap_df_2['Sample_ID'].unique(), cell_number_cutoff) # untreated
            
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
    
    test_final_df = bootstrap_GSTR_and_shrinkage(temp_input,cohort_1,cohort_2,cell_number_cutoff,control_gRNA_list,number_of_bootstrap,sgRNA_number)
    print(f"Bootstrapping steps have finished")
    
    if args.o2:
        test_final_df.to_csv(args.o2,index = False)
    else:
        print(f"No intermediate file output")
        
    if number_of_bootstrap!=0:
        # generate summary statistics
        temp_trait_list = ['ScoreRTN', 'ScoreRGM', 'Shrinkage']
        temp_trait_list = list(set(temp_trait_list))
        
        Final_summary_df = Generate_Final_Summary_Dataframe(test_final_df,temp_trait_list) # gRNA level
        Final_gene_summary_df = Generate_Gene_Level_Summary_Dataframe(test_final_df,temp_trait_list)
        
        Final_summary_df.to_csv(output_address+'.csv',index = False)
        Final_gene_summary_df.to_csv(output_address+'_gene_level.csv',index = False)
    else:
        test_final_df.to_csv(output_address+'.csv',index = False)
    print(f"All steps finished") 

if __name__ == "__main__":
    main() 