"""
    Generate two null distributions of ScoreRTN
- want to use null distribution of no GSTR to find p-value during bs
    1. null where all gRNAs have no GSTR and we treat all gRNAs as one gRNA. Then find one ScoreRTN and ScoreRGM
    2. null where gRNA i has no GSTR and we treated gRNA i seperately. This gives us ScoreRTN and ScoreRGM for each gRNA
- Chuan's method:
    - upsample the untreated mice with replacement to match mouse no. in the treated group
    - shrink the untreated tumors by S estimated from the treated group
    - apply the same cotoff as the treated group or the adjusted_cutoff given by me
- alternative mthod:
    - still upsample the untreated mice with replacement to match mouse no. in the treated group
    - use the basal_cutoff for the untreated group
        basal_cutoff = adjusted_cutoff / S
    - shrinkage is reflected by the cell number cutoff
"""
import pandas as pd
from GSTR_helpers import calculate_ScoreRTN

class NullScoreRTNCalculator():
    def __init__(self, treated_df, untreated_df, input_control_gRNA_list):
        self.treated_df = treated_df
        self.untreated_df = untreated_df
        self.input_control_gRNA_list = input_control_gRNA_list
    
    def calculate_ScoreRTN_group_non_inerts(self):
        """ for ScoreRTN null
        Args: 
            treated_df: usually synthetic treated df with TTN sum
            untreatd_df: usually observed untreated df with TTN sum
            input_control_gRNA_list: inert gRNA list
        Returns:
            df with TTN, RTN, ScoreRTN to save as intermediate df
        """
        TTN_inert_treatment = self.treated_df[self.treated_df['gRNA'].isin(self.input_control_gRNA_list)]['TTN'].sum()
        TTN_inert_untreated = self.untreated_df[self.untreated_df['gRNA'].isin(self.input_control_gRNA_list)]['TTN'].sum()
        
        # Calculate RTN for treatment and untreated
        treated_RTN_df = self.calculate_RTN_group_non_inerts(self.treated_df, TTN_inert_treatment)
        treated_RTN_df['TTN_inert'] = TTN_inert_treatment
        
        untreated_RTN_df = self.calculate_RTN_group_non_inerts(self.untreated_df, TTN_inert_untreated)
        untreated_RTN_df['TTN_inert'] = TTN_inert_untreated
        
        df_merged = treated_RTN_df.merge(untreated_RTN_df, on='gRNA', how='outer', suffixes=('_treated', '_untreated'))
        df_merged['ScoreRTN'] = np.log2(df_merged['RTN_treated'] / df_merged['RTN_untreated'])
            
        return df_merged
    
    def calculate_RTN_group_non_inerts(self, df, tumor_number_inert):
        TTN_sum_non_inert = df.loc[~df['gRNA'].isin(self.input_control_gRNA_list), 'TTN'].sum() # TTN sum for non-inerts
        RTN_non_inert = TTN_sum_non_inert / tumor_number_inert # RTN for non-inert gRNAs
        
        # Create a DataFrame for the combined non-inert row
        non_inert_row = pd.DataFrame({
            'gRNA': ['non_inerts'],
            'TTN': [TTN_sum_non_inert],
            'RTN': [RTN_non_inert]
        })
        
        # Filter out the inert gRNAs
        inert_rows = df[df['gRNA'].isin(self.input_control_gRNA_list)].copy(deep=True)
        # Calculate RTN for inert gRNAs
        inert_rows.loc[:, 'RTN'] = inert_rows['TTN'] / tumor_number_inert
        
        # Concatenate the non-inert row with the inert rows
        result_df = pd.concat([non_inert_row, inert_rows], ignore_index=True)
        return result_df
    
    def find_ScoreRTN_null(treated_df, untreated_df, input_control_gRNA_list):
        """ find null ScoreRTN
        Args:
            treated_df: treated_df with TTN sum
            untreated_df: untreated_df with TTN_sum
        """
        # first null group non-inerts
        ScoreRTN_df_group_non_inert = calculate_ScoreRTN_group_non_inerts(treated_df, untreated_df, input_control_gRNA_list)
        ScoreRTN_df_group_non_inert['Type'] = 'group'
        ScoreRTN_dict1 = df_to_dict(ScoreRTN_df_group_non_inert, 'gRNA', 'ScoreRTN')
        # second null treat each non-inert individually
        ScoreRTN_df = calculate_ScoreRTN(treated_df, untreated_df, input_control_gRNA_list)
        ScoreRTN_df['Type'] = 'indiv'
        ScoreRTN_dict2 = df_to_dict(ScoreRTN_df, 'gRNA', 'ScoreRTN')
        
        result_df = pd.concat([ScoreRTN_df_group_non_inert, ScoreRTN_df], ignore_index=True)
        return result_df, ScoreRTN_dict1, ScoreRTN_dict2
    
    def calculate_ScoreRTN_group_gRNAs(treated_df, untreated_df, input_control_gRNA_list):
        """ for ScoreRTN null
        Args: 
            treated_df: usually synthetic treated df with TTN sum
            untreatd_df: usually observed untreated df with TTN sum
            input_control_gRNA_list: inert gRNA list
        Returns:
            df with TTN, RTN, ScoreRTN to save as intermediate df
        """
        TTN_inert_treatment = treated_df[treated_df['gRNA'].isin(input_control_gRNA_list)]['TTN'].sum()
        TTN_inert_untreated = untreated_df[untreated_df['gRNA'].isin(input_control_gRNA_list)]['TTN'].sum()
        
        # Calculate RTN for treatment and untreated
        treated_RTN_df = calculate_RTN_group_gRANs(treated_df, TTN_inert_treatment, input_control_gRNA_list)
        treated_RTN_df['TTN_inert'] = TTN_inert_treatment
        
        untreated_RTN_df = calculate_RTN_group_gRANs(untreated_df, TTN_inert_untreated, input_control_gRNA_list)
        untreated_RTN_df['TTN_inert'] = TTN_inert_untreated
        
        df_merged = treated_RTN_df.merge(untreated_RTN_df, on='gRNA', how='outer', suffixes=('_treated', '_untreated'))
        df_merged['ScoreRTN'] = np.log2(df_merged['RTN_treated'] / df_merged['RTN_untreated'])
        
        return df_merged

    def calculate_RTN_group_gRANs(df, tumor_number_inert, input_control_gRNA_list):
        TTN_sum_non_inert = df.loc[~df['gRNA'].isin(input_control_gRNA_list), 'TTN'].sum() # TTN sum for non-inerts
        RTN_non_inert = TTN_sum_non_inert / tumor_number_inert # RTN for non-inert gRNAs
        
        RTN_inert = tumor_number_inert / tumor_number_inert # RTN for inert gRNAs
        # Create a DataFrame for the combined non-inert row
        result_df = pd.DataFrame({
            'gRNA': ['non_inerts', 'inerts'],
            'TTN': [TTN_sum_non_inert, tumor_number_inert],
            'RTN': [RTN_non_inert, RTN_inert]})
        return result_df
    
    def find_ScoreRTN_null_v2(treated_df, untreated_df, input_control_gRNA_list):
        """ find null ScoreRTN. this includes grouping of inerts and non-inerts
        Args:
            treated_df: treated_df with TTN sum
            untreated_df: untreated_df with TTN_sum
        """
        # first null group non-inerts
        ScoreRTN_df_group_gRNAs = calculate_ScoreRTN_group_gRNAs(treated_df, untreated_df, input_control_gRNA_list)
        # second null treat each non-inert individually
        ScoreRTN_df = calculate_ScoreRTN(treated_df, untreated_df, input_control_gRNA_list)

        result_df = pd.concat([ScoreRTN_df_group_gRNAs, ScoreRTN_df], ignore_index=True)
        ScoreRTN_dict = df_to_dict(result_df, 'gRNA', 'ScoreRTN')
        return result_df, ScoreRTN_dict