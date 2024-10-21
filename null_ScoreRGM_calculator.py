"""
    Generate two null distributions of ScoreRGM
- want to use null distribution of no GSTR to find p-value during bs
    1. null where all gRNAs have no GSTR and we treat all gRNAs as one gRNA. Then find one ScoreRGM for inerts and one for non-inerts
    2. null where gRNA i has no GSTR and we treated gRNA i seperately. This gives us ScoreRGM for each gRNA
- Chuan's method:
    - upsample the untreated mice with replacement to match mouse no. in the treated group
    - shrink the untreated tumors by S estimated from the treated group
    - apply the same cotoff as the treated group or the adjusted_cutoff given by me
- alternative mthod -> I use this:
    - still upsample the untreated mice with replacement to match mouse no. in the treated group
    - use the basal_cutoff for the untreated group
        basal_cutoff = adjusted_cutoff / S
    - shrinkage is reflected by the cell number cutoff
"""

import pandas as pd
from GSTR_helpers import prepare_RGM_data, compute_geo_mean, df_to_dict
from scipy.stats.mstats import gmean
import numpy as np

class NullScoreRGMCalculator():
    def __init__(self, treated_df, untreated_df, input_control_gRNA_list, ratio_dict):
        self.treated_df = treated_df # post cutoff
        self.untreated_df = untreated_df # post cutoff
        self.input_control_gRNA_list = input_control_gRNA_list
        self.ratio_dict = ratio_dict
        self.top_N_treated_df = None
        self.top_N_untreated_df = None
    
    def calculate_ScoreRGM_null(self):
        """
        Args:
            treated_df_cut: post cutoff
            untreated_df_cut: post cutoff
            ratio_dict: relative guide representation ratio to inert in each mouse

        Returns:
            df: concat gm_df grouped inert vs non-inert and indiv
        """
        # Step 1: Prepare Data
        self.top_N_treated_df, self.top_N_untreated_df = prepare_RGM_data(self.treated_df, self.untreated_df, self.ratio_dict, self.input_control_gRNA_list)
        
        # Step 2: Compute Geometric Means
        ## 1. group inerts and non-inerts respectively
        gm_df_grouped = self.calculate_ScoreRGM_null_group_gRNAs(self.top_N_treated_df, self.top_N_untreated_df)
        ## 2. treat each indiv
        gm_df_indiv = self.calculate_ScoreRGM_null_indiv_gRNAs(self.top_N_treated_df, self.top_N_untreated_df)

        # Step 3: Concat the indiv and grouped df
        final_df = pd.concat([gm_df_grouped, gm_df_indiv])
        gm_dict = df_to_dict(final_df,'gRNA', 'ScoreRGM')

        return final_df, gm_dict
    
    def compute_geo_mean_group_gRNAs(self, df, col_name):
        """ group df by inert vs non-inert gRNAs and calcualte geo mean respectively

        Args:
            df: any df with Cell_number and gRNA cols
            group_cols: gRNA, or [gRNA, Targeted_gene_name, Numbered_gene_name]
            col_name: normally geo_mean_treated or geo_mean_untreated
        Returns:
            gm_df with inert and non-inert gmean
        """
        gm_inert = gmean(df[df['gRNA'].isin(self.input_control_gRNA_list)]['Cell_number'])
        non_inert_df = df[~df['gRNA'].isin(self.input_control_gRNA_list)]
        gm_non_inert = gmean(non_inert_df['Cell_number'])
        gm_df = pd.DataFrame({'gRNA': ['non_inerts', 'inerts'],
                            col_name: [gm_non_inert, gm_inert],
                            f'{col_name}_inert': [gm_inert, gm_inert]})

        return gm_df

    def calculate_ScoreRGM_null_group_gRNAs(self, treated_df, untreated_df):
        """calculate grouped RGM for grouped inerts and non-inerts for treated and untreated

        Args:
            treated_df: top N tumors of each genotype in each mouse of the treated group
            untreated_df: top N tumors of each genotype in each mouse of the untreated group
        
        Returns:
            gm_treated_df: calculate RGM and ScoreRGM for grouped gRNAs
        """
        gm_treated_df = self.compute_geo_mean_group_gRNAs(treated_df, 'Geo_mean_treated')
        gm_untreated_df = self.compute_geo_mean_group_gRNAs(untreated_df, 'Geo_mean_untreated')
        gm_treated_df = gm_treated_df.merge(gm_untreated_df, on='gRNA')

        gm_treated_df['RGM_treated'] = gm_treated_df['Geo_mean_treated'] / gm_treated_df['Geo_mean_treated_inert']
        gm_treated_df['RGM_untreated'] = gm_treated_df['Geo_mean_untreated'] / gm_treated_df['Geo_mean_untreated_inert']

        gm_treated_df['ScoreRGM'] = np.log2(gm_treated_df['RGM_treated'] /  gm_treated_df['RGM_untreated'])

        return gm_treated_df

    def calculate_ScoreRGM_null_indiv_gRNAs(self, treated_df, untreated_df):
        """ calculate ScoreRGM for indiv gRNA
        Args:
            treated_df: top N tumors of each genotype in each mouse in the treated group
            untreated_df: top N tumors of each genotype in each mouse in the untreated group
        
        Returns:
            gm_treated_df: df with Geo_mean, RGM, ScoreRGM
        """
        gm_treated_df = compute_geo_mean(treated_df, self.input_control_gRNA_list, ['gRNA'], 'Geo_mean_treated')
        gm_untreated_df = compute_geo_mean(untreated_df, self.input_control_gRNA_list,  ['gRNA'], 'Geo_mean_untreated')

        gm_treated_df = gm_treated_df.merge(gm_untreated_df, on='gRNA')
        gm_treated_df['RGM_treated'] = gm_treated_df['Geo_mean_treated'] / gm_treated_df['Geo_mean_treated_inert']
        gm_treated_df['RGM_untreated'] = gm_treated_df['Geo_mean_untreated'] / gm_treated_df['Geo_mean_untreated_inert']
        gm_treated_df['ScoreRGM'] = np.log2(gm_treated_df['RGM_treated'] / gm_treated_df['RGM_untreated'])

        return gm_treated_df
