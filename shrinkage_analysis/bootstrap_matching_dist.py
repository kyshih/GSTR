import pandas as pd
import numpy as np
import importlib
from typing import List, Dict, Tuple
from scipy.optimize import minimize_scalar
from utils.bootstrapping_helpers import (
    Generate_ref_input_df,
    Generate_Index_Dictionary,
    Nested_Bootstrap_Index_single,
    Nested_Bootstrap_Index_Special_single
)
from shrinkage_analysis.find_S_cdf import find_optimal_S_grid
from shrinkage_analysis.find_S_max_pdf_auc import find_optimal_S_max_auc

class BootstrapAnalyzer:
    def __init__(self, raw_treated_df: pd.DataFrame, raw_untreated_df: pd.DataFrame,
                 cell_number_cutoff: int):
        self.treated_df = raw_treated_df
        self.untreated_df = raw_untreated_df
        self.cell_number_cutoff = cell_number_cutoff
        self.v1_treated = raw_treated_df[raw_treated_df['Variant']=='v1']
        self.v3_treated = raw_treated_df[raw_treated_df['Variant']=='v3']
        self.v1_untreated = raw_untreated_df[raw_untreated_df['Variant']=='v1']
        self.v3_untreated = raw_untreated_df[raw_untreated_df['Variant']=='v3']

    def generate_bootstrap_samples(self, n_replicates: int, total_gRNA_number: int) -> List[Dict]:
        results = self._calculate_shrinkage()
        
        if n_replicates == 0:
            return results

        treated_indices = Generate_Index_Dictionary(self.treated_df)
        untreated_indices = Generate_Index_Dictionary(self.untreated_df)

        for i in range(n_replicates):
            bootstrap_results = self._perform_single_bootstrap(
                treated_indices, untreated_indices, total_gRNA_number, i
            )
            results.append(bootstrap_results)

        return results

    def _calculate_shrinkage(self) -> List[Dict]:   
        best_S, min_ks, ks_distances, S_values = find_optimal_S_grid(treated_df=self.treated_df, untreated_df=self.untreated_df,
                                                  cutoff_tr=self.cell_number_cutoff)
        return [{
            'Shrinkage': best_S,
            'Bootstrap_id': 'Real',
            'cutoff_treated': self.cell_number_cutoff,
            'distance': min_ks
        }]
    
    # def _calculate_shrinkage(self) -> List[Dict]:   
        # best_S_v1, min_ks_v1, ks_distances_v1, S_values_v1 = find_optimal_S_grid(treated_df=self.v1_treated, untreated_df=self.v1_untreated,
        #                                           cutoff_tr=self.cell_number_cutoff)
        # best_S_v3, min_ks_v3, ks_distances_v3, S_values_v3 = find_optimal_S_grid(treated_df=self.v3_treated, untreated_df=self.v3_untreated,
        #                                           cutoff_tr=self.cell_number_cutoff)
        # return [{
        #     'Bootstrap_id': 'Real',
        #     'cutoff_treated': self.cell_number_cutoff,
        #     'Shrinkage_v1': best_S_v1,
        #     'distance_v1': min_ks_v1,
        #     'Shrinkage_v3': best_S_v3,
        #     'distance_v3': min_ks_v3
        # }]
        # best_S_v1, max_overlap_v1, overlaps_v1 = find_optimal_S_max_auc(treated_df=self.v1_treated, vehicle_df=self.v1_untreated,
        #                                                         cutoff_tr=self.cell_number_cutoff)
        # best_S_v3, max_overlap_v3, overlaps_v3 = find_optimal_S_max_auc(treated_df=self.v3_treated, vehicle_df=self.v3_untreated,
        #                                                         cutoff_tr=self.cell_number_cutoff)
        # return [{
        #     'Bootstrap_id': 'Real',
        #     'cutoff_treated': self.cell_number_cutoff,
        #     'Shrinkage_v1': best_S_v1,
        #     'overlap_v1': max_overlap_v1,
        #     'Shrinkage_v3': best_S_v3,
        #     'overlap_v3': max_overlap_v3
        # }]
                
    def _perform_single_bootstrap(self, treated_indices: Dict, untreated_indices: Dict, 
                                  total_gRNA_number: int, bs_cycle='') -> Dict:
        # ks cdf
        treated_sample = Nested_Bootstrap_Index_single(treated_indices)
        untreated_sample = Nested_Bootstrap_Index_Special_single(untreated_indices, self.untreated_df, total_gRNA_number)
        
        temp_treated = self.treated_df.loc[treated_sample]
        temp_untreated = self.untreated_df.loc[untreated_sample]
        
        best_S, min_ks, ks_distances, S_values = find_optimal_S_grid(treated_df=temp_treated, untreated_df=temp_untreated,
                                                  cutoff_tr=self.cell_number_cutoff)
        return {
            'Shrinkage': best_S,
            'Bootstrap_id': f'B{bs_cycle}',
            'cutoff_treated': self.cell_number_cutoff,
            'distance': min_ks
        }
    
    # def _perform_single_bootstrap(self, treated_indices: Dict, untreated_indices: Dict, 
    #                               total_gRNA_number: int, bs_cycle='') -> Dict:
    #     # ks cdf
    #     treated_sample = Nested_Bootstrap_Index_single(treated_indices)
    #     untreated_sample = Nested_Bootstrap_Index_Special_single(untreated_indices, self.untreated_df, total_gRNA_number)
        
    #     temp_treated = self.treated_df.loc[treated_sample]
    #     temp_untreated = self.untreated_df.loc[untreated_sample]
        
    #     temp_treated_v1 = temp_treated[temp_treated['Variant']=='v1']
    #     temp_treated_v3 = temp_treated[temp_treated['Variant']=='v3']
    #     temp_untreated_v1 = temp_untreated[temp_untreated['Variant']=='v1']
    #     temp_untreated_v3 = temp_untreated[temp_untreated['Variant']=='v3']
        
        # best_S_v1, min_ks_v1, ks_distances_v1, S_values_v1 = find_optimal_S_grid(treated_df=temp_treated_v1, untreated_df=temp_untreated_v1,
        #                                           cutoff_tr=self.cell_number_cutoff)
        # best_S_v3, min_ks_v3, ks_distances_v3, S_values_v3 = find_optimal_S_grid(treated_df=temp_treated_v3, untreated_df=temp_untreated_v3,
        #                                           cutoff_tr=self.cell_number_cutoff)
        # best_S_v1, max_overlap_v1, overlaps_v1 = find_optimal_S_max_auc(treated_df=temp_treated_v1, vehicle_df=temp_untreated_v1,
        #                                                         cutoff_tr=self.cell_number_cutoff)
        # best_S_v3, max_overlap_v3, overlaps_v3 = find_optimal_S_max_auc(treated_df=temp_treated_v3, vehicle_df=temp_untreated_v3,
        #                                                         cutoff_tr=self.cell_number_cutoff)
        # return {
        #     'Bootstrap_id': f'B{bs_cycle}',
        #     'cutoff_treated': self.cell_number_cutoff,
        #     'Shrinkage_v1': best_S_v1,
        #     'distance_v1': min_ks_v1,
        #     'Shrinkage_v3': best_S_v3,
        #     'distance_v3': min_ks_v3
        # }
        # return {
        #     'Bootstrap_id': f'B{bs_cycle}',
        #     'cutoff_treated': self.cell_number_cutoff,
        #     'Shrinkage_v1': best_S_v1,
        #     'overlap_v1': max_overlap_v1,
        #     'Shrinkage_v3': best_S_v3,
        #     'overlap_v3': max_overlap_v3
        # }

    # def _calculate_shrinkage(self) -> List[Dict]:   
    #     "max auc"
    #     best_S, max_overlap, overlaps = find_optimal_S_max_auc(treated_df=self.treated_df, vehicle_df=self.untreated_df,
    #                                                cutoff_tr=self.cell_number_cutoff)
    #     return [{
    #         'Shrinkage': best_S,
    #         'Bootstrap_id': 'Real',
    #         'cutoff_treated': self.cell_number_cutoff,
    #         'max_overlap': max_overlap
    #     }]
        
    # def _perform_single_bootstrap(self, treated_indices: Dict, untreated_indices: Dict, 
    #                               total_gRNA_number: int, bs_cycle='') -> Dict:
    #     # max auc
    #     np.randome.seed(2025)
    #     treated_sample = Nested_Bootstrap_Index_single(treated_indices)
    #     untreated_sample = Nested_Bootstrap_Index_Special_single(untreated_indices, self.untreated_df, total_gRNA_number)
        
    #     temp_treated = self.treated_df.loc[treated_sample]
    #     temp_untreated = self.untreated_df.loc[untreated_sample]
        
    #     best_S, max_overlap, overlaps = find_optimal_S_max_auc(treated_df=temp_treated, vehicle_df=temp_untreated,
    #                                                cutoff_tr=self.cell_number_cutoff)
    #     return {
    #         'Shrinkage': best_S,
    #         'Bootstrap_id': f'B{bs_cycle}',
    #         'cutoff_treated': self.cell_number_cutoff,
    #         'max_overlap': max_overlap
    #     }