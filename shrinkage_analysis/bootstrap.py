# shrinkage_analysis/bootstrap.py

import pandas as pd
import numpy as np
import importlib
from typing import List, Dict, Tuple
from bootstrapping_helpers import (
    Generate_ref_input_df,
    Generate_Index_Dictionary,
    Nested_Bootstrap_Index_single,
    Nested_Bootstrap_Index_Special_single
)
#from shrinkage_analysis.find_S_matching_distribution import find_S

class BootstrapAnalyzer:
    def __init__(self, raw_treated_df: pd.DataFrame, raw_untreated_df: pd.DataFrame,
                 control_gRNA_list: List[str], cell_number_cutoff: int, find_S_version):
        self.treated_df = raw_treated_df
        self.untreated_df = raw_untreated_df
        self.control_gRNA_list = control_gRNA_list
        self.cell_number_cutoff = cell_number_cutoff
        self.find_S = find_S_version

    def generate_bootstrap_samples(self, n_replicates: int, total_gRNA_number: int, obj: str, method: str) -> List[Dict]:
        results = self._calculate_shrinkage(obj=obj, method=method)
        
        if n_replicates == 0:
            return results

        treated_indices = Generate_Index_Dictionary(self.treated_df)
        untreated_indices = Generate_Index_Dictionary(self.untreated_df)

        for i in range(n_replicates):
            bootstrap_results = self._perform_single_bootstrap(
                treated_indices, untreated_indices, total_gRNA_number, obj, method
            )
            results.append(bootstrap_results)

        return results

    def _calculate_shrinkage(self, obj, method) -> List[Dict]:
        S, adjusted_cutoff, error = self.find_S(treated_df=self.treated_df, untreated_df=self.untreated_df,
                                                gRNAs=self.control_gRNA_list, base_cutoff=self.cell_number_cutoff, 
                                                method=method, loss=obj)
        return [{
            'Shrinkage': S,
            'Bootstrap_id': 'Real',
            'basal_cutoff': self.cell_number_cutoff,
            'adjusted_cutoff': adjusted_cutoff,
            'error': error
        }]

    def _perform_single_bootstrap(self, treated_indices: Dict, untreated_indices: Dict, 
                                  total_gRNA_number: int, obj: str, method: str) -> Dict:
        treated_sample = Nested_Bootstrap_Index_single(treated_indices)
        untreated_sample = Nested_Bootstrap_Index_Special_single(untreated_indices, self.untreated_df, total_gRNA_number)
        
        temp_treated = self.treated_df.loc[treated_sample]
        temp_untreated = self.untreated_df.loc[untreated_sample]
        
        #S, adjusted_cutoff, SE = find_optimal_S(temp_treated, temp_untreated, self.control_gRNA_list, self.cell_number_cutoff)
        S, adjusted_cutoff, error = self.find_S(temp_treated, temp_untreated, self.control_gRNA_list, self.cell_number_cutoff, method, obj)
        return {
            'Shrinkage': S,
            'Bootstrap_id': f'B{len(treated_sample)}',  # Using len(treated_sample) as a unique ID
            'basal_cutoff': self.cell_number_cutoff,
            'adjusted_cutoff': adjusted_cutoff,
            'error': error
        }
