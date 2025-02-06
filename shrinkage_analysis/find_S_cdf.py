'''Combine all mice to make one "super" mouse for treated and untreated groups. Focus on inert tumors to establish E
    1. I apply some cutoff in the treatment arm
    2. shrinkage untreated group by S
    3. apply cutoff in the shrunk untreated group
    4. evalutae treated and shrunk untreated similarity
''' 
import numpy as np
from scipy.stats import ks_2samp

def find_optimal_S_grid(treated_df, untreated_df, cutoff_tr, n_points=100, method='ks'):
    """
    Find shrinkage S using grid search.
    """
    S_values = np.logspace(-1, 0.5, n_points)  # try values from 0.01 to 3.16
    if method == 'ks':
        ks_distances = []
        for S in S_values:
            #overlap = calculate_density_overlap(treated_df, untreated_df, cutoff_tr, S)
            ks_distance = calculate_ks_stats(treated_df, untreated_df, cutoff_tr, S)
            ks_distances.append(ks_distance)
        
        min_ks = min(ks_distances)
        
        # Find all S values with the minimum KS statistic
        tied_S = [S_values[i] for i, res in enumerate(ks_distances) if res == min_ks]
        # Select the highest S among the ties
        #best_S = max(tied_S)
        best_S = np.mean(tied_S)
        return best_S, min_ks, ks_distances, S_values
    
    elif method == 'cdf_auc':
        overlaps = []
        for S in S_values:
            overlap = compute_cdf_overlap(treated_df, untreated_df, cutoff_tr, S)
            overlaps.append(overlap)
        max_overlap = np.max(overlaps)
        tied_S = [S_values[i] for i in range(len(overlaps)) if overlaps[i] == max_overlap]
        best_S = np.mean(tied_S)
    #return S_values[best_idx], -overlaps[best_idx]
        return best_S, max_overlap, overlaps, S_values

def calculate_ks_stats(treated_df, untreated_df, cutoff_tr, S):
    """
    Calculate KS statistic between treated and scaled untreated group.
    """
    treated_data = treated_df[(treated_df['Cell_number'] > cutoff_tr)]['Cell_number']
    untreated_df_copy = untreated_df.copy(deep=True)
    untreated_df_copy['Cell_number'] *= S
    untreated_data = untreated_df_copy[(untreated_df_copy['Cell_number'] > cutoff_tr)]['Cell_number']

    if len(treated_data) == 0 or len(untreated_data) == 0:
        print(untreated_data)
        print(f"Empty dataset found for S={S}")
        return np.inf
    
    # Calculate KS statistic
    ks_stat, _ = ks_2samp(treated_data, untreated_data)
    
    return ks_stat 
    
def compute_cdf_overlap(treated_df, untreated_df, cutoff_tr, S):
    treated_data = treated_df[treated_df['Cell_number']>cutoff_tr]['Cell_number']
    untreated_df_copy = untreated_df.copy()
    untreated_df_copy = untreated_df_copy['Cell_number'] * S
    untreated_data = untreated_df_copy[untreated_df_copy['Cell_number']>cutoff_tr]['Cell_number']
    
    # Combine all points for evaluation
    all_points = np.sort(np.concatenate([treated_data, untreated_data]))
    
    # Compute empirical CDFs at all points
    cdf_treated = stats.ecdf(treated_data)
    cdf_untreated = stats.ecdf(untreated_data)
    
    f1 = cdf_treated(all_points)
    f2 = cdf_untreated(all_points)
    
    # Compute overlap using trapezoidal integration
    overlap = np.trapz(np.minimum(f1, f2), all_points)
    
    return overlap