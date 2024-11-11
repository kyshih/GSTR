'''Combine all mice to make one "super" mouse and calculate its percentile tumor sizes.
''' 
import numpy as np
from scipy.optimize import minimize_scalar

def calculate_percentiles(df, cutoff, gRNAs, percentiles, group_col):
    """
    Calculate specified percentile tumor size for a given cutoff.
    """
    filtered_df = df[(df['gRNA'].isin(gRNAs))&(df['Cell_number']>cutoff)]
    #return np.percentile(filtered_df.groupby(group_col)['Cell_number'], percentiles)
    return np.percentile(filtered_df['Cell_number'], percentiles)

def find_S(treated_df, untreated_df, gRNAs, base_cutoff, percentiles=[50, 70, 80, 90, 95], group_col='Sample_ID'):
    """
    Find the optimal S value to match the percentile distribution of treated and vehicle groups.
    """
    
    def objective(S):
        # Adjust cutoff for treated group
        cutoff_tr = base_cutoff * S
        
        # Calculate percentiles for both treated and vehicle groups
        treated_percentiles = calculate_percentiles(treated_df, cutoff_tr, gRNAs, percentiles, group_col)
        vehicle_percentiles = calculate_percentiles(untreated_df, base_cutoff, gRNAs, percentiles, group_col)
        
        # Calculate the sum of squared differences between treated and vehicle percentiles
        score = np.sum((treated_percentiles - vehicle_percentiles) ** 2)
        return score
    
    # Use minimize_scalar to find the optimal S
    result = minimize_scalar(
        objective,
        bounds=(0.01, 200),  # Adjust bounds as necessary for your data
        method='bounded',
        options={'xatol': 1e-4}
    )
    
    if result.success:
        optimal_S = result.x
        optimal_cutoff_tr = base_cutoff * optimal_S
        SE = result.fun
        #print(f'Optimal S: {optimal_S}')
        #print(f'Optimal Treated Cutoff: {optimal_cutoff_tr}')
        #print(f'Minimized SE Error: {SE}')
        return optimal_S, optimal_cutoff_tr, SE
    else:
        raise ValueError("Optimization failed. Check the input data and parameters.")    