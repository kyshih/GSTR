'''Combine all mice to make one "super" mouse and calculate its percentile tumor sizes.
''' 
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import ks_2samp

def calculate_percentiles(df, cutoff, gRNAs, percentiles, group_col):
    """
    Calculate specified percentile tumor size for a given cutoff.
    """
    filtered_df = df[(df['gRNA'].isin(gRNAs))&(df['Cell_number']>cutoff)]
    
    if filtered_df.empty:
        print("cutoff is too high. df is empty.")
        return None
    if group_col:
        # Calculate percentiles within each group and take the median across groups
        return filtered_df.groupby(group_col)['Cell_number'].apply(lambda x: np.percentile(x, percentiles))
    else:
        return np.percentile(filtered_df['Cell_number'], percentiles)

def ks_loss(treated_df, untreated_df, cutoff_tr, base_cutoff, input_control_gRNA_list):
    treated_filtered_df = treated_df[(treated_df['Cell_number']>cutoff_tr)&(treated_df['gRNA'].isin(input_control_gRNA_list))]['Cell_number']
    untreated_filtered_df = untreated_df[(untreated_df['Cell_number']>base_cutoff)&(untreated_df['gRNA'].isin(input_control_gRNA_list))]['Cell_number']
    ks_stat, _ = ks_2samp(treated_filtered_df, untreated_filtered_df)
    print(f'ks stats is {ks_stat} and {_}')
    return ks_stat

def objective(S, treated_df, untreated_df, input_control_gRNA_list, base_cutoff, percentiles, group_col, loss):
    """
    Objective function to minimize to match the percentile distribution of treated and vehicle groups.
    """
    # Adjust cutoff for treated group
    cutoff_tr = base_cutoff * S
    #cutoff_tr = base_cutoff / S
    
    if loss == 'ks':
        score = ks_loss(treated_df, untreated_df, cutoff_tr, base_cutoff, input_control_gRNA_list)
        return score
    # Calculate percentiles for both treated and vehicle groups
    treated_percentiles = calculate_percentiles(treated_df, cutoff_tr, input_control_gRNA_list, percentiles, group_col)
    vehicle_percentiles = calculate_percentiles(untreated_df, base_cutoff, input_control_gRNA_list, percentiles, group_col)
    
    if treated_percentiles is None or vehicle_percentiles is None:
        print("Skipping this S value due to empty untreated group after filtering.")
        return np.inf
    
    # Calculate loss
    if loss == "L2":
        score = np.sum((treated_percentiles - vehicle_percentiles) ** 2)
    elif loss == "L1":
        score = np.sum(abs(treated_percentiles - vehicle_percentiles))    
    return score

def minimize_scalar_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, percentiles, group_col, loss="ks"):
    # Use minimize_scalar to find the optimal S
    result = minimize_scalar(
        objective,
        bounds=(0.01, 200),  # Adjust bounds as necessary for your data
        args=(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, percentiles, group_col, loss),
        method='bounded',
        options={'xatol': 1e-4}
    )
    if result.success:
        optimal_S = result.x
        optimal_cutoff_tr = base_cutoff * optimal_S
        error = result.fun
        return optimal_S, optimal_cutoff_tr, error
    else:
        raise ValueError("Optimization failed. Check the input data and parameters.")    

def GD_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, percentiles, group_col,
                       learning_rate=0.001, max_iter=10000, tolerance=300, epsilon=0.5, loss="L1"):
    S = 1
    for i in range(max_iter):
        error = objective(S, treated_df, untreated_df, input_control_gRNA_list, base_cutoff, percentiles, group_col, loss)
        #print(f'error is {error}')
        cutoff_tr = base_cutoff / S
        if error <= tolerance:
            return S, cutoff_tr, error
    
        # Approx gradient
        error_with_epsilon = objective(S + epsilon, treated_df, untreated_df, input_control_gRNA_list, base_cutoff, percentiles, group_col, loss)
        gradient = (error_with_epsilon - error) / epsilon
        print(error)
        print(error_with_epsilon)
        print(S)
        print(f'Gradient is {gradient}')
        S -= learning_rate * gradient
    return S, cutoff_tr, error
    
def find_S(treated_df, untreated_df, gRNAs, base_cutoff, percentiles=[80, 90, 95], group_col=None,
           method='GD', loss="L1"):
    """
    Find the optimal S value to match the percentile distribution of treated and vehicle groups.
    """
    if method == "minimize_scalar":
        return minimize_scalar_S(treated_df, untreated_df, gRNAs, base_cutoff, percentiles, group_col, loss)
    elif method == "GD":
        return GD_S(treated_df, untreated_df, gRNAs, base_cutoff, percentiles, group_col, loss)
    else:
        raise ValueError("Invalid method. Choose ''gradient_descent', or 'minimize_scalar'.")
