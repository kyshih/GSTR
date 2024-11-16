from scipy.optimize import minimize_scalar
import numpy as np

''' This module calculates the treatment response S by matching median inert TTN per mouse in treated and untreated. 
    - L (base cutoff) is given, iteratively find L' (adjusted cutoff) in treated.
    - Note: GD params might need to change to get values within desired tolerance and precision but balanced with computation time.
        - For GD, learning rate=0.05 and tolerance = 200 (for SE) or 20 (or 50 for L1) seem to work well
        - if wanna match median inert tumors across mice, the assumption is that titers (hence TTN/mouse) and mouse number
        per treatment arm is the same
'''

def calculate_error(treated_df, untreated_df, input_control_gRNA_list, cutoff_tr, base_cutoff, group_col, loss):
    """
    Calculate the median TTN difference (error) and the corresponding loss error (SE or L1) between treated and untreated groups.
    """
    untreated_filtered = untreated_df[
        (untreated_df['gRNA'].isin(input_control_gRNA_list)) & 
        (untreated_df['Cell_number'] > base_cutoff)
    ]
    treated_inert_df= treated_df[
        (treated_df['gRNA'].isin(input_control_gRNA_list))]
    
    # if treated_filtered.empty or untreated_filtered.empty:
    #     return None, float('inf')  # Return high error if no samples meet the criteria

    N_control_untreated = untreated_filtered.groupby(group_col)['Clonal_barcode'].count().median()
    N_control_treated = treated_inert_df[treated_inert_df['Cell_number']>cutoff_tr].groupby(group_col)['Clonal_barcode'].count().median()
    
    error = N_control_treated - N_control_untreated
    if loss == "SE":
        loss_error = error ** 2
    elif loss == "L1":
        loss_error = abs(error)
    else:
        raise ValueError("Invalid loss function. Choose 'SE' or 'L1'.")
    
    return error, loss_error, N_control_untreated

def objective(S, treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss):
    """
    Objective function to compute SE or L1 error for a given S using specified loss function.
    """
    cutoff_tr = base_cutoff * S
    _, loss_error, _ = calculate_error(treated_df, untreated_df, input_control_gRNA_list, cutoff_tr, base_cutoff, group_col, loss)
    return loss_error

def binary_search_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss="L1",
                    precision=0.001, tolerance=5, max_iter=100):
    lower, upper = 0.01, 20  # Define bounds
    i = 0
    while upper-lower >= precision and i < max_iter:
        S = (upper + lower) / 2
        cutoff_tr = base_cutoff * S
        error, loss_error, _ = calculate_error(treated_df, untreated_df, input_control_gRNA_list, cutoff_tr, base_cutoff, group_col, loss)
        if abs(error) <= tolerance: # error = N_control_untreated - N_control_treated
            return S, cutoff_tr, loss_error
        elif error > 0: # N_control_treated is too low. cutoff_tr is too high. decrease S to decrease cutoff_tr
            lower = S
        else:
            upper = S
        i += 1
    return S, cutoff_tr, loss_error

def gradient_descent_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss,
                       learning_rate=0.05, max_iter=10000, tolerance=5):
    S = 1  # Initial guess
    for i in range(max_iter):
        cutoff_tr = base_cutoff * S
        error, loss_error, N_control_untreated = calculate_error(treated_df, untreated_df, input_control_gRNA_list, cutoff_tr, base_cutoff, group_col, loss)
        
        if loss_error <= tolerance:
            return S, cutoff_tr, loss_error
        
        if loss == "SE":
            gradient = -2 * error
        elif loss == "L1":
            gradient = -(error / N_control_untreated)
            #gradient = -1
        else:
            raise ValueError("Invalid loss function. Choose 'SE' or 'L1'.")
        
        S -= learning_rate * gradient
    return S, cutoff_tr, loss_error

def minimize_scalar_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss):
    result = minimize_scalar(
        objective,
        bounds=(0.01, 100),
        args=(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss),
        method='bounded',
        options={'xatol': 1e-4}
    )
    
    if result.success:
        S = result.x
        cutoff_tr = base_cutoff * S
        loss_error = result.fun
        return S, cutoff_tr, loss_error
    else:
        raise ValueError("Optimization failed. Check the input data and parameters.")

def find_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col=['Sample_ID'],
           method="gradient_descent", loss="L1", tolerance=50):
           #learning_rate=0.003, max_iter=10000):
    """
    Wrapper to find optimal S using specified method and loss function.
    """
    if method == "binary_search":
        return binary_search_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col)
    elif method == "gradient_descent":
        return gradient_descent_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss)
    elif method == "minimize_scalar":
        return minimize_scalar_S(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col, loss)
    else:
        raise ValueError("Invalid method. Choose 'binary_search', 'gradient_descent', or 'minimize_scalar'.")