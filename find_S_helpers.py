""" This module includes different versions of find_S functions to estimate GxE
    test out these functions for your experiment to decide which one to use. You might need to modify learning rate
    to avoid oscillations of errors esp in SE loss and in L1 loss to improve rate of convergence
"""
from scipy.optimize import minimize_scalar

# gRNA level median inert TTN
def find_S_given_base_cutoff(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, 
           upper=20, lower=0.01, precision=0.001, max_iter=100, tolerance=5, group_col=['Numbered_gene_name']):
    """iterative function to find S based on given base cutoff (L) in untreated group
    move L' in the treated group to match median TTN of inert tumors in untreated group
    adjusted cutoff L' = base cutoff L * S
    Args:
        treated_df: raw treated df pre cutoff
        untreated_df: raw untreated df pre cutoff
        base_cutoff: base cutoff L in untreated group
        upper: upper limit for Shrinkage S. Defaults to 20.
        lower: lower limit for Shrinkage S. Defaults to 0.01.
        precision (float, optional): upper - lower. Defaults to 0.001.
        tolerance: stopping condition for median(inert TTN|untreated) - median(inert TTN|treated). Defaults to 2.

    Returns:
        S: Shrinkage S
        cutoff_tr: cutoff for treated group
    """
    print(f'in find_S_given_base_cutoff, inert guides are {input_control_gRNA_list}')
    untreated_inert_df = untreated_df[(untreated_df['gRNA'].isin(input_control_gRNA_list))&(untreated_df['Cell_number']>base_cutoff)]
    treated_inert_df = treated_df[(treated_df['gRNA'].isin(input_control_gRNA_list))]
    N_control_untreated = untreated_inert_df.groupby(group_col).Clonal_barcode.count().median()
    print(f'N inert untreated is {N_control_untreated}')
    i = 0
    if upper < lower:
        print(f"upper bound must be higher than the lower bound")
        return -1
    while upper-lower >= precision and i < max_iter:
        S = (upper + lower) / 2
        cutoff_tr = base_cutoff * S
        N_control_treated = treated_inert_df[treated_inert_df['Cell_number']>cutoff_tr].groupby(group_col).Clonal_barcode.count().median()
        print(f'N inert treated is {N_control_treated}')
        error = N_control_treated - N_control_untreated
        
        if abs(error) <= tolerance:
            return S, cutoff_tr
        if N_control_treated > N_control_untreated:
            # cutoff for treated is too low
            # increase S to increase cutoff_tr
            lower = S
        else:
            # cutoff for treated is too high
            # decrease S to decrease cutoff_tr
            upper = S
        i += 1
    return S, cutoff_tr, error # cutoff in untreated will be basal cutoff

def find_S_gradient_descent_L1_loss_given_base_cutoff(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, 
                            learning_rate=0.05, max_iter=10000, tolerance=5, group_col=['Numbered_gene_name']):
    """GD to find the optimal shrinkage factor S using L1 loss.

    Args:
        treated_df: raw treated df pre cutoff
        untreated_df: raw untreated df pre cutoff
        base_cutoff: Cell number cutoff for the untreated group given by me
        learning_rate: Step size for gradient descent
        max_iter: Maximum number of iterations
        tolerance: Minimum error difference to stop iteration

    Returns:
        S, cutoff_unt: Optimal shrinkage factor, Cutoff for the untreated group
    """
    untreated_inert_df = untreated_df[(untreated_df['gRNA'].isin(input_control_gRNA_list)) & (untreated_df['Cell_number'] > base_cutoff)]
    treated_inert_df = treated_df[treated_df['gRNA'].isin(input_control_gRNA_list)]
    N_control_untreated = untreated_inert_df.groupby(group_col).Clonal_barcode.count().median()

    S = 1  # Initial guess
    for i in range(max_iter):
        cutoff_tr = base_cutoff * S
        N_control_treated = treated_inert_df[treated_inert_df['Cell_number'] > cutoff_tr].groupby(group_col).Clonal_barcode.count().median()
        
        error = N_control_untreated - N_control_treated
        print(f'Iteration {i+1}, S = {S}, L1 error = {error}')
        if abs(error) <= tolerance:
            return S, cutoff_tr
        
        gradient = error / N_control_untreated  # Adjust step size based on relative error
        S -= learning_rate * gradient  # Adjust S based on gradient (step size)
    return S, cutoff_tr

def find_S_gradient_descent_se_given_base_cutoff(treated_df, untreated_df, input_control_gRNA_list, base_cutoff,
                                learning_rate=0.005, max_iter=5000, tolerance=2500, group_col=['Numbered_gene_name']):
    """
    Gradient descent to find optimal shrinkage factor S using squared error (SE).
    
    Args:
        treated_df: raw treated df pre cutoff
        untreated_df: raw untreated df pre cutoff
        adjusted_cutoff: Cell number cutoff for the treated group
        learning_rate: Step size for gradient descent
        max_iter: Maximum number of iterations
        tolerance: Minimum error difference to stop iteration
    
    Returns:
        S: Optimal shrinkage factor
        cutoff_unt: Cutoff for untreated group (L = L' / S)
    """
    untreated_inert_df = untreated_df[(untreated_df['gRNA'].isin(input_control_gRNA_list)) & (untreated_df['Cell_number'] > base_cutoff)]
    treated_inert_df = treated_df[treated_df['gRNA'].isin(input_control_gRNA_list)]
    
    # Calculate the treated group's inert tumors median
    N_control_untreated = untreated_inert_df.groupby(group_col).Clonal_barcode.count().median()
    print(f'median inert TTN untreated: {N_control_untreated}')
    
    S = 1  # Initial guess for shrinkage factor
    for i in range(max_iter):
        cutoff_tr = base_cutoff * S
        N_control_treated = treated_inert_df[treated_inert_df['Cell_number'] > cutoff_tr].groupby(group_col).Clonal_barcode.count().median()
        print(f'Median inert TTN treated: {N_control_treated}')
        
        # Compute Squared Error (SE)
        error = (N_control_untreated - N_control_treated)
        se_error = error ** 2  # Squared error
        
        print(f'Iteration {i}: S = {S}, SE Error = {se_error}')
        # If the error is within the tolerance, we can stop
        if se_error <= tolerance:
            return S, cutoff_tr, se_error
        
        # Gradient is proportional to the error
        gradient = -2 * error  # Derivative of (N_control_treated - N_control_untreated) ** 2
        print(f'Grandient: {gradient}')
        # Update S based on the gradient
        S += learning_rate * gradient
    
    return S, cutoff_tr, se_error  # Return the final value of S and cutoff for untreated group

def objective(S, treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col):
    """
    Objective function to compute squared error for a given S.
    """
    cutoff_tr = base_cutoff * S
    # Filter treated samples above cutoff_tr
    treated_filtered = treated_df[
        (treated_df['gRNA'].isin(input_control_gRNA_list)) & 
        (treated_df['Cell_number'] > cutoff_tr)
    ]
    if treated_filtered.empty:
        # Assign a high error if no samples exceed the cutoff
        return float('inf')
    # Compute median Clonal_barcode count for treated
    N_control_treated = treated_filtered.groupby(group_col)['Clonal_barcode'].count().median()
    
    # Filter untreated samples above base_cutoff
    untreated_filtered = untreated_df[
        (untreated_df['gRNA'].isin(input_control_gRNA_list)) & 
        (untreated_df['Cell_number'] > base_cutoff)
    ]
    N_control_untreated = untreated_filtered.groupby(group_col)['Clonal_barcode'].count().median()
    
    # Compute squared error
    error = N_control_untreated - N_control_treated
    se_error = error ** 2
    return se_error

def find_S_min_scalar(treated_df, untreated_df, input_control_gRNA_list, base_cutoff,group_col=['Sample_ID']):
    """
    Finds the optimal S using scipy's minimize_scalar. Default maxiter is 500
    """
    #print(f'group col is {group_col}')
    result = minimize_scalar(
        objective,
        bounds=(0.01, 100),  # Adjust bounds based on domain knowledge
        args=(treated_df, untreated_df, input_control_gRNA_list, base_cutoff, group_col),
        method='bounded',
        options={'xatol': 1e-4}
    )
    
    if result.success:
        optimal_S = result.x
        optimal_cutoff_tr = base_cutoff * optimal_S
        SE = result.fun
        print(f'Optimal S: {optimal_S}')
        print(f'Optimal Treated Cutoff (cutoff_tr): {optimal_cutoff_tr}')
        print(f'Minimized SE Error: {SE}')
        return optimal_S, optimal_cutoff_tr, SE
    else:
        raise ValueError("Optimization failed. Check the input data and parameters.")