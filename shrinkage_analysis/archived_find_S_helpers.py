'''Archived find_S functions
    These functions set adjusted cutoff in treated group, and find the matching base cutoff in the untreated group.
    find L (base cutoff) given L' (adjusted cutoff in treated)
'''

def find_S(treated_df, untreated_df, input_control_gRNA_list, adjusted_cutoff, 
           upper=20, lower=0.01, precision=0.001, max_iter=100, tolerance=5, group_col=['Numbered_gene_name']):
    """
        iterative
        move L in the untreated group to match median TTN of inert tumors in treated group
        binary search
        
        Args:
        treated_df: raw treated df pre cutoff
        untreated_df: raw untreated df pre cutoff
        adjusted_cutoff: cell number cutoff for treated group = L'
        tolerance: tolerance for matching median TTN of inert tumors in treated group wrt untreated
        
        Return:
        S: shrinkage, interaction term
        cutoff_unt: cutoff for untreated group = L (L = L' / S)
    """
    print(f'in find_S, inert guides are {input_control_gRNA_list}')
    treated_inert_df = treated_df[(treated_df['gRNA'].isin(input_control_gRNA_list))&(treated_df['Cell_number']>adjusted_cutoff)]
    untreated_inert_df = untreated_df[untreated_df['gRNA'].isin(input_control_gRNA_list)]
    N_control_treated = treated_inert_df.groupby(group_col).Clonal_barcode.count().median()
    print(f'N inert treated is {N_control_treated}')
    i = 0
    if upper < lower:
        print(f"upper bound must be higher than the lower bound")
        return -1
    while upper-lower >= precision and i < max_iter:
        S = (upper + lower) / 2
        cutoff_unt = adjusted_cutoff / S
        N_control_untreated = untreated_inert_df[untreated_inert_df['Cell_number']>cutoff_unt].groupby(group_col).Clonal_barcode.count().median()
        print(f'N inert untreated is {N_control_untreated}')
	
        if abs(N_control_treated - N_control_untreated) <= tolerance:
            return S, cutoff_unt
        if N_control_treated > N_control_untreated:
            # cutoff for untreated is too high
            # increase S to decrease cutoff_unt
            lower = S
        else:
            # cutoff for untreated is too low
            # decrease S to increase cutoff_unt
            upper = S
        i += 1
    return S, cutoff_unt # cutoff in untreated will be basal cutoff

def find_S_gradient_descent_L1_loss(treated_df, untreated_df, input_control_gRNA_list, adjusted_cutoff, 
                            learning_rate=0.01, max_iter=5000, tolerance=5, group_col=['Numbered_gene_name']):
    """GD to find the optimal shrinkage factor S using L1 loss.

    Args:
        treated_df: raw treated df pre cutoff
        untreated_df: raw untreated df pre cutoff
        input_control_gRNA_list: List of inert gRNAs
        adjusted_cutoff: Cell number cutoff for the treated group given by me
        learning_rate: Step size for gradient descent
        max_iter: Maximum number of iterations
        tolerance: Minimum error difference to stop iteration

    Returns:
        S, cutoff_unt: Optimal shrinkage factor, Cutoff for the untreated group
    """
    treated_inert_df = treated_df[(treated_df['gRNA'].isin(input_control_gRNA_list)) & (treated_df['Cell_number'] > adjusted_cutoff)]
    untreated_inert_df = untreated_df[untreated_df['gRNA'].isin(input_control_gRNA_list)]
    N_control_treated = treated_inert_df.groupby(group_col).Clonal_barcode.count().median()

    S = 1  # Initial guess
    for i in range(max_iter):
        cutoff_unt = adjusted_cutoff / S
        N_control_untreated = untreated_inert_df[untreated_inert_df['Cell_number'] > cutoff_unt].groupby(group_col).Clonal_barcode.count().median()
        
        error = N_control_treated - N_control_untreated
        print(f'Iteration {i+1}, S = {S}, L1 error = {error}')
        if abs(error) <= tolerance:
            return S, cutoff_unt
        
        gradient = error / N_control_treated  # Adjust step size based on relative error
        S += learning_rate * gradient  # Adjust S based on gradient (step size)
    return S, cutoff_unt

def find_S_gradient_descent_se(treated_df, untreated_df, input_control_gRNA_list, adjusted_cutoff, 
                                learning_rate=0.001, max_iter=5000, tolerance=25, group_col=['Numbered_gene_name']):
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
    treated_inert_df = treated_df[(treated_df['gRNA'].isin(input_control_gRNA_list)) & (treated_df['Cell_number'] > adjusted_cutoff)]
    untreated_inert_df = untreated_df[untreated_df['gRNA'].isin(input_control_gRNA_list)]
    
    # Calculate the treated group's inert tumors median
    N_control_treated = treated_inert_df.groupby(group_col).Clonal_barcode.count().median()
    
    S = 1  # Initial guess for shrinkage factor
    for i in range(max_iter):
        cutoff_unt = adjusted_cutoff / S
        N_control_untreated = untreated_inert_df[untreated_inert_df['Cell_number'] > cutoff_unt].groupby(group_col).Clonal_barcode.count().median()

        # Compute Squared Error (SE)
        error = (N_control_treated - N_control_untreated)
        se_error = error ** 2  # Squared error
        
        print(f'Iteration {i}: S = {S}, SE Error = {se_error}')

        # If the error is within the tolerance, we can stop
        if se_error <= tolerance:
            return S, cutoff_unt
        
        # Gradient is proportional to the error
        gradient = -2 * error  # Derivative of (N_control_treated - N_control_untreated) ** 2
        S -= learning_rate * gradient
    
    return S, cutoff_unt  # Return the final value of S and cutoff for untreated group

