import numpy as np
from scipy.optimize import minimize_scalar
from scipy import stats

def calculate_density_overlap(treated_df, vehicle_df, cutoff_tr, S, bins=100):
    """
    Calculate overlap between density distributions.
    """
    treated_data = treated_df[
        (treated_df['Cell_number'] > cutoff_tr)
    ]['Cell_number']

    # Shrink vehicle data
    vehicle_df_copy = vehicle_df.copy()
    vehicle_df_copy['Cell_number'] = vehicle_df_copy['Cell_number'] * S
    # apply cutoff
    vehicle_data = vehicle_df_copy[
        (vehicle_df_copy['Cell_number'] > cutoff_tr)
    ]['Cell_number']
    
    if len(treated_data) == 0 or len(vehicle_data) == 0:
        print(f"Empty dataset found for S={S}")
        return -np.inf
    
    # Print some debug info
    # print(f"\nCalculating overlap for S={S}")
    # print(f"Treated data range: {treated_data.min():.2f} to {treated_data.max():.2f}")
    # print(f"Vehicle data range (after scaling): {vehicle_data.min():.2f} to {vehicle_data.max():.2f}")
    
    # Get the range for density estimation
    min_val = min(treated_data.min(), vehicle_data.min())
    max_val = max(treated_data.max(), vehicle_data.max())
    x_range = np.logspace(np.log10(min_val), np.log10(max_val), bins)
    dx = np.diff(np.log10(x_range))[0]
    
    # Calculate kernel density estimates with adaptive bandwidth
    treated_kde = stats.gaussian_kde(np.log10(treated_data), bw_method='silverman')
    vehicle_kde = stats.gaussian_kde(np.log10(vehicle_data), bw_method='silverman')
    
    # Evaluate densities
    treated_density = treated_kde(np.log10(x_range))
    vehicle_density = vehicle_kde(np.log10(x_range))
    
    # Calculate overlap
    overlap = np.minimum(treated_density, vehicle_density).sum() * dx
    # print(f"Calculated overlap: {overlap}")
    
    return -overlap

def find_optimal_S_max_auc(treated_df, vehicle_df, cutoff_tr, bins=100, post_clip=False):
    """
    Find optimal shrinkage factor S with more detailed debugging.
    """
    def objective(S):
        if not post_clip:
            return calculate_density_overlap(treated_df, vehicle_df, cutoff_tr, S, bins)
        elif post_clip:
            return calculate_density_overlap_post_clip(treated_df, vehicle_df, cutoff_tr, S, bins)
    
    # Try multiple starting points with different bounds
    S_candidates = []
    overlaps = []
    
    initial_points = [0.2, 0.5, 1.0, 1.5]
    for init_S in initial_points:
        # print(f"\nTrying optimization with initial S={init_S}")
        
        # Use tighter bounds around the initial point
        lower_bound = max(0.01, init_S * 0.5)
        upper_bound = min(5.0, init_S * 2.0)
        
        result = minimize_scalar(
            objective,
            bounds=(lower_bound, upper_bound),
            method='bounded',
            options={'xatol': 1e-4, 'maxiter': 100}
        )
        
        # print(f"Optimization result for initial S={init_S}:")
        # print(f"  Final S: {result.x}")
        # print(f"  Final overlap: {-result.fun}")
        # print(f"  Success: {result.success}")
        # print(f"  Number of iterations: {result.nfev}")
        
        if result.success:
            S_candidates.append(result.x)
            overlaps.append(result.fun)
    
    if not S_candidates:
        raise ValueError("Optimization failed for all starting points")
    
    # Print all results
    # print("\nAll optimization results:")
    # for i, (S, overlap) in enumerate(zip(S_candidates, overlaps)):
    #     print(f"Candidate {i+1}: S={S:.4f}, overlap={-overlap:.4f}")
    
    # Select the best result
    best_idx = np.argmin(overlaps)
    return S_candidates[best_idx], -overlaps[best_idx], overlaps

