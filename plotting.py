"""This module plot trait comparions
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec

def generate_symmetric_list(input_n):
    # generate a symmetric list of numbers centered around 0, used to distribute points symmetrically around an axis when plotting
    if input_n % 2 == 1:
        return(np.array(range(-int((input_n-1)/2),int((input_n+1)/2)))*0.1)
    else:
        return(np.array(range(-int(input_n/2),int(input_n/2)))*0.1+0.05)

def plot_trait_comparison(df1, df2, df3, df4, trait, title1, title2, title3, title4, y_label, fig_size=(22, 20),
                          logy=True, ymin=None, ymax=None, fig_name=None, output_path=None):
    """
    Plot trait comparison across three datasets with error bars and scatter plots.

    Parameters:
    - df1, df2, df3: DataFrames containing the data to be plotted.
    - trait: The trait/column of interest for plotting on the y-axis.
    - title1, title2, title3: Titles for the three subplots.
    - y_label: Label for the y-axis.
    - fig_name: Optional. The name to save the figure as (include file extension, e.g., 'plot.png').
    - output_path: Optional. The path to save the figure.
    """
    
    # Configurations
    group_trait = 'Targeted_gene_name'
    legend_row = 4
    direction = 'up'
    
    def prepare_data(df, direction, trait):
        """Prepare the data by sorting and setting categorical order."""
        df_sorted = df.sort_values(by='Numbered_gene_name', ascending=(direction == 'up')).copy()
        
        # Set categorical order for the group and legend traits
        df_sorted[group_trait] = pd.Categorical(df_sorted[group_trait], ordered=True, categories=df_sorted[group_trait].unique())
        
        # Generate new axis values for the plot
        gene_counts = df_sorted.groupby('Numbered_gene_name')[trait].count()
        gene_counts = gene_counts.reindex(df_sorted['Numbered_gene_name'].unique(), fill_value=0)
        
        new_axis = []
        simple_axis = []
        for i, count in enumerate(gene_counts):
            x = i + 1
            temp_list = np.linspace(x - 0.5, x + 0.5, count)
            new_axis.extend(temp_list)
            simple_axis.append(x)
        
        df_sorted['New_axis'] = new_axis
        
        return df_sorted, simple_axis, gene_counts.index

    # Setup the figure and grid
    fig = plt.figure(figsize=fig_size)
    gs = gridspec.GridSpec(22, 12)
    
    # Prepare and plot data for the first subplot
    df1_prepared, x_ticks1, gene_order1 = prepare_data(df1, direction, trait)
    ax1 = fig.add_subplot(gs[:5, 0:10])
    plot_trait_panel(ax1, df1_prepared, trait, title1, y_label, x_ticks1, gene_order1, logy, ymin, ymax)
    
    # Prepare and plot data for the second subplot
    if df2 is not None:
        df2_prepared, x_ticks2, gene_order2 = prepare_data(df2, direction, trait)
        ax2 = fig.add_subplot(gs[6:11, 0:10])
        plot_trait_panel(ax2, df2_prepared, trait, title2, y_label, x_ticks2, gene_order2, logy, ymin, ymax)
    
    # Prepare and plot data for the third subplot
    if df3 is not None:
        df3_prepared, x_ticks3, gene_order3 = prepare_data(df3, direction, trait)
        ax3 = fig.add_subplot(gs[12:17, 0:10])
        plot_trait_panel(ax3, df3_prepared, trait, title3, y_label, x_ticks3, gene_order3, logy, ymin, ymax)
    
    if df4 is not None:
        df4_prepared, x_ticks4, gene_order4 = prepare_data(df4, direction, trait)
        ax4 = fig.add_subplot(gs[18:, 0:10])
        plot_trait_panel(ax4, df4_prepared, trait, title4, y_label, x_ticks4, gene_order4, logy, ymin, ymax)
    
    # Save the figure if required
    if fig_name and output_path:
        fig.savefig(f'{output_path}/{fig_name}', bbox_inches='tight')

def plot_trait_panel(ax, df, trait, title, y_label, x_ticks, gene_order, logy, ymin, ymax):
    """
    Helper function to plot a single panel with error bars and scatter plot.
    Parameters:
    - ax: The axis to plot on.
    - df: The DataFrame containing the data.
    - trait: The column to be plotted on the y-axis.
    - title: The title of the plot.
    - y_label: The label for the y-axis.
    - x_ticks: The positions of the x-ticks.
    - gene_order: The order of the genes along the x-axis.
    - logy: Boolean. If True, use a logarithmic scale for the y-axis.
    - y_min: Optional. Minimum value for the y-axis.
    - y_max: Optional. Maximum value for the y-axis.
    """
    # Plot error bars
    error_bar = [df[trait] - df[f'{trait}_2.5P'], df[f'{trait}_97.5P'] - df[trait]]
    contains_neg = (error_bar[0] < 0).any() or (error_bar[1] < 0).any()
    if contains_neg: # if the observed value error bar has negative value
        print(f'observed value error bar contains negative value')
        error_bar_median = [df[trait+'_bootstrap_median'] - df[f'{trait}_2.5P'], df[f'{trait}_97.5P'] - df[trait+'_bootstrap_median']]
        contains_neg_in_error_bar_median = (error_bar_median[0] < 0).any() or (error_bar_median[1] < 0).any()
        if contains_neg_in_error_bar_median: # if bs median error has negative value use bootstrap mean
            print(f'bootstrap mean error bar is used')
            error_bar_mean = [df[trait+'_bootstrap_mean'] - df[f'{trait}_2.5P'], df[f'{trait}_97.5P'] - df[trait+'_bootstrap_mean']]
            ax.errorbar(df['New_axis'], df[trait], 
            yerr=error_bar_mean, 
            linestyle='', color='black', elinewidth=0.5, capsize=1.5, capthick=0.5)
        else:   # if bs mean error bar has negative value use bootstrap median
            print(f'bootstrap median error bar is used')
            ax.errorbar(df['New_axis'], df[trait], 
            yerr=error_bar_median, 
            linestyle='', color='black', elinewidth=0.5, capsize=1.5, capthick=0.5)
    else:
        ax.errorbar(df['New_axis'], df[trait], 
                yerr=error_bar, 
                linestyle='', color='black', elinewidth=0.5, capsize=1.5, capthick=0.5)
    
    # Plot scatter plot
    sns.scatterplot(x='New_axis', y=trait, hue='Targeted_gene_name', data=df, ax=ax, s=25)
    
    # Add a control line at y=1
    if 'log' in trait or trait == 'ScoreRTN':
        ax.axhline(y=0, color='black', linestyle='--')
    else:
        ax.axhline(y=1, color='black', linestyle='--')
    
    # Draw vertical lines to separate Targeted_gene_name groups
    unique_genes = df['Targeted_gene_name'].unique()
    for i in range(1, len(unique_genes)):
        last_pos = df[df['Targeted_gene_name'] == unique_genes[i-1]]['New_axis'].max()
        ax.axvline(x=last_pos + 0.5, color='lightgrey', linestyle='-', linewidth=0.8)
    
    # Set titles and labels
    ax.set_title(title, loc='left')
    ax.set_ylabel(y_label)
    ax.set_xlabel(None)
    
    # Rotate x-axis labels for readability
    ax.tick_params(axis='x', labelrotation=90)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(gene_order)
    
    # Set y-axis to log scale if required
    # important to set the scale before defining the limit
    # This ensures that the limits are interpreted correctly in log scale
    if logy:
        ax.set_yscale('log', base=2)
        ax.set_yticklabels(ax.get_yticks())
    
    # Set y-axis limits if provided
    if ymin is not None or ymax is not None:
        ax.set_ylim(ymin, ymax)
        
    # Hide the legend
    ax.legend([], [], frameon=False)
