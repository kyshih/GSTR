import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_trait_comparison_gene_level(df1, df2, df3, title1, trait, ylabel, title2='', title3='', big_title='',logy=True,
                                     fig_name=None, output_path=None): 
    fig = plt.figure(figsize=(12,10))
    gs = fig.add_gridspec(3,1)
    fig.suptitle(big_title)

    ax1 = fig.add_subplot(gs[0,0])
    plot_trait_panel_barplot(ax1, df1, trait, title1, ylabel, logy)
    if df2 is not None:
        ax2 = fig.add_subplot(gs[1,0])
        plot_trait_panel_barplot(ax2, df2, trait, title2, ylabel, logy)
    if df3 is not None:
        ax3 = fig.add_subplot(gs[2,0])
        plot_trait_panel_barplot(ax3, df3, trait, title3, ylabel, logy)
        
    # Adjust the space between the subplots
    plt.subplots_adjust(hspace=0.5)  # Add more space between the subplots
    
    # Save the figure if required
    if fig_name and output_path:
        fig.savefig(f'{output_path}/{fig_name}', bbox_inches='tight')
    
def plot_trait_panel_barplot(ax, df, trait, title, ylabel, logy):
    sns.barplot(data=df, y=trait, x='Targeted_gene_name', ax=ax, hue='Targeted_gene_name')
    error_bar = [df[trait] - df[f'{trait}_2.5P'], df[f'{trait}_97.5P'] - df[trait]]
    ax.errorbar(df['Targeted_gene_name'], df[trait], 
                yerr=error_bar, 
                linestyle='', color='black', elinewidth=0.5, capsize=1.5, capthick=0.5)
    
    if trait == 'ScoreRTN':
        ax.axhline(y=0, linestyle='--', color='black')
    else:
        ax.axhline(y=1, linestyle='--', color='black')
    ax.set_xlabel(None)
    ax.set_ylabel(ylabel)
    ax.set_title(title,loc='left',fontsize=10)
    # Rotate x-axis labels for readability
    ax.tick_params(axis='x', labelrotation=90)
    if logy:
        ax.set_yscale('log', base=2)
        
        # Find the range considering the error bars
        min_value = np.min(df[f'{trait}_2.5P'])
        max_value = np.max(df[f'{trait}_97.5P'])
        
        # Set y-ticks dynamically based on the range of the data and error bars
        min_tick = np.floor(np.log2(min_value))
        max_tick = np.ceil(np.log2(max_value))
        ticks = [2**i for i in range(int(min_tick), int(max_tick) + 1)]
        
        # Set y-ticks and their labels to 2, 4, 8, 16, etc.
        ax.set_yticks(ticks)
        ax.set_yticklabels([str(tick) for tick in ticks])  # Set y-tick labels as 2, 4, 8, etc

def plot_trait_comparison_gene_level_zero_baseline(df1, df2, df3, title1, trait, ylabel, title2='', title3='', big_title='',hue_group='Targeted_gene_name',
                                                   set_legend='auto',logy=True,fig_name=None, output_path=None): 
    fig = plt.figure(figsize=(12,10))
    gs = fig.add_gridspec(3,1)
    fig.suptitle(big_title)

    ax1 = fig.add_subplot(gs[0,0])
    plot_trait_panel_barplot_zero_baseline(ax1, df1, trait, title1, ylabel, logy, hue_group, set_legend)
    if df2 is not None:
        ax2 = fig.add_subplot(gs[1,0])
        plot_trait_panel_barplot_zero_baseline(ax2, df2, trait, title2, ylabel, logy, hue_group, set_legend)
    if df3 is not None:
        ax3 = fig.add_subplot(gs[2,0])
        plot_trait_panel_barplot_zero_baseline(ax3, df3, trait, title3, ylabel, logy, hue_group, set_legend)
        
    # Adjust the space between the subplots
    plt.subplots_adjust(hspace=0.5)  # Add more space between the subplots
    
    # Save the figure if required
    if fig_name and output_path:
        fig.savefig(f'{output_path}/{fig_name}', bbox_inches='tight')

def plot_trait_panel_barplot_zero_baseline(ax, df, trait, title, ylabel, logy, hue_group='Targeted_gene_name', set_legend='auto'):
    # Make sure all values are centered around 0
    df = df.copy()
    df['Trait_value'] = df[trait] - 0  # Normalize so 1 is the baseline
    
    # Bar plot with zero-centered bars
    sns.barplot(data=df, y='Trait_value', x='Targeted_gene_name', ax=ax, hue=hue_group, legend=set_legend)
    
    # Error bars
    error_bar = [df['Trait_value'] - (df[f'{trait}_2.5P'] - 0), 
                 (df[f'{trait}_97.5P'] - 0) - df['Trait_value']]
    
    ax.errorbar(df['Targeted_gene_name'], df['Trait_value'], 
                yerr=error_bar, 
                linestyle='', color='black', elinewidth=0.5, capsize=1.5, capthick=0.5)

    # Zero reference line
    ax.axhline(y=0, linestyle='--', color='black')

    # Labels and formatting
    ax.set_xlabel(None)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left', fontsize=10)
    ax.tick_params(axis='x', labelrotation=90)
    
    # Log scale with symmetric positive and negative ticks
    if logy:
        min_value = np.min(df[f'{trait}_2.5P'])
        max_value = np.max(df[f'{trait}_97.5P'])

        # Compute dynamic tick range
        min_tick = np.floor(np.log2(abs(min_value))) if min_value < 0 else 0
        max_tick = np.ceil(np.log2(max_value)) if max_value > 0 else 0
        ticks = sorted([-2**i for i in range(int(abs(min_tick)), -1, -1)] + 
                       [2**i for i in range(0, int(max_tick) + 1)])
        
        # Apply log scale if valid
        if min_tick < 0 or max_tick > 0:
            ax.set_yscale('symlog', linthresh=1, base=2)  # Symmetric log scale
            ax.set_yticks(ticks)
            ax.set_yticklabels([f'{tick:.0f}' for tick in ticks])

    # Remove extra white space
    ax.margins(y=0.1)