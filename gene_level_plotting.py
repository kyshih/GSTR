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


