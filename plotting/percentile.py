""" This module contains functions to plot percentile plot.
    TODO: refactor the percentile plotting code
    TODO: set ylim doesn't work prob cuz of tick setting messes up with it
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec

def plot_percentile_comparisons(input_df,input_df1=None,input_df2=None,title1='',title2='',title3='',
                                input_size=(12, 17),logy=True,ymax=None, ymin=None, grid=False, xlabel_fontsize=10,perc_legend=False,
                                fig_name=None,output_path=None):
    # KTHC df, KT df and the trait I am focusing on
    trait_of_interest = 'LN_mean_relative'
    order_by_trait_of_interest = False

    # percentile I am interested in
    input_percentile = [50, 70, 80, 90, 95]

    temp1 = [str(x) for x in input_percentile]

    temp = np.array_split((np.array(range(1,101))/100)[::-1],len(temp1))
    temp2 = [x[0] for x in temp][::-1]
    opacity_dic = dict(zip(temp1,temp2))

    # group_trait is the sgRNA first group, # legend group. trait is the trait the dot's color are based on 
    # I want to make sure legend_group to be the same as group_trait
    group_trait, legend_group_trait, direction = 'Targeted_gene_name', 'Targeted_gene_name', 'up',

    # output address and if y axis is log scaled
    # (width, height)
    lengend_row = 2

    # y label name
    temp_y_label_name = 'Tumor size at xth percentile \n'

    temp_x_label_name = 'Target gene'

    # grid specification
    gs = gridspec.GridSpec(17, 5)

    # set figure size
    fig1 = plt.figure(figsize=input_size)

    # set the first panel axis
    ax1=fig1.add_subplot(gs[:5, 0:5])

    # sort based on focal trait value first
    if direction == 'up':
        #temp_df = input_df.sort_values(by = trait_of_interest,ascending=False).copy()
        temp_df = input_df.sort_values(by = 'Numbered_gene_name',ascending=True).copy()
    if direction == 'down':
        temp_df = input_df.sort_values(by = trait_of_interest,ascending=True).copy()

    # record the categorical order of legend_group_trait. 
    # I do this such that in each legend group, the sgRNA is order based on their focal trait value
    order = temp_df[legend_group_trait].unique()
    temp_df[legend_group_trait] = pd.CategoricalIndex(temp_df[legend_group_trait], ordered=True, categories=order)

    # record the categorical order of group trait
    order2 = temp_df[group_trait].unique()
    temp_df[group_trait] = pd.CategoricalIndex(temp_df[group_trait], ordered=True, categories=order2)    

    # reorder the dataframe.
    if order_by_trait_of_interest:
        if direction == 'up':
            temp_df = temp_df.sort_values(by = [legend_group_trait,group_trait,trait_of_interest],ascending =[True,True,False])
        if direction == 'down':
            temp_df = temp_df.sort_values(by = [legend_group_trait,group_trait,trait_of_interest],ascending =[True,True,True])
        
    # Melt dataframe
    temp_q_list0 = [str(x) + '_percentile_relative' for x in input_percentile]
    temp_q_list1 = [x +'_2.5P' for x in temp_q_list0]
    temp_q_list2 = [x +'_97.5P' for x in temp_q_list0]
    temp_sub = pd.melt(temp_df , id_vars=['gRNA','Targeted_gene_name','Numbered_gene_name'],
                    value_vars=temp_q_list0+temp_q_list1+temp_q_list2)
    temp_sub['X_label'] = temp_sub['Numbered_gene_name']+'_'+temp_sub['variable']
    order_x = temp_sub['Numbered_gene_name'].unique()
    temp_sub['Numbered_gene_name'] = pd.CategoricalIndex(temp_sub['Numbered_gene_name'], ordered=True, categories=order_x)
    temp_sub = temp_sub.sort_values(by = ['Numbered_gene_name','X_label'])
    temp_sub['Percentile_label'] = temp_sub['variable'].apply(lambda x: x.split('_percentile')[0])
    temp_sub0 = temp_sub[temp_sub['variable'].isin(temp_q_list0)].copy()
    temp_sub1 = temp_sub[temp_sub['variable'].isin(temp_q_list1)].copy()
    temp_sub2 = temp_sub[temp_sub['variable'].isin(temp_q_list2)].copy()

    # create new axis
    temp_g = temp_sub0.groupby(['Numbered_gene_name'],as_index=False)['value'].count()
    order3 = list(dict.fromkeys(temp_sub0['Numbered_gene_name'].values))
    # reorder the gene order in temp_g to match with the temp_df
    temp_g['Numbered_gene_name'] = pd.CategoricalIndex(temp_g['Numbered_gene_name'], ordered=True, categories=order3) 
    temp_g = temp_g.sort_values(by =['Numbered_gene_name'] ,ascending =[True])

    # make new axis
    temp_new_axis = []
    simple_axis = []
    for n,sgRNA_count in enumerate(temp_g['value'].values):
        x = n+1
        temp_list = generate_symmetric_list(sgRNA_count)+x
        temp_new_axis = temp_new_axis + list(temp_list)
        simple_axis.append(x)
    temp_sub0['New_axis'] = temp_new_axis

    # generate color_opacity
    temp_op = [opacity_dic.get(x) for x in temp_sub0['Percentile_label']]

    # Plot
    ax1.errorbar(temp_sub0['New_axis'], temp_sub0['value'], 
                     yerr = [np.array(temp_sub0['value'].to_list())-np.array(temp_sub1['value'].to_list()),
                            np.array(temp_sub2['value'].to_list())-np.array(temp_sub0['value'].to_list())],
                    elinewidth=0.75, linestyle='', c = 'black')   

    # sns.scatterplot(x='New_axis', y='value', hue=legend_group_trait,
    #               data=temp_sub0, alpha = temp_op, ax = ax1)
    # sns.scatterplot(x='New_axis', y='value', hue=legend_group_trait,
    #               data=temp_sub0, alpha = 1, ax = ax1) 

    opacity_dic = dict(zip(temp_sub0['Percentile_label'], temp_op))
    # Loop through each unique Percentile_label
    for percentile in temp_sub0['Percentile_label'].unique():
        # Subset the data for the current percentile
        subset = temp_sub0[temp_sub0['Percentile_label'] == percentile]
        
        # Get the corresponding alpha value from your dictionary (or set a default)
        alpha_value = opacity_dic.get(str(percentile), 1.0)  # Default alpha to 1.0 if not found

        # Plot the subset with the specific alpha value
        sns.scatterplot(
            x='New_axis', 
            y='value', 
            hue='Targeted_gene_name', 
            data=subset, 
            alpha=alpha_value, 
            ax=ax1, 
            legend=False  # We'll handle the legend separately
        ) 
        
    # Calculate midpoints between groups of Targeted_gene_name
    targeted_gene_name_bound = temp_sub0.groupby(['Targeted_gene_name'], as_index=False).agg(bound=('New_axis',max))['bound']+0.3

    # Add vertical lines at the midpoints to separate the Targeted_gene_name groups
    for elem in targeted_gene_name_bound:
        ax1.axvline(x=elem, color='lightgray', linestyle='--', linewidth=1)

    # manipulate axes tick and label
    ax1.set_xticklabels(temp_g['Numbered_gene_name'])
    ax1.set_xticks(simple_axis)
    ax1.axhline(y=1,color='black', linestyle='--')
    # # log scale y axis
    # ax1.set_ylim(temp_min_y, temp_max_y)
        
    if logy == True:
        ax1.set_yscale('log',base = 2)
        #ax1.set_yscale('log',base = 10)
        temp_ticks_list=ax1.get_yticks()
        ax1.set_yticklabels(temp_ticks_list)
        
    # Set y-axis limits if provided
    if ymin is not None or ymax is not None:
        ax1.set_ylim(ymin, ymax)
     
    # set x,y axis label
    ax1.set_xlabel(None)
    ax1.set_ylabel(temp_y_label_name)
    # Rotate x ticklabel orientation
    ax1.tick_params(axis='x', labelrotation = 90)
    # change front size
    if xlabel_fontsize:
        ax1.tick_params(axis='x', which='major', labelsize=xlabel_fontsize)

    # change legend row number
    # temp_col = round(len(temp_sub['Percentile_label'].unique())/lengend_row) # column number to make sure the row for legend is 3
    temp_col = round(len(temp_sub['Targeted_gene_name'].unique())/lengend_row)
    if grid:
        plt.grid(True)
    ax1.legend([],[], frameon=False) 
    ax1.set_title(title1, loc ='Left')
    # ax1.legend(ncol=temp_col, loc='upper right', title = 'Gene')
    
    if input_df1 is not None:
        # set the first panel axis
        ax2=fig1.add_subplot(gs[6:11, 0:5])

        # sort based on focal trait value first
        if direction == 'up':
            #temp_df = input_df.sort_values(by = trait_of_interest,ascending=False).copy()
            temp_df = input_df1.sort_values(by = 'Numbered_gene_name',ascending=True).copy()
        if direction == 'down':
            temp_df = input_df1.sort_values(by = trait_of_interest,ascending=True).copy()

        # record the categorical order of legend_group_trait. 
        # I do this such that in each legend group, the sgRNA is order based on their focal trait value
        order = temp_df[legend_group_trait].unique()
        temp_df[legend_group_trait] = pd.CategoricalIndex(temp_df[legend_group_trait], ordered=True, categories=order)

        # record the categorical order of group trait
        order2 = temp_df[group_trait].unique()
        temp_df[group_trait] = pd.CategoricalIndex(temp_df[group_trait], ordered=True, categories=order2)    

        # reorder the dataframe.
        if order_by_trait_of_interest:
            if direction == 'up':
                temp_df = temp_df.sort_values(by = [legend_group_trait,group_trait,trait_of_interest],ascending =[True,True,False])
            if direction == 'down':
                temp_df = temp_df.sort_values(by = [legend_group_trait,group_trait,trait_of_interest],ascending =[True,True,True])
            
        # Melt dataframe
        temp_q_list0 = [str(x) + '_percentile_relative' for x in input_percentile]
        temp_q_list1 = [x +'_2.5P' for x in temp_q_list0]
        temp_q_list2 = [x +'_97.5P' for x in temp_q_list0]
        temp_sub = pd.melt(temp_df , id_vars=['gRNA','Targeted_gene_name','Numbered_gene_name'],
                        value_vars=temp_q_list0+temp_q_list1+temp_q_list2)
        temp_sub['X_label'] = temp_sub['Numbered_gene_name']+'_'+temp_sub['variable']
        order_x = temp_sub['Numbered_gene_name'].unique()
        temp_sub['Numbered_gene_name'] = pd.CategoricalIndex(temp_sub['Numbered_gene_name'], ordered=True, categories=order_x)
        temp_sub = temp_sub.sort_values(by = ['Numbered_gene_name','X_label'])
        temp_sub['Percentile_label'] = temp_sub['variable'].apply(lambda x: x.split('_percentile')[0])
        temp_sub0 = temp_sub[temp_sub['variable'].isin(temp_q_list0)].copy()
        temp_sub1 = temp_sub[temp_sub['variable'].isin(temp_q_list1)].copy()
        temp_sub2 = temp_sub[temp_sub['variable'].isin(temp_q_list2)].copy()

        # create new axis
        temp_g = temp_sub0.groupby(['Numbered_gene_name'],as_index=False)['value'].count()
        order3 = list(dict.fromkeys(temp_sub0['Numbered_gene_name'].values))
        # reorder the gene order in temp_g to match with the temp_df
        temp_g['Numbered_gene_name'] = pd.CategoricalIndex(temp_g['Numbered_gene_name'], ordered=True, categories=order3) 
        temp_g = temp_g.sort_values(by =['Numbered_gene_name'] ,ascending =[True])

        # make new axis
        temp_new_axis = []
        simple_axis = []
        for n,sgRNA_count in enumerate(temp_g['value'].values):
            x = n+1
            temp_list = generate_symmetric_list(sgRNA_count)+x
            temp_new_axis = temp_new_axis + list(temp_list)
            simple_axis.append(x)
        temp_sub0['New_axis'] = temp_new_axis

        # generate color_opacity
        temp_op = [opacity_dic.get(x) for x in temp_sub0['Percentile_label']]

        # Plot
        ax2.errorbar(temp_sub0['New_axis'], temp_sub0['value'], 
                    yerr = [np.array(temp_sub0['value'].to_list())-np.array(temp_sub1['value'].to_list()),
                            np.array(temp_sub2['value'].to_list())-np.array(temp_sub0['value'].to_list())],
                    elinewidth=0.75, linestyle='', c = 'black')

        # sns.scatterplot(x='New_axis', y='value', hue=legend_group_trait,
        #               data=temp_sub0, alpha = temp_op, ax = ax2)
        # sns.scatterplot(x='New_axis', y='value', hue=legend_group_trait,
        #               data=temp_sub0, alpha = 1, ax = ax2)  

        opacity_dic = dict(zip(temp_sub0['Percentile_label'], temp_op))
        # Loop through each unique Percentile_label
        for percentile in temp_sub0['Percentile_label'].unique():
            # Subset the data for the current percentile
            subset = temp_sub0[temp_sub0['Percentile_label'] == percentile]
            
            # Get the corresponding alpha value from your dictionary (or set a default)
            alpha_value = opacity_dic.get(str(percentile), 1.0)  # Default alpha to 1.0 if not found

            # Plot the subset with the specific alpha value
            sns.scatterplot(
                x='New_axis', 
                y='value', 
                hue='Targeted_gene_name', 
                data=subset, 
                alpha=alpha_value, 
                ax=ax2, 
                legend=False  # We'll handle the legend separately
            )
        
        # Calculate midpoints between groups of Targeted_gene_name
        targeted_gene_name_bound = temp_sub0.groupby(['Targeted_gene_name'], as_index=False).agg(bound=('New_axis',max))['bound']+0.3

        # Add vertical lines at the midpoints to separate the Targeted_gene_name groups
        for elem in targeted_gene_name_bound:
            ax2.axvline(x=elem, color='lightgray', linestyle='--', linewidth=1)

        # manipulate axes tick and label
        ax2.set_xticklabels(temp_g['Numbered_gene_name'])
        ax2.set_xticks(simple_axis)
        ax2.axhline(y=1,color='black', linestyle='--')
        # # log scale y axis
        # ax2.set_ylim(temp_min_y, temp_max_y)
            
        if logy == True:
            ax2.set_yscale('log',base = 2)
            #ax2.set_yscale('log',base = 10)
            temp_ticks_list=ax2.get_yticks()
            ax2.set_yticklabels(temp_ticks_list)
            
        # Set y-axis limits if provided
        if ymin is not None or ymax is not None:
            ax2.set_ylim(ymin, ymax)
            
        # set x,y axis label
        ax2.set_xlabel(None)
        ax2.set_ylabel(temp_y_label_name)
        # Rotate x ticklabel orientation
        ax2.tick_params(axis='x', labelrotation = 90)
        # change front size
        if xlabel_fontsize:
            ax2.tick_params(axis='x', which='major', labelsize=xlabel_fontsize)

        # change legend row number
        # temp_col = round(len(temp_sub['Percentile_label'].unique())/lengend_row) # column number to make sure the row for legend is 3
        temp_col = round(len(temp_sub['Targeted_gene_name'].unique())/lengend_row)
        ax2.legend([],[], frameon=False) 
        ax2.set_title(title2, loc ='Left')
    
    if input_df2 is not None:
        # set the first panel axis
        ax3=fig1.add_subplot(gs[12:17, 0:5])

        # sort based on focal trait value first
        if direction == 'up':
            #temp_df = input_df.sort_values(by = trait_of_interest,ascending=False).copy()
            temp_df = input_df2.sort_values(by = 'Numbered_gene_name',ascending=True).copy()
        if direction == 'down':
            temp_df = input_df2.sort_values(by = trait_of_interest,ascending=True).copy()

        # record the categorical order of legend_group_trait. 
        # I do this such that in each legend group, the sgRNA is order based on their focal trait value
        order = temp_df[legend_group_trait].unique()
        temp_df[legend_group_trait] = pd.CategoricalIndex(temp_df[legend_group_trait], ordered=True, categories=order)

        # record the categorical order of group trait
        order2 = temp_df[group_trait].unique()
        temp_df[group_trait] = pd.CategoricalIndex(temp_df[group_trait], ordered=True, categories=order2)    

        # reorder the dataframe.
        if order_by_trait_of_interest:
            if direction == 'up':
                temp_df = temp_df.sort_values(by = [legend_group_trait,group_trait,trait_of_interest],ascending =[True,True,False])
            if direction == 'down':
                temp_df = temp_df.sort_values(by = [legend_group_trait,group_trait,trait_of_interest],ascending =[True,True,True])
            
        # Melt dataframe
        temp_q_list0 = [str(x) + '_percentile_relative' for x in input_percentile]
        temp_q_list1 = [x +'_2.5P' for x in temp_q_list0]
        temp_q_list2 = [x +'_97.5P' for x in temp_q_list0]
        temp_sub = pd.melt(temp_df , id_vars=['gRNA','Targeted_gene_name','Numbered_gene_name'],
                        value_vars=temp_q_list0+temp_q_list1+temp_q_list2)
        temp_sub['X_label'] = temp_sub['Numbered_gene_name']+'_'+temp_sub['variable']
        order_x = temp_sub['Numbered_gene_name'].unique()
        temp_sub['Numbered_gene_name'] = pd.CategoricalIndex(temp_sub['Numbered_gene_name'], ordered=True, categories=order_x)
        temp_sub = temp_sub.sort_values(by = ['Numbered_gene_name','X_label'])
        temp_sub['Percentile_label'] = temp_sub['variable'].apply(lambda x: x.split('_percentile')[0])
        temp_sub0 = temp_sub[temp_sub['variable'].isin(temp_q_list0)].copy()
        temp_sub1 = temp_sub[temp_sub['variable'].isin(temp_q_list1)].copy()
        temp_sub2 = temp_sub[temp_sub['variable'].isin(temp_q_list2)].copy()

        # create new axis
        temp_g = temp_sub0.groupby(['Numbered_gene_name'],as_index=False)['value'].count()
        order3 = list(dict.fromkeys(temp_sub0['Numbered_gene_name'].values))
        # reorder the gene order in temp_g to match with the temp_df
        temp_g['Numbered_gene_name'] = pd.CategoricalIndex(temp_g['Numbered_gene_name'], ordered=True, categories=order3) 
        temp_g = temp_g.sort_values(by =['Numbered_gene_name'] ,ascending =[True])

        # make new axis
        temp_new_axis = []
        simple_axis = []
        for n,sgRNA_count in enumerate(temp_g['value'].values):
            x = n+1
            temp_list = generate_symmetric_list(sgRNA_count)+x
            temp_new_axis = temp_new_axis + list(temp_list)
            simple_axis.append(x)
        temp_sub0['New_axis'] = temp_new_axis

        # generate color_opacity
        temp_op = [opacity_dic.get(x) for x in temp_sub0['Percentile_label']]

        # Plot
        ax3.errorbar(temp_sub0['New_axis'], temp_sub0['value'], 
                    yerr = [np.array(temp_sub0['value'].to_list())-np.array(temp_sub1['value'].to_list()),
                            np.array(temp_sub2['value'].to_list())-np.array(temp_sub0['value'].to_list())],
                    elinewidth=0.75, linestyle='', c = 'black')

        # sns.scatterplot(x='New_axis', y='value', hue=legend_group_trait,
        #               data=temp_sub0, alpha = temp_op, ax = ax3)
        # sns.scatterplot(x='New_axis', y='value', hue=legend_group_trait,
        #               data=temp_sub0, alpha = 1, ax = ax3)  

        opacity_dic = dict(zip(temp_sub0['Percentile_label'], temp_op))
        # Loop through each unique Percentile_label
        for percentile in temp_sub0['Percentile_label'].unique():
            # Subset the data for the current percentile
            subset = temp_sub0[temp_sub0['Percentile_label'] == percentile]
            
            # Get the corresponding alpha value from your dictionary (or set a default)
            alpha_value = opacity_dic.get(str(percentile), 1.0)  # Default alpha to 1.0 if not found

            # Plot the subset with the specific alpha value
            sns.scatterplot(
                x='New_axis', 
                y='value', 
                hue='Targeted_gene_name', 
                data=subset, 
                alpha=alpha_value, 
                ax=ax3, 
                legend=False  # We'll handle the legend separately
            )
                
        # Calculate midpoints between groups of Targeted_gene_name
        targeted_gene_name_bound = temp_sub0.groupby(['Targeted_gene_name'], as_index=False).agg(bound=('New_axis',max))['bound']+0.3

        # Add vertical lines at the midpoints to separate the Targeted_gene_name groups
        for elem in targeted_gene_name_bound:
            ax3.axvline(x=elem, color='lightgray', linestyle='--', linewidth=1)
            
        # manipulate axes tick and label
        ax3.set_xticklabels(temp_g['Numbered_gene_name'])
        ax3.set_xticks(simple_axis)
        ax3.axhline(y=1,color='black', linestyle='--')
        # # log scale y axis
        # ax3.set_ylim(temp_min_y, temp_max_y)
            
        if logy == True:
            ax3.set_yscale('log',base = 2)
            #ax3.set_yscale('log',base = 10)
            temp_ticks_list=ax3.get_yticks()
            ax3.set_yticklabels(temp_ticks_list)
            
        # Set y-axis limits if provided
        if ymin is not None or ymax is not None:
            ax3.set_ylim(ymin, ymax)
            
        # set x,y axis label
        ax3.set_xlabel(None)
        ax3.set_ylabel(temp_y_label_name)
        # Rotate x ticklabel orientation
        ax3.tick_params(axis='x', labelrotation = 90)
        # change front size
        if xlabel_fontsize:
            ax3.tick_params(axis='x', which='major', labelsize=xlabel_fontsize)

        # change legend row number
        # temp_col = round(len(temp_sub['Percentile_label'].unique())/lengend_row) # column number to make sure the row for legend is 3
        temp_col = round(len(temp_sub['Targeted_gene_name'].unique())/lengend_row)
        ax3.legend([],[], frameon=False) 
        ax3.set_title(title3, loc ='Left')

        # Save the figure if required
    
    # Call the function to create the legend
    if perc_legend:
        legend_labels = [50, 70, 80, 90, 95]
        create_greyscale_dot_legend(ax1, len(legend_labels), legend_labels)
    
    if fig_name and output_path:
        fig1.savefig(f'{output_path}/{fig_name}', bbox_inches='tight')
        
def generate_symmetric_list(input_n):
    # generate a symmetric list of numbers centered around 0, used to distribute points symmetrically around an axis when plotting
    if input_n % 2 == 1:
        return(np.array(range(-int((input_n-1)/2),int((input_n+1)/2)))*0.1)
    else:
        return(np.array(range(-int(input_n/2),int(input_n/2)))*0.1+0.05)
        
# this gives vertical dots
def create_greyscale_dot_legend(ax, num_dots, legend_labels, dot_size=100):
    handles = []
    colors = plt.cm.gray_r(np.linspace(0.2, 0.7, num_dots))  # Adjust the range as needed

    for i, color in enumerate(colors):
        #handle = ax.scatter([], [], c=color, s=dot_size, marker='o')
        handle = ax.scatter(np.arange(num_dots), np.zeros(num_dots), c=color, s=dot_size, marker='o')
        handles.append(handle)

    ax.legend(handles, legend_labels, title="Percentiles", loc='upper right')

# this create horizontal dots
def create_greyscale_horitontal_dot_legend(num_dots, legend_labels, dot_size=100, width=2, height=1, title='Percentiles',
                                fig_name=None,output_path=None):
    fig, ax = plt.subplots(figsize=(width, height))
    fig.suptitle(title, fontsize=14, y=1.05)  # Adding a title above the legend
    ax.set_axis_off()

    colors = plt.cm.gray_r(np.linspace(0.2, 0.7, num_dots))  # Adjust the range as needed

    for i, color in enumerate(colors):
        ax.scatter(i / (num_dots - 4), 0.2, c=color, s=dot_size, marker='o')

    for i, label in enumerate(legend_labels):
        ax.text(i / (len(legend_labels) - 0.9), -0.1, label, fontsize=8, transform=ax.transAxes, ha='center', va='center')

    plt.show()
    if fig_name and output_path:
        fig.savefig(f'{output_path}/{fig_name}', bbox_inches='tight')
    
def generate_percentile_opacity(input_percentile):
    """Generate a dictionary mapping percentiles to opacity values."""
    temp = np.array_split((np.array(range(1, 101)) / 100)[::-1], len(input_percentile)) # split a [0.01, 1] in reverse into chunks. each chunk corresponds to a percentile
    return dict(zip([str(x) for x in input_percentile], [x[0] for x in temp][::-1])) # extract the first value of each chunk and assigns it as an opacity

def melt_dataframe(df, percentile_labels):
    """Melt the dataframe to prepare it for plotting."""
    temp_q_list0 = [f"{x}_percentile_relative" for x in percentile_labels]
    temp_q_list1 = [f"{x}_2.5P" for x in temp_q_list0]
    temp_q_list2 = [f"{x}_97.5P" for x in temp_q_list0]
    
    temp_sub = pd.melt(df, id_vars=['gRNA', 'Targeted_gene_name', 'Numbered_gene_name'],
                       value_vars=temp_q_list0 + temp_q_list1 + temp_q_list2)
    temp_sub['X_label'] = temp_sub['Numbered_gene_name'].astype(str) + '_' + temp_sub['variable']
    temp_sub['Percentile_label'] = temp_sub['variable'].apply(lambda x: x.split('_percentile')[0])
    
    temp_sub0 = temp_sub[temp_sub['variable'].isin(temp_q_list0)].copy()
    temp_sub1 = temp_sub[temp_sub['variable'].isin(temp_q_list1)].copy()
    temp_sub2 = temp_sub[temp_sub['variable'].isin(temp_q_list2)].copy()
    
    return temp_sub0, temp_sub1, temp_sub2

def create_new_axis(temp_g, value_column):
    """Create new axis positions for the plot."""
    temp_new_axis = []
    simple_axis = []
    for n, sgRNA_count in enumerate(temp_g[value_column].values):
        x = n + 1
        temp_list = generate_symmetric_list(sgRNA_count) + x
        temp_new_axis.extend(temp_list)
        simple_axis.append(x)
    return temp_new_axis, simple_axis

def add_legend(ax, handles, nrow):
    """Add a legend to the plot."""
    ax.legend(handles=handles, ncol=nrow, loc='upper right', title='Gene')