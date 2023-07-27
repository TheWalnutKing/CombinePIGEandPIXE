import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy
import os
import natsort as natsorted
import collections
from collections import Counter
import sys
import statistics
from scipy import optimize


# region Importing PIGE and PIXE Data



your_folder_path = input('What is the path to your folder with PIGE and PIXE in it? ')

con = float(input('Enter concentration as decimal please '))

triplicates = input('Are your samples triplicates? Yes or No ')

h_needed = input('Do you want to calculate h-factor from standards? Yes or No ')


your_folder_contains = os.listdir(your_folder_path)



your_PIXE_file = [x for x in your_folder_contains if 'PIXECON.CSV' in x]

if your_PIXE_file == []:
    sys.exit("Your PIXE file doesn't exist or is not formatted correctly (should include PIXECON.CSV)")


your_path_PIXE = your_folder_path + '\\' + your_PIXE_file[0]

your_PIGE_file = [x for x in your_folder_contains if 'PIGEaucRepv3' in x]

if your_PIGE_file == []:
    sys.exit("Your PIGE file doesn't exist or is not formatted correctly (should include PIGEaucRepv3.xlsx)")

your_path_PIGE = your_folder_path + '\\' + your_PIGE_file[0]

# df = pd.read_excel(your_path_PIXE,sheet_name=int(sheet_number),index_col=0)
df = pd.read_csv(your_path_PIXE)
df_PIGE = pd.read_excel(your_path_PIGE, sheet_name=0, header=0)
# endregion

# region Formatting data sheets

# Removes last 4 digits of sample labels
df['Filename'] = df['Filename'].str[:-4]

df_PIGE = df_PIGE.rename(columns={'Sample': 'Filename'})

if df_PIGE['Filename'].iloc[1][-3:] == 'Spe':
    df_PIGE['Filename'] = df_PIGE['Filename'].str[:-4]

df_PIGE.drop(df_PIGE.columns[[1,2,3]], axis = 1, inplace = True)
# endregion

# region Combines PIGE and PIXE data
# This combines the PIGE and PIXE files by their filename
# from the PIGE file


df_combined = df.merge(
    df_PIGE)
combined_cols = df_combined.shape[1]

for i in range(7, combined_cols-1): #This loops over all elements besides Ar to normalize to 770 counts
    df_combined.iloc[:, i] = df_combined.iloc[:, i] * \
        df_combined['770 keV counts'].mean()/df_combined['770 keV counts']

df_combined_1 = df_combined.copy()
# endregion

# region Finds which samples are the triplicate samples

if triplicates == 'Yes':
    # Take Filename and remove last 2 letters
    # removes triplicate number
    df_combined['Filename'] = df_combined['Filename'].str[:-2]

    # Finds occurances of each sample name
    sample_names = Counter(df_combined['Filename']).most_common()
    replicate_names = []


    for x in range(len(sample_names)):  # Finds all triplicates
        if sample_names[x][1] == 3:
            replicate_names.append(sample_names[x][0])

file_names = df_combined['Filename']

# endregion

# region Finds the row indices for each triplicate
replicate_rows = []
blank_rows =[]
nist_2781_rows = []
am_soils_rows = []

if triplicates == "Yes":
    


    for x in replicate_names:
        replicate_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])
        if 'Blank' in x or 'blank' in x or 'BLK' in x:
            blank_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])
        if 'NIST'in x or '2781' in x: #Creates own list for NIST 2781 STD rows
            nist_2781_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])
        if 'SOIL' in x:
            am_soils_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])

if triplicates == "No":
    


    for x in file_names:
        
        if 'Blank' in x or 'blank' in x or 'BLK' in x:
            blank_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])
        if 'NIST'in x or '2781' in x: #Creates own list for NIST 2781 STD rows
            nist_2781_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])
        if 'SOIL' in x:
            am_soils_rows.append(df_combined.index[df_combined['Filename'].str.contains(x)])

# endregion

# region Find H value to scale all PIGE and PIXE elements by according to NIST 2781

nist_2781_official = [3100,750,28000,628,1272]
nist_2781_elements = ['TiK','MnK','FeKA','CuK','ZnK']
am_soil_reference = [1600,200,8400]
am_soil_elements = ['TiK','MnK','FeKA']
ratios = []


def standard_ratios(rows,elements,reference): #creates function to perform standardizations to be used for different standards if available
    df_standards = pd.DataFrame()
    for x in rows: #Same as statistics calculations below, this is purely for correcting to NIST standards
        std_means = []
        for i in range(7,df_combined.shape[1]):
            std_mean = df_combined.iloc[x,i].mean()
            std_means.append(std_mean)
        df_standards = pd.concat([df_standards,pd.DataFrame(std_means)],axis=1)
        
    df_standards = pd.DataFrame.transpose(df_standards)
    df_standards.columns = df_combined.columns[7:]

    measured = []
    for x in elements: #iterates over favored elements to pick out those averages
        chosen_elements = df_standards.filter(like = x, axis = 1).iloc[0][x]
        measured.append(chosen_elements)

    for x in range(len(reference)): #creates list of ratios
        ratios.append(measured[x]/reference[x])

    ratio_mean = statistics.fmean(ratios) #finds the mean of the ratios

    return 1.2/ratio_mean #factor to multiply to correct to standards

h_factor = 0

if nist_2781_rows != []:
    h_factor = standard_ratios(nist_2781_rows,nist_2781_elements,nist_2781_official)

if nist_2781_rows == [] and am_soils_rows != []:
   h_factor = standard_ratios(am_soils_rows,am_soil_elements,am_soil_reference)



# endregion

# region This Scales all elemental values besides the 770 peak with the H value

if h_factor != 0 and h_needed == 'Yes':

    for i in range(7, combined_cols-1):
        df_combined_1.iloc[:, i] = df_combined.iloc[:, i] * h_factor * con
        
df_combined_2 = df_combined_1.copy()

#endregion 


# region Calculates stats for elements and creates new Dataframe to hold it

#Sets up new pandas dataframes to be used below
df_stats = pd.DataFrame()
df_p = pd.DataFrame()
df_uncertainties = pd.DataFrame()
df_averages = pd.DataFrame()

if triplicates == "Yes":
    for x in replicate_rows:  # For loop to loop over replicate rows and calculate important stats for each triplicate
        element_means = []
        element_stds = []
        element_errors = []
        element_p_values = []

        for i in range(7, df_combined.shape[1]): #This loops over all elements
            element_mean = df_combined_2.iloc[x, i].mean()
            element_std = df_combined_2.iloc[x, i].std()

            if blank_rows[0].shape == x.shape: #In case no blanks are found, it doesn't calculate p values in comparison to blanks
                element_p_value = scipy.stats.ttest_ind(df_combined_2.iloc[x,i],df_combined_2.iloc[blank_rows[0],i])

            if blank_rows[0].shape != x.shape: #Related to above
                element_p_value = 2

            if element_mean != 0: #These two if statements are to deal with zeros in your uncertainties due to means equaling zero, which PMF does not like
                element_errors.append(int(element_std/element_mean*100))
            if element_mean == 0:
                element_errors.append(100)

            if element_p_value[1] < 100: #These two if statements are if the T test results in an nan
                element_p_values.append(round(element_p_value[1],3))
            if (element_p_value[1] < 100) == False:
                element_p_values.append('nan')

            if element_std <= 2: #In case element std is zero, sets it to arbitrary number, again to satisfy PMF
                element_stds.append(10)
            if element_std > 2:
                element_stds.append(int(element_std))

            element_means.append(int(element_mean))
            
        
            
        df_uncertainties = pd.concat([df_uncertainties, pd.DataFrame(element_stds)], axis=1) #This is to give the STD as uncertainty for each replicate
        df_averages = pd.concat([df_averages,pd.DataFrame(element_means)],axis=1) #This to give new dataframe for means. This and above are used for PMF
        df_stats = pd.concat([df_stats, pd.DataFrame(element_means), pd.DataFrame(element_stds), pd.DataFrame(element_errors), pd.DataFrame(element_p_values)], axis=1)
        # df_p = pd.concat([df_p,pd.DataFrame(element_p_values)])

    # Transpose because dataframe is flipped
    df_stats = pd.DataFrame.transpose(df_stats)
    df_uncertainties = pd.DataFrame.transpose(df_uncertainties)
    df_averages = pd.DataFrame.transpose(df_averages)
    # Relabels the columns to the elements
    df_stats.columns = df_combined.columns[7:]
    df_uncertainties.columns = df_combined.columns[7:]
    df_averages.columns = df_combined.columns[7:]


    stats_names = []
    uncertainty_names = []
    for x in replicate_names:  # Creates labels list for the stats using replicate names
        stats_names.append(x + ' Average')
        stats_names.append(x + ' Standard Deviation')
        stats_names.append(x + ' Percent Error')
        stats_names.append(x + ' p Value compare to Blank')
        uncertainty_names.append(x)

    #These three below insert the label lists created above into dataframe
    df_stats.insert(0,'Labels',stats_names)
    df_uncertainties.insert(0,'Labels', uncertainty_names)
    df_averages.insert(0,'Labels',uncertainty_names)

    df_stats = df_stats.set_index('Labels', drop=False).rename_axis(
        None)  # renames rows in dataframe

    df_uncertainties = df_uncertainties.set_index('Labels', drop=False).rename_axis(None)

    df_averages = df_averages.set_index('Labels', drop=False).rename_axis(None)

# endregion

# region Writes out the Dataframes to an excel file 

folder_name = os.path.basename(os.path.dirname(your_path_PIXE))
folder_path = os.path.dirname(os.path.abspath(your_path_PIXE))

df_combined_1 = df_combined_1.astype('int32', errors= 'ignore')


with pd.ExcelWriter(folder_path + "\\" + folder_name + "-CombinedPIGEPIXEv8.xlsx", engine= 'xlsxwriter') as writer:

    df_combined_1.to_excel(writer, sheet_name = "AR Normed Data", index = False)
    if triplicates == "Yes":
        df_stats.to_excel(writer, sheet_name = "Ar Normed Stats", index = False)
        df_averages.to_excel(writer, sheet_name = "Averages", index = False)
        df_uncertainties.to_excel(writer, sheet_name = "Uncertainties", index = False)


#endregion