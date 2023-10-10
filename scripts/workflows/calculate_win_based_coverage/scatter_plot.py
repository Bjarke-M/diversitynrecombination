import pandas as pd
import os


# Define the input path
input_file = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/PDGP_metadata.txt'
data_dir = '/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/mask/'

def generate_input_files(metadata, window_size):
    with open(metadata, 'r') as f:
        next(f)
        list_of_files = {}
        for line in f:
            PDGP_ID, Genus, Species, FROH, Sex, ref_assembly = line.strip().split(',')
            if ref_assembly == 'Daubentonia_madagascariensis':
                file_path = os.path.join(data_dir, ref_assembly, PDGP_ID, 'csv' , f'{PDGP_ID}_masked_{window_size}.csv')
                # Include sex, species, and ID in the dictionary
                list_of_files[PDGP_ID] = {
                    'Genus': Genus,
                    'Species': Species,
                    'FROH': FROH,
                    'Sex': Sex,
                    'ref_assembly': ref_assembly,
                    'file_path': file_path
                }
    return list_of_files

def extract_data_from_files(list_of_files, window_size):
    dfs = {}
    for ind, info in list_of_files.items():  # Iterate over the dictionary items
        #print(ind)
        file_path = info['file_path']
        #print(file_path)
    #     file_path = info['file_path']
       # Read the CSV file
        data = pd.read_csv(file_path, header=0)  # Assuming the first row contains column headers
        print(data[:5:])
        if ind not in dfs:
            dfs[ind] = {}
        median_freq_by_chr = data.groupby('chrA')['freq'].median().to_dict()
        #print(median_freq_by_chr)
        # Add the result to the list of DataFrames
        dfs[ind][window_size]=median_freq_by_chr
        print(dfs)
    return dfs



import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def create_scatter_plot(data_dict):
    # Convert the dictionary to a pandas DataFrame for easier manipulation
    df = pd.DataFrame(data_dict)

    # Reset index and rename columns for better visualization
    df = df.reset_index()
    df = df.rename(columns={"index": "Window Size"})

    # Melt the DataFrame to make it suitable for scatter plot and faceting
    df = pd.melt(df, id_vars=["Window Size"], var_name="Chromosome", value_name="Frequency")

    # Create a scatter plot using seaborn
    sns.set(style="whitegrid")
    g = sns.relplot(
        data=df,
        x="Window Size",
        y="Frequency",
        hue="Chromosome",
        kind="scatter",
        facet_kws=dict(sharex=False),
        height=5,
        aspect=1.5
    )

    # Customize plot labels and title
    g.set_axis_labels("Window Size", "Frequency")
    g.set_titles("Chromosome {col_name}")

    # Show the plot
    plt.show()

create_scatter_plot(extract_data_from_files(generate_input_files(input_file, 2000000),2000000))