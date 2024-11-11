import pandas as pd

class DataLoader:
    def __init__(self, parent_dir='labs/mwinslow/Karen/Bootstrapping_analysis', exp_name, outfile):
        self.raw_df_path = f'{parent_dir}/exp_name/Output/outfile'
        self.raw_df = pd.read_csv(self.raw_df_path)