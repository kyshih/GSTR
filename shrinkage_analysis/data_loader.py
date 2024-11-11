# shrinkage_analysis/data_loader.py
import pandas as pd
import argparse
from shrinkage_analysis.data_utils import (
    get_excluded_samples_from_file,
    filter_df,
    find_unique_values
)

class DataLoader:
    def __init__(self, raw_df_path, discard_samples_path, pattern='Safe|Neo|NT'):
        self.raw_df_path = raw_df_path
        self.discard_samples_path = discard_samples_path
        self.pattern = pattern

    def load_excluded_samples(self):
        return get_excluded_samples_from_file(self.discard_samples_path)
    
    def load_and_exclude_samples(self):
        excluded_samples = self.load_excluded_samples()
        df = pd.read_csv(self.raw_df_path)
        return filter_df(df, 'Sample_ID', excluded_samples, exclude=True)
    
    def get_gRNA_data(self, df):
        return filter_df(df, 'Identity', ['gRNA'], exclude=False)
    
    def find_control_gRNAs(self, df):
        from bootstrapping_helpers import Find_Controls
        return Find_Controls(df, self.pattern)
        
    @staticmethod
    def parse_args():
        parser = argparse.ArgumentParser(description='Load and clean data for analysis')
        parser.add_argument('--raw_df_path', required=True, help='Path to the raw data CSV file')
        parser.add_argument('--discard_samples_path', required=True, help='Path to the file with samples to exclude')
        parser.add_argument('--pattern', default='Safe|Neo|NT', help='Pattern to match control gRNAs')
        return parser.parse_args()

    @classmethod
    def from_args(cls):
        args = cls.parse_args()
        return cls(args.raw_df_path, args.discard_samples_path, args.pattern)

# Example usage when running as a script
if __name__ == "__main__":
    data_loader = DataLoader.from_args()
    cleaned_data = data_loader.load_and_clean_data()
    print("Loaded and cleaned data:")
    print(cleaned_data.head())

    # Example usage of control gRNA fetching
    control_gRNAs = data_loader.find_control_gRNAs(cleaned_data)
    print("Control gRNAs found:", control_gRNAs)
    total_gRNAs = data_loader.get_total_gRNA(cleaned_data)
    print("Total gRNA is:", total_gRNAs)