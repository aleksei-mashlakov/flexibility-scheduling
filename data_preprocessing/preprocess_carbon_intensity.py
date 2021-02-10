import pandas as pd

def main():
    df = pd.read_csv('./datasets/UK_df_fuel_ckan.csv')
    df['DATETIME'] = pd.DatetimeIndex(df['DATETIME'])
    df.set_index('DATETIME',inplace=True)
    df_2019 = df['2019-01-01 00:00:00':'2019-12-31 23:30:00']
    df_2019[['CARBON_INTENSITY']].to_csv('./datasets/Carbon_intensity_GB_2019.csv')

if __name__ == '__main__':
    main()
