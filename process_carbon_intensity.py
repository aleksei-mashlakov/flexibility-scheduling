import pandas as pd
df = pd.read_csv('./data/UK_df_fuel_ckan.csv')
df['DATETIME'] = pd.DatetimeIndex(df['DATETIME'])
df.set_index('DATETIME',inplace=True)
df_2019 = df['2019-01-01 00:00:00':'2019-12-31 23:30:00']
df_2019[['CARBON_INTENSITY']].to_csv('./data/Carbon_intensity_GB_2019.csv')
