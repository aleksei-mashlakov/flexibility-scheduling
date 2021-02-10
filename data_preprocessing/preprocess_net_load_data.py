#!/usr/bin/env python
# coding: utf-8

# ## Import

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import os
pd.set_option('display.max_columns', None)

# ## Data TC5

data_rep1 = './datasets/TC5'

# ## Data TC20IHD

data_rep2 = './datasets/TC20IHD'


def main():
    # ### Net load

    df_ = pd.read_csv('./datasets/TC5/{}.csv'.format('whole home power import'),
                      index_col=['timestamp'], parse_dates=True)
    # reverse measurements
    df_['25222'] = df_['25222']*(-1)


    df_2 = pd.read_csv('./datasets/TC20IHD/{}.csv'.format('whole home power import'),
                      index_col=['timestamp'], parse_dates=True)
    # reverse measurements
    df_2['95032'] = df_2['95032']*(-1)
    # from visual inspection the measurements can be corrupted by the datetime index and pattern can be restored ignoring the missing gaps but no clear way to prove the correct reorganisation is available

    # #### Concat TC5 and TC20IHD

    df_ = pd.concat([df_,df_2], axis=1, sort=True)
    # measurements are partilly reversed
    df_['95112']['2012-12-13 00:00:00':'2012-12-24 12:58:00'] = df_['95112']['2012-12-13 00:00:00':'2012-12-24 12:58:00']*(-1)
    data_rep = './datasets'
    df_.reset_index().to_csv(os.path.join(data_rep,'{} 1 min full dataset.csv'.format('whole home power import')), index=False)

    # #### Total number of customers
    print(f"Total number of customers: {len(df_.columns)}")


    # #### Percentage of missing data
    df = df_['2012-12-13 00:00:00':'2014-01-08 23:59:00']

    train = df
    total = train.isnull().sum().sort_values(ascending = False)
    percent = (train.isnull().sum()/train.isnull().count()*100).sort_values(ascending = False)
    missing__train_data  = pd.concat([total, percent], axis=1, keys=['Total', 'Percent'])
    # missing__train_data.head()


    # missing__train_data['Percent'].plot()

    print("Number of customers with less than 6 % missing data:")
    print(missing__train_data[missing__train_data['Percent']<6].shape[0])


    # #### Replace NaNs with mean minute values per column

    import numpy as np
    df = df[missing__train_data[missing__train_data['Percent']<6].index]
    for i, col in enumerate(df):
        a = np.reshape(df[col].values, (-1, 24*60))
        col_mean = np.nanmean(a, axis=0)
        inds = np.where(np.isnan(a))
        a[inds] = np.take(col_mean, inds[1])
        df.loc[:, col] = a.reshape(-1)

    data_rep = './datasets'
    meas = 'whole home power import 1 min 12 months averaged'
    df.reset_index().to_csv(os.path.join(data_rep,'{}.csv'.format(meas)), index=False)
    meas = 'whole home power import 30 min 12 months averaged'
    df.resample('30T').mean().reset_index().dropna().to_csv(os.path.join(data_rep,'{}.csv'.format(meas)), index=False)

if __name__ == '__main__':
    main()
