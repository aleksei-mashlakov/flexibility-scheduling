import glob
import pandas as pd
import numpy as np
import os, sys
import pickle

class BESS(object):
    """ Limited droop-curve model of BESS
    """
    def __init__(self, energy, power, init_soc, efficiency):
        self.soc = init_soc
        self.f = 50.0
        self.B = 49.75
        self.A_p = 0.09*power
        self.B_p = 0.4444444*power
        self.E = 50.25
        self.f_low_limit = 49.5
        self.dead_band = 0.05
        self.f_zero = 50.0
        self.f_high_limit = 50.5
        self.max_e_capacity = energy
        self.efficiency = efficiency
        self.energy_capacity = (self.max_e_capacity * self.soc)/100
        self.power_capacity = power
        self.slope_coeff = self.power_capacity/(self.f_high_limit-(self.f_zero+self.dead_band))
        self.slope_coeff_UE = (self.B_p-self.A_p)/(self.E-(self.f_zero+self.dead_band))
        self.slope_coeff_LE = (self.B_p+self.A_p)/(self.E-(self.f_zero+self.dead_band))
        self.timer_d = 0
        self.timer_u = 0
        self.p = 0

    def follow_frequency(self, f):
        ## IDLE
        if (f >= (self.f_zero - self.dead_band)) and (f <= (self.f_zero + self.dead_band)):
            #print('Zone A')
            self.timer_d = 0
            self.timer_u = 0
            power_LE = - self.A_p
            power_UE = self.A_p
            if self.p == 0:
                power = 0.0
            elif self.p > 0:
                power = max([0, self.p-0.01*self.power_capacity])
                power = min([power, self.A_p])
            elif self.p < 0:
                power = min([0, self.p+0.01*self.power_capacity])
                power = max([power, -self.A_p])
            else:
                ValueError('Frequency value is unidentified in zone A.')

        ## DOWN regulation
        elif (f > (self.f_zero + self.dead_band)):
            #print('Zone B-D')
            if f >= self.f_high_limit:
                power = -self.power_capacity
                self.timer_d += 1
                power_LE = power
                power_UE = power

            elif (f < self.f_high_limit) and (f > self.E):
                power = -(f-(self.f_zero+self.dead_band)) * self.slope_coeff
                self.timer_d = 0
                power_LE = power
                power_UE = power

            elif (f <= self.E):
                #print('Zone B: DOWN-regulation')
                self.timer_d = 0
                power_LE = - self.A_p -(f-(self.f_zero+self.dead_band)) * self.slope_coeff_UE
                power_UE = self.A_p -(f-(self.f_zero+self.dead_band)) * self.slope_coeff_LE
                power = self.p+(-(1/0.45)*(f-self.f))*self.power_capacity
                if not ((power <= power_UE) and (power >= power_LE)):
                    #print('Within boundaries')
                    if (power >= power_UE):
                        power = power_UE
                        #print('Above UE boundaries')
                    elif (power <= power_LE):
                        power = power_LE
                        #print('Below LE boundaries')
            else:
                ValueError('Frequency value is unidentified in zone B-Down.')

        ## UP regulation
        elif (f < (self.f_zero - self.dead_band)):
            #print('Zone B-U')
            if f <= self.f_low_limit:
                power = self.power_capacity
                self.timer_u += 1
                power_LE = power
                power_UE = power

            elif (f > self.f_low_limit) and (f < self.B):
                power = self.power_capacity - ((f-self.f_low_limit) * self.slope_coeff)
                self.timer_u = 0
                power_LE = power
                power_UE = power

            elif (f >= self.B):
                #print('Zone B: UP-regulation')
                self.timer_u = 0
                power_UE = self.B_p  - ((f-self.B) * self.slope_coeff_UE)
                power_LE = self.B_p  - ((f-self.B) * self.slope_coeff_LE)
                power = self.p + (-(1/0.45)*(f-self.f))*self.power_capacity

                if not ((power <= power_UE) and (power >= power_LE)):
                        #print('Within boundaries')
                    if (power >= power_UE):
                        power = power_UE
                        #print('Above UE boundaries')
                    elif (power <= power_LE):
                        power = power_LE
                        #print('Below LE boundaries')
            else:
                ValueError('Frequency value is unidentified in zone B-Up.')
        power_ch = 0
        power_dc = 0
        if power < 0:
            power_ch = power
        elif power >= 0:
            power_dc = power
        delta_energy = -(power_ch*self.efficiency + power_dc/self.efficiency) * (1.0/3600.0)
        self.soc = delta_energy

        if (self.timer_u > 60*15) or (self.timer_d > 60*15):
            print('MAX duration reached')
        self.f = f
        self.p = power
        return round(self.soc, 9), round(power, 6), round(power_LE, 6), round(power_UE, 6)

    def activate(self, filename, year):
        if year==2014:
            df = pd.read_csv(filename, sep=",", parse_dates=['dtm'],
                             infer_datetime_format=True,
                             na_values=['nan', '?'], index_col='dtm', dayfirst=True)
        else:
            df = pd.read_csv(filename, sep=",", parse_dates=['dtm'],
                             infer_datetime_format=True,
                             na_values=['nan', '?'], index_col='dtm', dayfirst=True,
                             date_parser=lambda x: pd.to_datetime(x.rsplit(' ', 1)[0]))

        df.rename(columns={​​​​​​​​'f': 'frequency'}​​​​​​​​, inplace=True)
        df['frequency'].replace('N/A', np.nan, inplace=True)
        df['frequency'].replace('INVA', np.nan, inplace=True)
        df['frequency'].fillna(method='ffill', inplace=True)
        soc_measurements=[]
        power_measurements=[]
        power_UE_list = []
        power_LE_list = []
        for index, row in df.iterrows():
            f = row['frequency']
            soc, power, power_LE, power_UE = self.follow_frequency(f)
            soc_measurements.append(soc)
            power_measurements.append(power)
            power_LE_list.append(power_LE)
            power_UE_list.append(power_UE)
        df['soc'] = soc_measurements
        df['power'] = power_measurements
        df['power_LE'] = power_LE_list
        df['power_UE'] = power_UE_list
        print(df.head())
        return df

years_list = [2019]
repository_address = r"./data"
project_address = r"./data"
delta_time = '30T'#'15T'

def get_date_indexes(start_year):
    date_list = []
    for i in range(1,13):
        date_list.append(str(start_year)+' '+str(i))
    return date_list

# prepare data
dfList_1s_y=[]
dfList_30m_y=[]
bess = BESS(7.5*0.9, 3.3, 0.0, 0.93)
# bess = BESS(13.5, 3.68, 0.0, 0.95)
for year in years_list:
    dfList_1s = []
    dfList_30m = []
    print(year)
    for month in get_date_indexes(year):
        print(month)
        filename = glob.glob(os.path.join(repository_address, str(year), month, 'f '+month+'.csv'))[0]
        df = bess.activate(filename, year)
        df['soc_up'] = df['soc']
        df['soc_up'][df['soc_up']>=0] = 0
        df['soc_down'] = df['soc']
        df['soc_down'][df['soc_down']<=0] = 0
        df2 = df['soc_up'].resample(delta_time).sum()
        df3 = df['soc_down'].resample(delta_time).sum()
        df5 = df[['power']].rename(columns={​​​​​​​​'power':'max_power'}​​​​​​​​).resample(delta_time).max()
        df6 = df[['power']].rename(columns={​​​​​​​​'power':'min_power'}​​​​​​​​).resample(delta_time).min()
        df3 = pd.concat([df2, df3, df5, df6], axis=1).fillna(method='ffill')
        dfList_1s.append(df.drop(['frequency','soc','soc_up','soc_down'], axis=1))
        dfList_30m.append(df3)

        with open('./1s_sim.pkl', 'wb') as b:
            pickle.dump(dfList_1s,b)
        with open('./30m_sim.pkl', 'wb') as b:
            pickle.dump(dfList_30m,b)

    concatDf_1s = pd.concat(dfList_1s, axis=0)
    concatDf_30m = pd.concat(dfList_30m, axis=0)
    dfList_1s_y.append(concatDf_1s)
    dfList_30m_y.append(concatDf_30m)


concatDf_1s = pd.concat(dfList_1s_y, axis=0)
concatDf_30m = pd.concat(dfList_30m_y, axis=0)
savefile_1s = os.path.join(repository_address, 'GB-EFR-1-s-2019-tesla135.csv')
savefile_30m = os.path.join(repository_address, 'GB-EFR-30-m-2019-tesla135.csv')
concatDf_1s.to_csv(savefile_1s, sep=',', encoding='utf-8', header='column_names')
concatDf_30m.to_csv(savefile_30m, sep=',', encoding='utf-8', header='column_names')
