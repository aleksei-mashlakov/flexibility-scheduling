print('we are here')
# ### Online battery validation

# In[342]:


import os
import glob
import pandas as pd
#import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# %matplotlib inline

import pickle


# In[425]:


class BESS(object):
    def __init__(self, max_energy, max_power, init_soc_proc, efficiency):
        self.soc = init_soc_proc
        self.max_e_capacity = max_energy
        self.efficiency = efficiency
        self.energy = self.max_e_capacity * (self.soc)/100
        self.power = max_power

    def online_control(self):
        """
        """

        return

    def calculate_NLF(self, net_load_day):
        """ Net load factor
        """
        df = pd.DataFrame(net_load_day).abs()
        NLF = df.mean()/df.max()
        return NLF[0]

    def calculate_SBSPM(self, NR, LE, UE, error=0.01):
        """
        Calculates second by second Service Performance Measure (SBSPM)
        """
        if (NR >= LE - error) and (NR <= UE + error):
            SBSPM = 1
        elif (NR > UE + error):
            SBSPM = max([1-abs(NR - UE), 0])
        elif (NR < LE - error):
            SBSPM = max([1-abs(NR - LE), 0])
        else:
            raise ValueError('The NR is undefined {}'.format(NR))
        return SBSPM

    def average_SPM_over_SP(self, SBSPM_list):
        """
        Averages SPM over Settlement period
        """
        SPM = sum(SBSPM_list)/1800
        return SPM

    def check_availability(self, SPM):
        """
        Returns availability factor
        """
        if SPM >= 0.95:
            availability_factor = 1
        elif (SPM >= 0.75) and (SPM < 0.95):
            availability_factor = 0.75
        elif (SPM >= 0.5) and (SPM < 0.75):
            availability_factor = 0.5
        elif (SPM < 0.5):
            availability_factor = 0
        return availability_factor


def save_to_pickle(name, list_to_save, save_path):
    with open(os.path.join(save_path, '{}.pkl'.format(name)), 'wb') as f:
        pickle.dump(list_to_save, f)
    return

def load_from_pickle(name, save_path):
    with open(os.path.join(save_path, '{}.pkl'.format(name)), 'rb') as f:
        p = pickle.load(f)
    return p



import os
import pandas as pd
path = "."
bess_name = "sonnen"
apath = os.path.join(path, 'simulations_{}'.format(bess_name), '{}'.format(1),'agent_{}.csv'.format(1))
nl = pd.read_csv(apath).loc[:,['nl_a{}'.format(1)]]
pb = pd.read_csv(apath).loc[:,['pb_a{}'.format(1)]]
c_reg = pd.read_csv(apath).loc[:,['c_reg_a{}'.format(1)]]
#flex_down = pd.read_csv(apath).loc[:,['flex_down_a{}'.format(1)]]
#flex_up = pd.read_csv(apath).loc[:,['flex_up_a{}'.format(1)]]


apath = os.path.join(path, 'forecasts', bess_name, 'min_power_forecast.csv')
min_efr_forecast = pd.read_csv(apath, sep=",")
min_efr_forecast[min_efr_forecast<0] = 0

apath = os.path.join(path, 'forecasts', bess_name, 'max_power_forecast.csv')
max_efr_forecast = pd.read_csv(apath, sep=",")
max_efr_forecast[max_efr_forecast<0] = 0


net_load_path = '.'
apath = os.path.join(net_load_path,'whole home power import 1 min 12 months averaged.csv')
net_load_true = pd.read_csv(apath, index_col=['timestamp'], parse_dates=True)
net_load_true = net_load_true[net_load_true.iloc[0:100*60*24+1,:].index[-1]:]


# In[51]:


# project_address = r"C:\Users\h17353\PycharmProjects\Frequency\frequency_measurements\GB"
if bess_name == 'sonnen':
    apath = os.path.join(net_load_path, 'GB-EFR-1-s-2019-sonnenEco75.csv')
elif bess_name == 'tesla':
    apath = os.path.join(net_load_path, 'GB-EFR-1-s-2019-tesla135.csv')

efr_true = pd.read_csv(apath, sep=",", parse_dates=['dtm'],infer_datetime_format=True,
                             na_values=['nan', '?'], index_col='dtm', dayfirst=True)
efr_true = efr_true['2019-04-11 00:00:00':]


# In[511]:

priority = 'schedule' #'efr' #
dataset_folder = '.'
sim_path = "."
save_path = "./save_results_09998"
len_betas = 1
len_agents = 150
len_days = 150
HH=48
M=30
T=60


## Initialize list of dataframes
dfs_SPM = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_avail_fact = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_bat_energy = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_bat_imb = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_nl_imb = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_frcst_imb = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))

dfs_SPM_std = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_avail_fact_std = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_bat_energy_std = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_bat_imb_std = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_nl_imb_std = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))
dfs_frcst_imb_std = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))

df_NLF_agg_real = pd.DataFrame(columns=range(1,len_days+1,1))
df_NLF_agg_frcst = pd.DataFrame(columns=range(1,len_days+1,1))
dfs_real_NLF = pd.DataFrame(index=range(1,len_days+1,1), columns=range(len_agents))
dfs_frcst_NLF = pd.DataFrame(index=range(1,len_days+1,1), columns=range(len_agents))
df_agg_imb = pd.DataFrame(index=range(HH), columns=range(1,len_days+1,1))


dfs_daily_nl_imb = pd.DataFrame(index=range(1,len_days+1,1), columns=range(len_agents))
dfs_daily_bat_imb = pd.DataFrame(index=range(1,len_days+1,1), columns=range(len_agents))
dfs_daily_frcst_imb = pd.DataFrame(index=range(1,len_days+1,1), columns=range(len_agents))

for day in range(1,len_days+1,1):
    print('day - {}'.format(day))
    ## Initialize daily profiles for all agents (EFR)
    efr_true_day = efr_true.iloc[(24*60*60)*(day-1):(24*60*60)*day,]

    ## Initialize daily profiles per agent
    net_load_true_day = net_load_true.iloc[(24*60)*(day-1):(24*60)*day,]

    ## Create dataframes to save Metric(agents)
    df_SPM = pd.DataFrame(index=range(HH), columns=range(len_agents))
    df_avail_fact = pd.DataFrame(index=range(HH), columns=range(len_agents))
    df_frcst_imbs = pd.DataFrame(index=range(HH), columns=range(len_agents))
    df_bat_energy = pd.DataFrame(index=range(HH), columns=range(len_agents))
    df_nl_imb = pd.DataFrame(index=range(HH), columns=range(len_agents))
    df_bat_imb = pd.DataFrame(index=range(HH), columns=range(len_agents))

    for beta in range(len_betas):
        print('beta - {}'.format(beta))
        ## Getting average int index of selected plans for all agents
        beta_i = 0.9998
        path = os.path.join(dataset_folder, 'run_beta_day{}'.format(day),'output-{}-beta-{}'.format(day, beta_i))
        file_loc = glob.glob(os.path.join(path, next(os.walk(path))[1][0],'*selected-plans.csv'), recursive=True)
        file_cost = glob.glob(os.path.join(path, next(os.walk(path))[1][0],'*global-cost.csv'), recursive=True)
        #file_termin = glob.glob(os.path.join(path, next(os.walk(path))[1][0],'*termination.csv'), recursive=True)
        #termin = pd.read_csv(file_termin[0]).iloc[:,1].max()
        run = pd.read_csv(file_cost[0]).groupby(['Iteration']).min().iloc[-1,2:].idxmin().split('-')[1]
        df = pd.read_csv(file_loc[0])
        list_selected = df.loc[(df['Run']==int(run))&(df['Iteration']==df['Iteration'].max())].iloc[0,2:].values
        #df = pd.read_csv(file_loc[0]).groupby(['Iteration']).mean().drop('Run', axis=1)
        #list_selected = df.iloc[(df.index[-1]),:].astype(int).values #termin

        agg_nl_real = pd.Series(data=np.zeros(48))
        agg_nl_frcst = pd.Series(data=np.zeros(48))

        for agent in range(len_agents):
            #print('agent - {}'.format(agent))
            ## Initialize the battery instance
            if bess_name == 'sonnen':
                bess = BESS(0.9*7.5, 3.3, 50, 0.93)
            elif bess_name == 'tesla':
                bess = BESS(13.5, 3.68, 50, 0.95)
            else:
                raise ValueError('Bess name is not given')

            ## Extract plan id related net load, battery schedule, and EFR capacity for the agent
            plan_id = list_selected[agent]
            #print('Plan id {}'.format(plan_id))

            apath = os.path.join(sim_path, 'simulations_{}'.format(bess_name), '{}'.format(day),'agent_{}.csv'.format(agent))
            df_sim = pd.read_csv(apath)
            nl = df_sim.loc[:,['nl_a{}'.format(plan_id+1)]]
            pb = df_sim.loc[:,['pb_a{}'.format(plan_id+1)]]
            c_reg = df_sim.loc[:,['c_reg_a{}'.format(plan_id+1)]]

            min_efr_forecast_day = min_efr_forecast.iloc[(day-1)*(48):(day)*(48),[plan_id]]
            max_efr_forecast_day = max_efr_forecast.iloc[(day-1)*(48):(day)*(48),[plan_id]]

            ## Initialize the metrics
            net_load_day = []
            ava_factors_day = []
            frcst_imbs_day = []
            SPMs = []
            bat_imbs = []
            bat_energy_day = []

            ## Iterate through half hourly intervals
            for hh in range(HH):
                #print('-------- {} ----------'.format(hh))
                net_load_hh = []
                spm_hh = []
                frcst_imbs_hh = []
                bat_imbs_hh = []
                bat_energy_hh = []

                ## Obtain scheduled battery power and minmax response powers
                Pbat_sch = pb.iloc[hh, 0]
                P_down = min_efr_forecast_day.iloc[hh, 0]
                P_up = max_efr_forecast_day.iloc[hh, 0]
                C_reg = c_reg.iloc[hh, 0]

                for m in range(M):
                    ## Get minute-based net power
                    P_real = net_load_true_day.iloc[(hh*30)+m, agent]

                    for t in range(T):
                        ## Check the schedule limits violation
                        P_bat_resp_req = efr_true_day.iloc[((hh*30)+m)*60+t,0]*C_reg  ## 0 - 'power'
                        LE = efr_true_day.iloc[((hh*30)+m)*60+t,1]*C_reg              ## 1 - 'LE'
                        UE = efr_true_day.iloc[((hh*30)+m)*60+t,2]*C_reg              ## 2 - 'UE'

                        #print('LE {} - P_bat_resp_req {} - UE {}: (C_reg) {}'.format(LE, P_bat_resp_req, UE, C_reg))

                        ## down-regulation
                        if P_bat_resp_req < 0:
                            P_bat_resp = max([P_bat_resp_req, -P_down*C_reg])
                            efr_frcst_imb = P_bat_resp_req - P_bat_resp               ## neg imbalance

                        ## up-regulation
                        elif P_bat_resp_req >= 0:
                            P_bat_resp = min([P_bat_resp_req, P_up*C_reg])
                            efr_frcst_imb = P_bat_resp_req - P_bat_resp               ## pos imbalance
                        else:
                            raise ValueError('Error in P_bat_resp')

                        ## normalize
                        if C_reg!= 0:
                            efr_frcst_imb = efr_frcst_imb/(bess.power*C_reg)

                        ## Forecast imbalance
                        #print('P_down; P_bat_resp_req; P_up; efr_frcst_imb; (C_reg)')
                        #print('{}; {}; {}; {}; {}'.format(-P_down, P_bat_resp_req, P_up, efr_frcst_imb, C_reg))
                        frcst_imbs_hh.append(efr_frcst_imb)

                        ## Check the energy limits violation
                        Pbat = Pbat_sch + P_bat_resp_req
                        Pbat_ch = 0.
                        Pbat_dc = 0.
                        ## down-regulation
                        if Pbat < 0.:
                            Pbat_ch = max([Pbat, -bess.power, (bess.energy-bess.max_e_capacity)/(bess.efficiency*(1./3600.))])
                            power_imb = Pbat - Pbat_ch                                 ## neg or zero
                            if power_imb < 0.:
                                if P_bat_resp_req < 0.:
                                    if priority=='schedule':
                                        ## priority of battery schedule
                                        if (Pbat_sch-Pbat_ch) <= 0.:
                                            NR = 0.
                                            ## Pbat_sch = Pbat_ch
                                        else:
                                            NR = max([P_bat_resp_req, Pbat_ch - Pbat_sch])
                                            ## Pbat_sch = Pbat_sch
                                    elif priority=='efr':
                                        ## priority of efr
                                        if (P_bat_resp_req-Pbat_ch) <= 0.:
                                            NR = Pbat_ch
                                            ## Pbat_sch = 0
                                        else:
                                            NR = P_bat_resp_req
                                            ## Pbat_sch = Pbat_ch - P_bat_resp_req
                                    else:
                                        raise ValueError('Priority is not given')
                                elif P_bat_resp_req >= 0.:
                                    NR = P_bat_resp_req
                                    ## Pbat_sch = Pbat_sch - power_imb
                            elif power_imb == 0.:
                                NR = P_bat_resp_req
                            else:
                                raise ValueError('Error in power_imb')

                        elif Pbat >= 0.:
                            Pbat_dc = min([Pbat, bess.power, bess.efficiency*(bess.energy-0)/(1./3600.)])
                            power_imb = Pbat - Pbat_dc                                ## pos or zero
                            if power_imb > 0.:
                                if P_bat_resp_req > 0.:
                                    if priority=='schedule':
                                        ## priority of battery schedule
                                        if (Pbat_sch-Pbat_dc)>=0.:
                                            NR = 0.
                                            ## Pbat_sch = Pbat_ch
                                        else:
                                            NR = min([P_bat_resp_req, Pbat_dc - Pbat_sch])
                                            ## Pbat_sch = Pbat_sch
                                    elif priority=='efr':
                                        ## priority of efr
                                        if (P_bat_resp_req-Pbat_dc) >= 0.:
                                            NR = Pbat_dc
                                            ## Pbat_sch = 0
                                        else:
                                            NR = P_bat_resp_req
                                            ## Pbat_sch = Pbat_dc - P_bat_resp_req
                                    else:
                                        raise ValueError('Priority is not given')

                                elif P_bat_resp_req <= 0:
                                    NR = P_bat_resp_req
                                    ## Pbat_sch = Pbat_sch - power_imb
                            elif power_imb == 0:
                                NR = P_bat_resp_req
                            else:
                                raise ValueError('Error in power_imb')
                        else:
                            raise ValueError('Error in Pbat')

                        #print('Pbat_sch {} P_bat_resp_req {} NR {} power_imb {} bess_energy {}'.format(Pbat_sch,P_bat_resp_req,
                        #                                                                NR, power_imb, bess.energy))
                        bat_imbs_hh.append(power_imb)

                        ## Change the energy state
                        bess.energy += -(Pbat_ch*bess.efficiency + Pbat_dc/bess.efficiency) * (1.0/3600.0)

                        Pb_real = Pbat_ch + Pbat_dc
                        P_nl = P_real - Pb_real
                        #print('P_real - Pb_real = P_nl')
                        #print('{} - {} = {}'.format(P_real,Pb_real,P_nl))
                        net_load_hh.append(P_nl)
                        if C_reg!= 0:
                            NR = NR/(bess.power*C_reg)
                            LE = LE/(bess.power*C_reg)
                            UE = UE/(bess.power*C_reg)

                        SBSPM = bess.calculate_SBSPM(NR, LE, UE, error=0.01)
                        spm_hh.append(SBSPM)
                        bat_energy_hh.append(bess.energy)

                bat_energy_day.append(pd.Series(bat_energy_hh).mean())
                net_load_day.append(pd.Series(net_load_hh).mean())
                bat_imbs.append(pd.Series(bat_imbs_hh).mean())
                SPM = bess.average_SPM_over_SP(spm_hh)
                SPMs.append(SPM)
                ava_factors_day.append(bess.check_availability(SPM))
                frcst_imbs_day.append(pd.Series(frcst_imbs_hh).mean())

            df_bat_energy.loc[:, agent] = pd.Series(bat_energy_day)
            df_bat_imb.loc[:, agent] = pd.Series(bat_imbs)
            df_frcst_imbs.loc[:, agent] = pd.Series(frcst_imbs_day)
            df_SPM.loc[:, agent] = pd.Series(SPMs)
            df_avail_fact.loc[:, agent] = pd.Series(ava_factors_day)

            dfs_real_NLF.loc[day, agent] = bess.calculate_NLF(net_load_day)
            dfs_frcst_NLF.loc[day, agent] = bess.calculate_NLF(nl.iloc[:,0].values)

            nl_real = pd.Series(net_load_day)
            nl_frcst = nl.iloc[:,0]

            nl_imbalance = nl_frcst - nl_real
            df_nl_imb.loc[:, agent] = nl_imbalance

            dfs_daily_nl_imb.loc[day, agent] = nl_imbalance.abs().sum()
            dfs_daily_bat_imb.loc[day, agent] = pd.Series(bat_imbs).abs().sum()
            dfs_daily_frcst_imb.loc[day, agent] = pd.Series(frcst_imbs_day).abs().sum()

            agg_nl_real += nl_real
            agg_nl_frcst += nl_frcst


        df_NLF_agg_real.loc[0, day] = bess.calculate_NLF(agg_nl_real)
        df_NLF_agg_frcst.loc[0, day] = bess.calculate_NLF(agg_nl_frcst)

        agg_nl_imb = agg_nl_frcst - agg_nl_real
        df_agg_imb.loc[:, day] = agg_nl_imb

        dfs_SPM.loc[:,day] = df_SPM.mean(axis=1)
        dfs_avail_fact.loc[:,day] = df_avail_fact.mean(axis=1)
        dfs_frcst_imb.loc[:,day] = df_frcst_imbs.mean(axis=1)
        dfs_bat_energy.loc[:,day] = df_bat_energy.mean(axis=1)
        dfs_bat_imb.loc[:,day] = df_bat_imb.mean(axis=1)
        dfs_nl_imb.loc[:,day] = df_nl_imb.mean(axis=1)

        dfs_SPM_std.loc[:,day] = df_SPM.std(axis=1)
        dfs_avail_fact_std.loc[:,day] = df_avail_fact.std(axis=1)
        dfs_frcst_imb_std.loc[:,day] = df_frcst_imbs.std(axis=1)
        dfs_bat_energy_std.loc[:,day] = df_bat_energy.std(axis=1)
        dfs_bat_imb_std.loc[:,day] = df_bat_imb.std(axis=1)
        dfs_nl_imb_std.loc[:,day] = df_nl_imb.std(axis=1)


df_NLF_agg_real.to_pickle(os.path.join(save_path, 'NLF_agg_real.pkl'))
df_NLF_agg_frcst.to_pickle(os.path.join(save_path, 'NLF_agg_frcst.pkl'))
df_agg_imb.to_pickle(os.path.join(save_path, 'df_agg_imb.pkl'))
save_to_pickle('real_NLF', dfs_real_NLF, save_path)
save_to_pickle('frcst_NLF', dfs_frcst_NLF, save_path)

save_to_pickle('SPM', dfs_SPM, save_path)
save_to_pickle('AF', dfs_avail_fact, save_path)
save_to_pickle('frcst_imb', dfs_frcst_imb, save_path)
save_to_pickle('bat_energy', dfs_bat_energy, save_path)
save_to_pickle('bat_power_imb', dfs_bat_imb, save_path)
save_to_pickle('nl_imb', dfs_nl_imb, save_path)

save_to_pickle('SPM_std', dfs_SPM_std, save_path)
save_to_pickle('AF_std', dfs_avail_fact_std, save_path)
save_to_pickle('frcst_imb_std', dfs_frcst_imb_std, save_path)
save_to_pickle('bat_energy_std', dfs_bat_energy_std, save_path)
save_to_pickle('bat_power_imb_std', dfs_bat_imb_std, save_path)
save_to_pickle('nl_imb_std', dfs_nl_imb_std, save_path)

save_to_pickle('dfs_daily_nl_imb', dfs_daily_nl_imb, save_path)
save_to_pickle('dfs_daily_bat_imb', dfs_daily_bat_imb, save_path)
save_to_pickle('dfs_daily_frcst_imb', dfs_daily_frcst_imb, save_path)

