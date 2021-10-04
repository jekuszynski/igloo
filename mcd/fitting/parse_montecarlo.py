import sys
import os
import pandas as pd
from matplotlib import gridspec, pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def parse_params(param_data,peak_count):
    all_peak_data = {}
    for n in range(1,peak_count+1):
        all_peak_data["peak_{0}_data".format(n)] = pd.DataFrame()
        print('Parsing Peak {0}...'.format(n))
        for i, param in param_data.iterrows():
            peak_num = i%3 + 1
            if peak_num == n:
                all_peak_data["peak_{0}_data".format(n)] = pd.concat([all_peak_data["peak_{0}_data".format(n)], param.to_frame().T])
    return all_peak_data
        
def plot_params(all_params,material):
    plt.clf()
    peak_count = (len(all_params))
    fig,ax=plt.subplots(4, peak_count, figsize=(8,8))
    for peak_num, peak in enumerate(all_params.keys()):
        y = all_params[peak]['redchi']
        for i in range(all_params[peak].shape[1]):
            for n in range(2,9,2):
                if n == i:
                    x = all_params[peak].iloc[:,n]
                    xerr = all_params[peak].iloc[:,n+1]
                    param_name = pd.Index(x).name
                    row = int(i/2 - 1)
                    ax[row,peak_num].errorbar(x, y, fmt='.')
                    ax[row,peak_num].xaxis.set_minor_locator(AutoMinorLocator())
                    ax[row,peak_num].yaxis.set_minor_locator(AutoMinorLocator())
                    ax[row,0].set_ylabel(param_name)
                    ax[0,peak_num].set_title('Peak {0}'.format(peak_num+1))
                    ax[-1,peak_num].set_xlabel(r'Reduced ${\chi}^{2}$')
    plt.tight_layout()
    plt.savefig(material + '_' + str(peak_count) + '_monte_carlo_params', dpi=300)

if __name__ == '__main__':

    material = '5-1CFS'
    peak_count = '3'

    working_path = '/home/jkusz/github/igloo/mcd/fitting/temp_param/'
    os.chdir(working_path)

    '''Read Data'''
    param_data_path = '/home/jkusz/github/igloo/mcd/fitting/temp_param/5-1CFS_parameter_testing.csv'
    param_data = pd.read_csv(param_data_path, index_col=0)

    '''Parse Data'''
    all_params = parse_params(param_data,int(peak_count))

    plot_params(all_params,material)

    '''Read and Organize Data'''

    print("All Done~!")
