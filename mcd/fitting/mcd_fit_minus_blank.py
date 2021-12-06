from mcd_full_fitting_routine_primary import *

def findFieldAbsorption(signal, blank):
    abs_df = {}
    for name, signal_values in signal.items():
        abs_df = pd.DataFrame()
        abs_df[name + '_absorbance'] = 2 - np.log10(100 * signal_values['chpx'] / blank[name]['chpx'])
    plot_abs(abs_df,'raw',)

def plotAllFieldAbsorption(absorption_data):
    fig,ax=plt.subplots(figsize=(4,2))
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(0,10,6),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='Greys',norm=norm) 
    fig.colorbar(sm,ticks=range(0,11,2),label='H (T)') #make color bar based on H (T) for plot
    for df in absorption_data.values():
        #Dr. Seaborn or: How I Learned to Stop Worrying and Love sns.lineplot. Such efficiency. Much wow.
        sns.lineplot(data=df,x='energy',y='mdeg', linewidth=0.6,
                    hue='field',hue_norm=(0,10),
                    palette=sns.color_palette('Greys',as_cmap=True),
                    legend=None)
    if x_axis=='Energy (eV)':
        plt.xlabel(x_axis)
    if op=='raw':
        plt.title("Raw MCD")
    if op=='avg':
        plt.title("Difference MCD")
    plt.ylabel('MCD (mdeg)')
    plt.xlim(3.2,.55)
    baseline = lines.Line2D(range(6),np.zeros(1),c='black',ls='--',lw=0.6) #draw baseline at 0T
    ax.add_line(baseline) #add baseline to plot
    plt.style.use('seaborn-paper')
    plt.savefig('diff_mcd',dpi=100,transparent=True,bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    '''parse data'''
    material = 'Au-PEG'
    signal_data = parse_mcd('C:\Users\roflc\Dropbox\Research\FSU\Strouse\Projects\MCD\Fringe Field Testing\11112021_Au-PEG_400_800nm_PMT_cuvette')
    blank_data = parse_mcd('C:\Users\roflc\Dropbox\Research\FSU\Strouse\Projects\MCD\Fringe Field Testing\11102021_Blank-400-800nm_PMT_cuvette')

    material = 'Au-PEG'

    for sig, blank in zip(signal_data, blank_data):
        print()




