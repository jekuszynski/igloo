import sys
import os
import re
import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly 
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as lines
import seaborn as sns
from scipy import signal,optimize
import math

plt.rcParams.update({'figure.max_open_warning': 0}) #Remove figure creation RuntimeWarning.

def parse_mcd(path):
    d={}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "test" not in str.lower(name): #remove any test files performed during data acqusition
                field_name = re.search('(.*)(?=T_)..',name).group(0) #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                f=field_name + str(num%3) #differentiate repeated scans at same field
                # print("Adding", f + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['field']=int(re.search('(.*)(?=T)',name).group(0)) #add column for field
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                d[f] = df #send dataframe to dictionary
    return d

def parse_mcd_tests(path):
    d={}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "test" in str.lower(name): #remove any test files performed during data acqusition
                field_name = re.search('(.*)(?=T_)..',name).group(0) #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                f=field_name + str(num%3) #differentiate repeated scans at same field
                # print("Adding", f + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['field']=int(re.search('(.*)(?=T)',name).group(0)) #add column for field
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                d[f] = df #send dataframe to dictionary
    return d

# uncomment below to test dict format
    # print(d['-4T_0']['field'])
    # sys.exit()

# old code - keeping for posterity
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # cmap = plt.get_cmap('coolwarm') #set colormap for pos/neg field plotting
    # norm = colors.Normalize(vmin=-10,vmax=10) #normalize colormap knowing that only -10T and 10T are possible
    # fig,ax=plt.subplots()
    # for name, df in d.items():
    #     x=d[name]['energy']
    #     y=d[name]['mdeg']
    #     c=d[name]['field']
    #     ax.scatter(x,y,label=name,color=cmap(norm(c.values)))
    #     ax.legend(ncol=5)

def plot_mcd(dic,op='avg',x_axis='Energy (eV)',title='[PH]',xdata='energy',ydata='mdeg'):
    fig,ax=plt.subplots(figsize=(4,2))
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(-10,10,11),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    fig.colorbar(sm,ticks=range(-10,11,2),label='H (T)') #make color bar based on H (T) for plot
    for df in dic.values():
        #Dr. Seaborn or: How I Learned to Stop Worrying and Love sns.lineplot. Such efficiency. Much wow.
        sns.lineplot(data=df,x=xdata,y=ydata, linewidth=0.6,
                    hue='field',hue_norm=(-10,10),
                    palette=sns.color_palette('coolwarm_r',as_cmap=True),
                    legend=None)
    if x_axis=='Energy (eV)':
        plt.xlabel(x_axis)
    if op=='raw':
        plt.title("Raw MCD " + title)
    if op=='avg':
        plt.title("Averaged MCD " + title)
    plt.ylabel('MCD (mdeg)')
    plt.xlim(3.2,.55)
    plt.ylim(-50,50)
    plt.style.use('seaborn-paper')
    plt.savefig(op + '_mcd_' + title,dpi=200,transparent=False,bbox_inches='tight')
    plt.show()

def calc_raw_avg_mcd(dic): #need to define this before finding the mcd difference
    df_avgs={}
    for name, df in dic.items():
        field = re.search('(.*)(?=T)',name).group(0) #set variable 'field' equal to int(field) from original dictionary
        if field not in df_avgs: 
            df_avgs[field] = pd.DataFrame() #if field is not in new dictionary, create an empty dataframe
        if field in df_avgs:
            df_concat=pd.concat([df_avgs[field], df]).groupby(['wavelength'], as_index=False) #concatenate field entry with new df entry, so long as field is matching
            df_avgs[field] = df_concat #update dictionary entry with newly concatenated one
        df_avgs[field]=df_avgs[field].mean() #take the average of the concatenated df
    return df_avgs

def calc_diff_mcd(dic,op='sub'):
    df_diff={}
    for name, df in dic.items():
        df_diff[name] = pd.DataFrame()
        if name == '0': 
            # if op=='add':
            #     df_diff[name] = df + dic[list(name)[list(name).index(name)+3]]
            # elif op=='sub':
            #     df_diff[name] = df - dic[list(name)[list(name).index(name)+3]]
            del df_diff[name] #placeholder for now. In future would like to plot 0T field difference, but previous function looks like it deletes a version of 0 so no pair to subtract.
            pass
        elif '-' not in name: #loop only positive half of dictionary
            if op=='add':
                df_diff[name] = df + dic['-' + name] #add positive and negative dictionary entries
                df_diff[name]['energy'] = df['energy'] #fix energy back to original values
                df_diff[name]['field'] = df['field'] #fix field back to original values
            elif op=='sub':
                df_diff[name] = df - dic['-' + name] #subtract positive and negative dictionary entries
                df_diff[name]['energy'] = df['energy'] #fix field back to original values
                df_diff[name]['field'] = df['field']
        else:
            del df_diff[name]
            continue
    return df_diff

def mcd_blank_subtraction(dic_mcd,dic_blank):
    df_blank_subtracted={}
    for name, df in dic_mcd.items():
        df_blank_subtracted[name] = df - dic_blank[name]
        df_blank_subtracted[name]['energy'] = df['energy'] #fix field back to original values
        df_blank_subtracted[name]['field'] = df['field']

def plot_diff_mcd(dic,op='avg',x_axis='Energy (eV)'):
    fig,ax=plt.subplots(figsize=(4,2))
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(0,10,6),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='Greys',norm=norm) 
    fig.colorbar(sm,ticks=range(0,11,2),label='H (T)') #make color bar based on H (T) for plot
    for df in dic.values():
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

def parse_abs(path):
    d_abs = {}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "blank" in str.lower(name): #select for blank and ems file
                dic_name="Blank"
                # print("Adding", dic_name + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                d_abs[dic_name] = df #send dataframe to dictionary
            if "ems" or "sample" in str.lower(name):
                dic_name="Ems"
                # print("Adding", dic_name + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                d_abs[dic_name] = df #send dataframe to dictionary
    df_abs=pd.DataFrame(data=(d_abs['Ems']['wavelength'], d_abs['Ems']['energy'], d_abs['Ems']['chpx'], d_abs['Blank']['chpx'])).transpose() #make dataframe from ems/blank dictionary
    df_abs.columns=['wavelength','energy','ems','blank'] #setup columns
    df_abs['absorbance']=(2-np.log10(100 * df_abs['ems'] / df_abs['blank'])) #calculate absorbance from emission and blank data
    df_abs['smoothed_absorbance']=signal.savgol_filter(df_abs['absorbance'],25,2) #smooth absorbance plot using Savitzky-Golay
    return df_abs

def plot_abs(df,op='smooth',x_axis='energy'):
    fig,ax=plt.subplots(figsize=(4,2))
    if x_axis=='energy':
        plt.xlabel('Energy (eV)')
    if op=='raw':
        plt.title("Raw Absorbance")
        sns.lineplot(data=df,x='energy',y='absorbance',color='Black')
    if op=='smooth':
        plt.title("Smoothed Absorbance")
        sns.lineplot(data=df,x='energy',y='smoothed_absorbance',color='Black')
    plt.ylabel('Absorbance (a.u.)')
    plt.xlim(4,.55)
    plt.style.use('seaborn-paper')
    plt.savefig(op + '_abs',dpi=200,transparent=False,bbox_inches='tight')
    plt.show()

def plot_CP_diff(x,y,ev=0.04): #function to visually show separation of LCP and RCP from base abs
    coeff_L=poly.polyfit([x+ev for x in x],y,9) #LCP poly fit coeffs
    coeff_R=poly.polyfit([x-ev for x in x],y,9) #RCP poly fit coeffs
    fit_L=poly.polyval(x,coeff_L) #LCP line fit
    fit_R=poly.polyval(x,coeff_R) #RCP line fit
    fit_diff=(fit_L-fit_R)/(np.max(x))*1000 #calculate LCP-RCP normalized to absorbance max. Will incorporate based on field later.
    # x = x.values.tolist()
    plt.figure(figsize=(6,6),dpi=80)

    plt.subplot(2,1,1)
    plt.ylabel('Absorbance (a.u.)')
    plt.xlim(3.2,.55)
    plt.scatter(x,y,s=1.3,c='Black')
    plt.plot(x,fit_L,c='Blue')
    plt.plot(x,fit_R,c='Red')
    plt.legend(('LCP','RCP','Raw'))

    plt.subplot(2,1,2)
    plt.ylabel('Absorbance (a.u.)')
    plt.xlabel('Energy (eV)')
    plt.xlim(3.0,.55)
    plt.plot(x,fit_diff,c='Purple')
    plt.legend(('Simulated MCD'))
    plt.savefig('Simulated MCD.png',dpi=200,transparent=False,bbox_inches='tight')
    plt.show()

    return fit_diff

def nm_to_eV(nm):
    eV = 1240 / nm
    return ["%.3f" % z for z in eV]

def func(x,ev,yoffset): #define simulated mcd function from absorbance spectrum
    coeffL=poly.polyfit(df_abs['energy']+ev,df_abs['absorbance'],9) #find polynomial coeffs from original absorption spectra
    coeffR=poly.polyfit(df_abs['energy']-ev,df_abs['absorbance'],9) #find polynomial coeffs from original absorption spectra
    LCP=poly.polyval(x,coeffL) #find y from +ev shifted LCP spectrum
    RCP=poly.polyval(x,coeffR) #find y from -ev shifted RCP spectrum
    return LCP-RCP #return y from LCP-RCP, Note: may need to flip this depending on spectrum

    # df_abs['energy'],df_abs['absorbance']

def calc_effective_mass_and_plot(abs_fit,diff_dic):
    ev_list=[]
    m_list=[]
    for field in diff_dic.keys():
        if field is not '0':
            xdata=diff_dic[field]['energy'] 
            ydata=diff_dic[field]['zero_subtracted'] / 32982  # divided by mdeg conversion to obtain deltaA
            
            # all_pem_channels_added_diff_df[name]['sub_mod_zero_subtracted'] ? 

            # Perhaps I'm using either wrong dictionary, wrong ydata, OR I'm just using the func() wrong. Maybe mimic like plotted above?

            ydata_normalized=ydata/(np.max(df_abs['absorbance']))
            ydata_normalized=np.nan_to_num(ydata_normalized, nan=0.0)
            popt,pcov = optimize.curve_fit(func,xdata,ydata_normalized,p0=[0.0005,0.1]) #lsf optimization to spit out zeeman split mev, guessing ~0.05 meV
            
            # print(popt)
            # print(pcov) #list of residuals
            ev=np.absolute(popt[0]) #return absolute val of minimzed ev to variable
            ev_list.append(ev*1000) #add ev to list
            c=299792458 #speed of light (m/s)
            e=1.60217662E-19 #charge of electron (C)
            m_e=9.10938356E-31 #mass of electron (kg)
            w_c=c/(1240/(ev)*(10**-9)) #cyclotron resonance frequency
            effective_mass=e/w_c/2/m_e/math.pi #effective mass (m*/m_e), removed field scaling for now
            m_list.append(effective_mass) #add m* to list

            fig,ax=plt.subplots(figsize=(2,4))
            plt.title(str(field) + 'T Fit')
            plt.ylabel('MCD (deltaA/A_max*B) (T^-1) (x 10^-3)')
            plt.xlabel('Energy (eV)')
            plt.xlim(3.2,.55)
            plt.plot(xdata,ydata_normalized,label='experiment_data',c='Black')
            plt.plot(xdata,[func(x,*popt) for x in xdata],label='simulated_fit',c='Red')
            baseline = lines.Line2D(range(6),np.zeros(1),c='black',ls='--',lw=0.6) #draw baseline at 0T
            ax.add_line(baseline) #add baseline to plot
            plt.legend(loc=0)
            # plt.text(1.2,-0.5,'%.3f meV\n%.3f m*' % (ev,effective_mass),fontweight='bold',bbox={'facecolor':'white','alpha':0.5,'pad':0.1})
            plt.text(0,0,'%.3f meV\n%.3f m*' % (ev*1000,effective_mass),fontweight='bold',bbox={'facecolor':'white','alpha':0.5,'pad':0.1})
            plt.savefig(str(field) + "T_fit",dpi=100,transparent=False,bbox_inches='tight')
    average_ev = np.mean(ev_list)
    std_dev_ev = np.std(ev_list)
    average_m = np.mean(m_list)
    std_dev_m = np.std(m_list)
    return average_ev, std_dev_ev, average_m, std_dev_m, ev_list, m_list

def openHTML(f,title):
    f.write("<!DOCTYPE html>\n")
    f.write("<html lang='en'>\n")
    f.write("<head>\n")
    f.write('<base target="_blank"/>\n')
    f.write("<title>%s</title>\n" % title)
    f.write("</head>\n")
    f.write("<body>\n")
    f.write("<h1>%s</h1>\n" % title)

def writeHTMLimage(f,title,imgpath):
	# f.write('<p>%s</p>\n' % title)
	f.write('<img src="%s" />\n' % imgpath)

def writeHTMLspacer(f,spacer):
    f.write('%s' % spacer)

def closeHTML(f):
	f.write("</body>\n")
	f.write("</html>\n")
	f.close()	

def writeHTMLfile(file_name):
    f=open(file_name,'w')
    openHTML(f,'MCD 20201111 Report')
    writeHTMLspacer(f,'<div>\n')
    writeHTMLimage(f,'raw_mcd','raw_mcd.png')
    writeHTMLimage(f,'avg_mcd','avg_mcd.png')
    writeHTMLimage(f,'diff_mcd','diff_mcd.png')
    writeHTMLimage(f,'smooth_abs','smooth_abs.png')
    writeHTMLspacer(f,'</div>\n<div>')
    #These next few most certainly deserve a loop... I'll get around to it eventually...
    writeHTMLimage(f,'2T_fit','2T_fit.png')
    writeHTMLimage(f,'4T_fit','4T_fit.png')
    writeHTMLimage(f,'6T_fit','6T_fit.png')
    writeHTMLimage(f,'8T_fit','8T_fit.png')
    writeHTMLimage(f,'10T_fit','10T_fit.png')
    writeHTMLspacer(f,'</div>\n')
    f.write('<p><u>From the above data:</u></p>')
    f.write('<p>The average Zeeman splitting energy is <b>%.3f</b> \u00B1 %.4f meV.</p>' % (average_ev,std_dev_ev))
    f.write('<p>The average effective mass (<i>m*</i>) is  <b>%.3f</b> \u00B1 %.4f.</p>' % (average_m,std_dev_m))
    closeHTML(f)

def writeHTMLfile_difference(file_name,report_date):
    f=open(file_name,'w')
    openHTML(f,'MCD ' + report_date + ' Report')
    writeHTMLspacer(f,'<div>\n')
    f.write('<p><b>Raw MCD Spectra</b></p>')
    writeHTMLimage(f,'raw_mcd','raw_mcd_sample.png')
    writeHTMLimage(f,'raw_mcd','raw_mcd_blank.png')
    writeHTMLspacer(f,'</div>\n<div>')
    f.write('<p><b>Average MCD Spectra</b></p>')
    writeHTMLimage(f,'avg_mcd','avg_mcd_sample.png')
    writeHTMLimage(f,'avg_mcd','avg_mcd_blank.png')
    writeHTMLspacer(f,'</div>\n<div>')
    f.write('<p><b>Diff MCD Spectra & Absorbance Spectra</b></p>')
    writeHTMLimage(f,'S-B_diff_mcd','avg_mcd_diff.png')
    writeHTMLimage(f,'S-0T_diff_mcd','avg_mcd_Diff_no_blank_0T_subbed.png')
    writeHTMLimage(f,'S-B-0T','avg_mcd_diff_0T_subbed.png') 
    writeHTMLimage(f,'raw_abs','raw_abs.png')
    writeHTMLimage(f,'smooth_abs','smooth_abs.png')
    writeHTMLspacer(f,'</div>\n<div>')
    f.write('<p><b>Diff Modulus MCD Spectra</b></p>')
    writeHTMLimage(f,'sample_mcd_modulus','avg_mcd_sample_modulus.png')
    writeHTMLimage(f,'blank_mcd_modulus','avg_mcd_blank_modulus.png')
    writeHTMLimage(f,'sample-0T','avg_mcd_sample-0T.png')  
    writeHTMLimage(f,'diff_mcd_modulus','avg_mcd_sub_modulus.png')
    writeHTMLimage(f,'diff_mcd_modulus_0T_subtracted','avg_mcd_sub_modulus_zero_subtracted.png')
    writeHTMLspacer(f,'</div>\n<div>')
    #These next few most certainly deserve a loop... I'll get around to it eventually...
    writeHTMLimage(f,'2T_fit','2T_fit.png')
    writeHTMLimage(f,'4T_fit','4T_fit.png')
    writeHTMLimage(f,'6T_fit','6T_fit.png')
    writeHTMLimage(f,'8T_fit','8T_fit.png')
    writeHTMLimage(f,'10T_fit','10T_fit.png')
    writeHTMLspacer(f,'</div>\n')
    f.write('<p><u>From the above data:</u></p>')
    f.write('<p>The average Zeeman splitting energy is <b>%.3f</b> \u00B1 %.4f meV.</p>' % (average_ev,std_dev_ev))
    f.write('<p>The average effective mass (<i>m*</i>) is  <b>%.3f</b> \u00B1 %.4f.</p>' % (average_m,std_dev_m))
    closeHTML(f)



'''-------------------------------FUNCTIONAL CODE BELOW-------------------------------'''


'''parse all data files'''
#Change these pathways if using from GitHub.
raw_mcd_dic = parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/MCD 04-02-21/") #raw mcd data in dictionary
df_abs = parse_abs("/mnt/c/Users/roflc/Desktop/MCD DATA/ABS 03-29-21/") #calculated abs data in dataframe

# raw_mcd_dic = parse_mcd("") #USB
# df_abs = parse_abs("") #USB

'''fit raw and avg mcd straight from datafile - no workup'''
plot_mcd(raw_mcd_dic,'raw',title='sample') #plot raw experimental mcd data
df_avgs = calc_raw_avg_mcd(raw_mcd_dic) #
plot_mcd(df_avgs,'avg',title='sample')

'''mcd difference (no blank)'''
for name, df in df_avgs.items():
    df_avgs[name]['vac_mod'] = np.sqrt((df_avgs[name]['pemx'])**2 + (df_avgs[name]['pemy'])**2)
    df_avgs[name]['vdc_mod'] = np.sqrt((df_avgs[name]['chpx'])**2 + (df_avgs[name]['chpy'])**2)
    df_avgs[name]['total_mod_test'] = df_avgs[name]['vac_mod'] / df_avgs[name]['vdc_mod'] * 32982
for name, df in df_avgs.items():
    df_avgs[name]['0T_total_mod_sub'] = df_avgs[name]['total_mod_test'] - df_avgs['0']['total_mod_test']
    df_avgs[name]['avg-0T'] = df_avgs[name]['mdeg'] - df_avgs['0']['mdeg']
plot_mcd(df_avgs,'avg',title='total_mod_test',ydata='total_mod_test')
plot_mcd(df_avgs,'avg',title='0T_total_mod_sub',ydata='0T_total_mod_sub')
plot_mcd(df_avgs,'avg',title='Diff_no_blank_0T_subbed',ydata='avg-0T')

'''mcd difference with blank'''
raw_mcd_dic_blank = parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/MCD 03-30-21 Blank/")
plot_mcd(raw_mcd_dic_blank,'raw',title='blank')
df_blank_avgs = calc_raw_avg_mcd(raw_mcd_dic_blank)
plot_mcd(df_blank_avgs,'avg',title='blank')

# make this a function, finds diff between sample and blank only.
diff_df={}
for name, df in df_avgs.items():
    for df_blank in df_blank_avgs.values():
        diff_df[name] = pd.DataFrame()
        diff_df[name] = df - df_blank #take difference of all values
        diff_df[name]['energy'] = df['energy'] #fix energy
        diff_df[name]['field'] = df['field'] #fix field
        diff_df[name]['wavelength'] = df['wavelength'] #fix wavelength
for name, df in df_avgs.items():
    for df_blank in df_blank_avgs.values():
        diff_df[name]['zero_subtracted'] = diff_df[name]['mdeg'] - diff_df['0']['mdeg']
plot_mcd(diff_df,'avg',title='diff',ydata='mdeg')
plot_mcd(diff_df,'avg',title='diff_0T_subbed',ydata='zero_subtracted')

# make this a function, finds diff between sample and blank using total x/y signal.
all_pem_channels_added_diff_df={}
for name, df in df_avgs.items():
    for df_blank in df_blank_avgs.values():
        all_pem_channels_added_diff_df[name] = pd.DataFrame()
        all_pem_channels_added_diff_df[name] = df - df_blank #take difference of all values
        all_pem_channels_added_diff_df[name]['energy'] = df['energy'] #fix energy
        all_pem_channels_added_diff_df[name]['field'] = df['field'] #fix field
        all_pem_channels_added_diff_df[name]['wavelength'] = df['wavelength'] #fix wavelength
        all_pem_channels_added_diff_df[name]['sample_mod'] = np.sqrt(df['pemx']**2 + df['pemy']**2) / np.sqrt(df['chpx']**2 + df['chpy']**2) * 32982
        all_pem_channels_added_diff_df[name]['blank_mod'] = np.sqrt(df_blank['pemx']**2 + df_blank['pemy']**2) / np.sqrt(df_blank['chpx']**2 + df_blank['chpy']**2) * 32982
for name, df in df_avgs.items():
    for df_blank in df_blank_avgs.values():
        all_pem_channels_added_diff_df[name]['sample-0T'] = (all_pem_channels_added_diff_df[name]['sample_mod'] - all_pem_channels_added_diff_df['0']['sample_mod']) 
        all_pem_channels_added_diff_df[name]['modulus_subtracted'] = all_pem_channels_added_diff_df[name]['sample_mod'] - all_pem_channels_added_diff_df[name]['blank_mod'] 
for name, df in df_avgs.items():
    for df_blank in df_blank_avgs.values():
        all_pem_channels_added_diff_df[name]['sub_mod_zero_subtracted'] = (all_pem_channels_added_diff_df[name]['modulus_subtracted'] - all_pem_channels_added_diff_df['0']['modulus_subtracted']) 
plot_mcd(all_pem_channels_added_diff_df,'avg',title='sample_modulus',ydata='sample_mod')
plot_mcd(all_pem_channels_added_diff_df,'avg',title='blank_modulus',ydata='blank_mod')
plot_mcd(all_pem_channels_added_diff_df,'avg',title='sample-0T',ydata='sample-0T')
plot_mcd(all_pem_channels_added_diff_df,'avg',title='sub_modulus',ydata='modulus_subtracted')
plot_mcd(all_pem_channels_added_diff_df,'avg',title='sub_modulus_zero_subtracted',ydata='sub_mod_zero_subtracted')

# '''perform absorbance simulation data fitting operations'''
plot_abs(df_abs,op='raw')
plot_abs(df_abs)

fit_diff=plot_CP_diff(df_abs['energy'],df_abs['absorbance'])
average_ev, std_dev_ev, average_m, std_dev_m, ev_list, m_list = calc_effective_mass_and_plot(fit_diff,diff_df)

# '''write HTML file report'''
# # writeHTMLfile('mcd.html','11-11-2020')
writeHTMLfile_difference('mcd_difference.html','03-30-2021, Both Max Signals')

print("...\nDone!")