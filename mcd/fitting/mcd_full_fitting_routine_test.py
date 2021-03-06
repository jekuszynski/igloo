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

def plot_mcd(dic,op='avg',x_axis='Energy (eV)'):
    fig,ax=plt.subplots(figsize=(4,2))
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(-10,10,11),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    fig.colorbar(sm,ticks=range(-10,11,2),label='H (T)') #make color bar based on H (T) for plot
    for df in dic.values():
        #Dr. Seaborn or: How I Learned to Stop Worrying and Love sns.lineplot. Such efficiency. Much wow.
        sns.lineplot(data=df,x='energy',y='mdeg', linewidth=0.6,
                    hue='field',hue_norm=(-10,10),
                    palette=sns.color_palette('coolwarm_r',as_cmap=True),
                    legend=None)
    if x_axis=='Energy (eV)':
        plt.xlabel(x_axis)
    if op=='raw':
        plt.title("Raw MCD")
    if op=='avg':
        plt.title("Averaged MCD")
    plt.ylabel('MCD (mdeg)')
    plt.xlim(1.55,.75)
    plt.style.use('seaborn-paper')
    plt.savefig(op + '_mcd',dpi=100,transparent=False,bbox_inches='tight')
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
    plt.xlim(1.55,.75)
    baseline = lines.Line2D(range(6),np.zeros(1),c='black',ls='--',lw=0.6) #draw baseline at 0T
    ax.add_line(baseline) #add baseline to plot
    plt.style.use('seaborn-paper')
    plt.savefig('diff_mcd',dpi=100,transparent=True,bbox_inches='tight')
    plt.show()

def parse_abs(path):
    d_abs = {}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "blank" or "ems" in str.lower(name): #remove any test files performed during data acqusition
                f = re.search('.*(?=_A)',name).group(0) #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                # print("Adding", f + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                d_abs[f] = df #send dataframe to dictionary
    df_abs=pd.DataFrame(data=(d_abs['Ems']['wavelength'], d_abs['Ems']['energy'], d_abs['Ems']['chpx'], d_abs['Blank']['chpx'])).transpose() #make dataframe from ems/blank dictionary
    df_abs.columns=['wavelength','energy','ems','blank'] #setup columns
    df_abs['absorbance']=(2-np.log10(100 * df_abs['ems'] / df_abs['blank'])) #calculate absorbance from emission and blank data
    df_abs['smoothed_absorbance']=signal.savgol_filter(df_abs['absorbance'],25,2) #smooth absorbance plot using Savitzky-Golay
    return df_abs

def plot_abs(df,op='smooth',x_axis='energy'):
    fig,ax=plt.subplots(figsize=(4,2))
    if x_axis=='energy':
        plt.xlabel('Energy (eV)')
    # if op=='raw':
    #     plt.title("Raw Absorbance")
    #     sns.lineplot(data=df,x='energy',y='absorbance',c='Black')
    if op=='smooth':
        plt.title("Smoothed Absorbance")
        sns.lineplot(data=df,x='energy',y='smoothed_absorbance',color='Black')
    plt.ylabel('Absorbance (a.u.)')
    plt.xlim(2.0,0.75)
    plt.style.use('seaborn-paper')
    plt.savefig(op + '_abs',dpi=100,transparent=True,bbox_inches='tight')
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
    plt.xlim(2.1,0.55)
    plt.scatter(x,y,s=1.3,c='Black')
    plt.plot(x,fit_L,c='Blue')
    plt.plot(x,fit_R,c='Red')
    plt.legend(('LCP','RCP','Raw'))

    plt.subplot(2,1,2)
    plt.ylabel('Absorbance (a.u.)')
    plt.xlabel('Energy (eV)')
    plt.xlim(2.1,0.55)
    plt.plot(x,fit_diff,c='Purple')
    plt.legend(('Simulated MCD'))
    plt.show()

    return fit_diff

def func(x,ev): #define simulated mcd function from absorbance spectrum
    coeff=poly.polyfit(df_abs['energy'],df_abs['absorbance'],9) #find polynomial coeffs from original absorption spectra
    LCP=poly.polyval(x+ev,coeff) #find y from +ev shifted LCP spectrum
    RCP=poly.polyval(x-ev,coeff) #find y from -ev shifted RCP spectrum
    return LCP-RCP #return y from LCP-RCP

def calc_effective_mass_and_plot(abs_fit,diff_dic):
    ev_list=[]
    m_list=[]
    for field in diff_dic.keys():
        xdata=diff_dic[field]['energy'] 
        ydata=diff_dic[field]['mdeg']
        ydata_normalized=ydata/np.max(np.absolute(ydata))
        popt,pcov = optimize.curve_fit(func,xdata,ydata_normalized,bounds=(0,1)) #lsf optimization to spit out ev

        ev=popt[0] #return minimzed ev to variable
        ev_list.append(ev) #add ev to list
        c=299792458 #speed of light (m/s)
        e=1.60217662E-19 #charge of electron (C)
        m_e=9.10938356E-31 #mass of electron (kg)
        w_c=c/(1240/(ev/1000)*(10**-9)) #cyclotron resonance frequency
        effective_mass=e/w_c/2/m_e/math.pi #effective mass (m*/m_e), removed field scaling for now
        m_list.append(effective_mass) #add m* to list

        fig,ax=plt.subplots(figsize=(2,4))
        plt.title(str(field) + 'T Fit')
        plt.ylabel('MCD (deltaA/A_max*B) (T^-1) (x 10^-3)')
        plt.xlabel('Energy (eV)')
        plt.xlim(1.55,0.75)
        plt.plot(xdata,ydata_normalized,label='experiment_data',c='Black')
        plt.plot(xdata,[func(x,*popt) for x in xdata],label='simulated_fit',c='Red')
        baseline = lines.Line2D(range(6),np.zeros(1),c='black',ls='--',lw=0.6) #draw baseline at 0T
        ax.add_line(baseline) #add baseline to plot
        plt.legend(loc='lower right')
        plt.text(1.2,-0.5,'%.3f meV\n%.3f m*' % (ev,effective_mass),fontweight='bold',bbox={'facecolor':'white','alpha':0.5,'pad':0.1})
        plt.savefig(str(field) + "T_fit",dpi=100,transparent=True,bbox_inches='tight')
    average_ev = np.mean(ev_list)
    std_dev_ev = np.std(ev_list)
    average_m = np.mean(m_list)
    std_dev_m = np.std(m_list)
    return average_ev, std_dev_ev, average_m, std_dev_m

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

    ##-----Functional Code Below-----##

'''parse all data files'''
#Change these pathways if using from GitHub.
raw_mcd_dic = parse_mcd("/mnt/c/users/roflc/Desktop/MCD 11-11-20/") #raw mcd data in dictionary
df_abs = parse_abs("/mnt/c/users/roflc/Desktop/Abs 11-11-20/") #calculated abs data in dataframe

'''perform mcd experimental data operations'''
plot_mcd(raw_mcd_dic,'raw') #plot raw experimental mcd data
df_avgs = calc_raw_avg_mcd(raw_mcd_dic) #
plot_mcd(df_avgs,'avg')
df_diff = calc_diff_mcd(df_avgs)
plot_diff_mcd(df_diff)

'''perform mcd difference calculation'''
raw_mcd_dic_blank = parse_mcd("") 
for field, df in raw_mcd_dic.items():
    diff_df = DataFrame()
    diff_df[field] = raw_mcd_dic[field][df]['deltaA'] - raw_mcd_dic_blank[field][df]['deltaA']
    print (diff_df)

'''perform absorbance simulation data fitting operations'''
plot_abs(df_abs)
fit_diff=plot_CP_diff(df_abs['energy'],df_abs['absorbance'])
average_ev, std_dev_ev, average_m, std_dev_m = calc_effective_mass_and_plot(fit_diff,df_diff)

'''write HTML file report'''
writeHTMLfile('mcd.html')
