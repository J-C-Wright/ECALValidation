# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.ticker import NullFormatter
from   matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as dates
import pandas as pd
import glob, os, re
import datetime
import re
import string
from math import sqrt
import matplotlib.cm as cm

region_labels = {
    "EB"      : "EB",
    "EB-gold" : "EB $R_{9} > 0.94$",
    "EB-bad"  : "EB $R_{9} < 0.94$",
    "EB-absEta_0_1" : "EB $\mid\eta\mid < 1$",
    "EB-absEta_1_1.4442" : "EB $1 < \mid\eta\mid < 1.4442$",
    "EE"      : "EE",
    "EE-gold" : "EE $R_{9} > 0.94$",
    "EE-bad" : "EE $R_{9} < 0.94$",
    "EE-absEta_1.566_2" : "EE $1.566 < \mid\eta\mid < 2$",
    "EE-absEta_2_2.5" : "EE $2 < \mid\eta\mid < 2.5$",
    "inclusive" : "All $\mid\eta\mid$ and $R_{9}$",
    "gold" : "All $\mid\eta\mid$  $R_{9} > 0.94$",
    "bad" : "All $\mid\eta\mid$ $R_{9} < 0.94$"
}

var_labels = {
    "DeltaM_data"           : "$\Delta m$",
    "DeltaP"                : "$\Delta p$",
    "width_data"            : "$\sigma_{CB}$",
    "rescaledWidth_data"    : "$\sigma_{CB}$ (Rescaled)",
    "additionalSmearing"    : "Additional Smearing",
    "chi2data"              : "$\chi^{2}$",
    "events_lumi"           : "events/lumi",
    "sigmaeff_data"         : "$\sigma_{eff}$",
    "sigmaeff_data_thirty"  : "$\sigma_{eff30}$",
    "sigmaeff_data_fifty"   : "$\sigma_{eff50}$",
    "nPV"                   : "$N_{pv}$"
}

var_ranges = {
    "DeltaM_data"           : [-4,4],
    "DeltaP"                : [-8,8],
    "width_data"            : [ 0,8],
    "rescaledWidth_data"    : [ 0,8],
    "additionalSmearing"    : [ 0,8],
    "chi2data"              : [ 0,8],
    "events_lumi"           : [-8,8],
    "sigmaeff_data"         : [-8,8],
    "sigmaeff_data_thirty"  : [-8,8],
    "sigmaeff_data_fifty"   : [-8,8],
    "nPV"                   : [-8,8]
}

format_fig_output = ['pdf','png']

def read_regions_from_table(path = "",tableFile= "",xVar=""):

    regions = []
    with open(path+'/'+tableFile) as f:
        for line in f.read().split('\n'):
            if line == '': continue
            if line[0] == '#': continue
            if 'category' in line: continue

            region = line.split('&')[0]
            if xVar == '':
                if region.split('-runNumber')[0] not in regions:
                    regions.append(region.split('-runNumber')[0])
            else:
                part = re.findall('[%s]+' % string.ascii_letters,region.replace(xVar,''))
                if len(part) == 0:
                    part = 'inclusive'
                    if part not in regions:
                        regions.append(part)
                else:
                    if part[0] not in regions:
                        regions.append(part[0])

    return regions


def read_regions_from_regionsfile(path = "",regionsFile=""):

    assert regionsFile != "", "The regions file must be specified"

    regions = []
    with open(path+regionsFile) as f:
        for line in f.read().split('\n'):
            if line != '' and line[0] != '#':
                regions.append(line.split()[0])

    return regions

def get_tables_from_config(path = '',config = ''):

    tables = []
    runRanges = []
    names = []
    colours = []

    with open(path+config) as f:
        for line in f.read().split('\n'):
            if line != '' and line[0] != '#':
                tables.append(line.split()[0])
                runRanges.append(line.split()[1])
                names.append(line.split()[2])
                colours.append(line.split()[3])

    return [names,runRanges,tables,colours]

def parse_table_over_regions(path = "", tableFile = "",category="",xVar=""):

    assert xVar != "", "The x variable name must be specified!"

    if category == "inclusive":
        cats = ['EB','EE','gold','bad'] #if inclusive specify the general eta and R9 regions to check against
    else:
        cats = [category]

    lines = []
    varmins = []
    varmaxes = []
    varmids = []

    text = open(path+tableFile).read()
    text = re.sub('\\pm','~',text)
    text = re.sub(r'[\\£@#${}^]','',text)

    assert xVar in text, "The specified variable label must be in the stability file!"

    #Parse top line for variables
    variables = text.split("\n")[0].split(' & ')

    for line in text.split('\n'):

        if line == '': continue
        if line[0] == '#': continue
        if 'category' in line: continue

        if category == 'inclusive':
            if not any(c in line.split('&')[0] for c in cats):
                lines.append(line)
                varmins.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[0]))
                varmaxes.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[1]))
                varmids.append((float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[1]) 
                                + float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[0]))/2.0)
        else:
            if category in line.split('&')[0]:
                lines.append(line)
                varmins.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[0]))
                varmaxes.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[1]))
                varmids.append((float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[1]) + float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[0]))/2.0)

    #If there's nothing for this category, return an empty dataframe
    if len(lines) == 0: return pd.DataFrame()
    varmins  = np.array(varmins)
    varmaxes = np.array(varmaxes)
    data = pd.DataFrame()
    data[xVar + '_min'] = varmins
    data[xVar + '_max'] = varmaxes
    data[xVar + '_mid'] = varmids
    data[xVar + '_err'] = (varmaxes - varmins)/2.0
    data[xVar ]         = (varmaxes + varmins)/2.0
    variables = [xVar+'min',xVar+'max',xVar+'mid'] + variables[2:]

    for i in range(3,len(variables),1):

        values = []
        errors = []
        for line in lines:
            value = line.split('&')[i].replace(' ','')

            if '--' in value: value = '0'

            if len(value.split('~')) == 2:
                values.append(float(value.split('~')[0]))
                errors.append(float(value.split('~')[1]))
            elif len(value.split('~')) == 1:
                values.append(float(value))
                errors.append(0.0)

            if 'mc' in variables[i]:
                variables[i] = variables[i].replace("mc","MC")
            if '/' in variables[i]:
                variables[i] = variables[i].replace("/","_")

        data[variables[i]] = values
        data[variables[i]+'_err'] = errors

    return data


def read_run_range(path = "", file = ""):

    raw_    = {
        'run_number': [],
        'Nevents'   : [],
        'UnixTime'  : []
    }
    data_ = pd.read_csv(path+'/'+file,
                sep="\t",
                names = ["run_number","Nevents","UnixTime"])
    #Transform the entries
    data_.Nevents = [float(x) for x in data_.Nevents]
    data_[ "UnixTime_min" ] = [ int(x.split('-')[0]) for x in data_.UnixTime ]
    data_[ "UnixTime_max" ] = [ int(x.split('-')[1]) for x in data_.UnixTime ]
    data_[ "run_min" ] = [ int(x.split('-')[0]) for x in data_.run_number ]
    data_[ "run_max" ] = [ int(x.split('-')[1]) for x in data_.run_number ]
    data_[ "date_min"] = [ datetime.datetime.fromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S') for x in data_.UnixTime_min ]
    data_[ "date_max"] = [ datetime.datetime.fromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S') for x in data_.UnixTime_max ]
    data_['time'] = pd.to_datetime(data_['date_max'], format='%Y-%m-%d')

    return data_

def append_variables(path='',file='',data=None,category=''):

    data_ = data

    text = open(path+'/'+file).read()
    text = re.sub('\\pm','~',text)
    text = re.sub(r'[\\£@#${}^]','',text)



    #Parse top line for variables
    variables = text.split("\n")[0].split(' & ')

    for i in range(2,len(variables),1):
        values = []
        errors = []
        for line in text.split('\n'):

            if line.split('-runNumber')[0] == category:

                value = line.split(' & ')[i]

                if '--' in value: value = '0'

                if len(value.split('~')) == 2:
                    values.append(float(value.split('~')[0]))
                    errors.append(float(value.split('~')[1]))
                elif len(value.split('~')) == 1:
                    values.append(float(value))
                    errors.append(0.0)

        if 'mc' in variables[i]:
            variables[i] = variables[i].replace("mc","MC")
        if '/' in variables[i]:
            variables[i] = variables[i].replace("/","_")

        data_[variables[i]]        = values
        data_[variables[i]+'_err'] = errors

    return data_

def draw_iov(ax, xData = None, iovs=[]):
    for v in iovs:
        for j in range(0, len(xData)-1):
            if v >= xData[j] and v < xData[j+1]:
                print "iov at ==", xData[j]
                ax.axvline(x=xData[j], color='red')

def plot_stability( xData = None, xData_err=None,
                    data_datasets = None, mc_datasets = None,
                    data_errorsets = None, mc_errorsets = None,
                    label = '', category = '', path = "",
                    evenX = False, xVar = '', iovs = None, 
                    showMC = False,names = None,style='ggplot',
                    oldStyle=False,colours=None):

    plt.style.use(style)

    left, width    = 0.1, 1.0
    bottom, height = 0.1, 0.5
    rect_hist = [left+width+0.01, bottom, 0.2, height]
    rect_plot = [left, bottom, width, height]

    nullfmt = NullFormatter()

    fig = plt.figure()

    ax_plot = plt.axes(rect_plot)
    ax_hist = plt.axes(rect_hist)
    for k, spine in ax_plot.spines.items():
        spine.set_zorder(10)
    for k, spine in ax_hist.spines.items():
        spine.set_zorder(10)

    ax_hist.yaxis.set_major_formatter(nullfmt)

    xPlaceholder = range(1,1+len(xData),1)
 
    for data_dataset,data_errorset,colour in zip(data_datasets,data_errorsets,colours):
        if evenX:
            if oldStyle:
                ax_plot.errorbar(xPlaceholder,data_dataset,yerr=data_errorset,capthick=0,marker='o',ms=4,c=colour,ls='none')
            else:
                ax_plot.errorbar(xPlaceholder,data_dataset,yerr=data_errorset,capthick=0,marker='.',ms=4,ls='solid',c=colour)
        else:
            if oldStyle:
                ax_plot.errorbar(xData,data_dataset,yerr=data_errorset,capthick=0,marker='o',ms=4,c=colour,ls='none')
            else:
                ax_plot.errorbar(xData,data_dataset,yerr=data_errorset,capthick=0,marker='.',ms=4,ls='solid',c=colour)

    # customise the axes
    xDataVar = xData.name
    if xDataVar == 'time':
        ax_plot.xaxis.set_minor_locator(dates.DayLocator(interval=10))
        ax_plot.xaxis.set_minor_formatter(dates.DateFormatter('%d\n%a'))
        ax_plot.xaxis.set_major_locator(dates.MonthLocator())
        ax_plot.xaxis.set_major_formatter(dates.DateFormatter('\n\n\n%b\n%Y'))
    elif (xDataVar == 'run_max' or xDataVar == 'run_min') and not evenX:
        majorLocator = MultipleLocator(125)
        minorLocator = MultipleLocator(62.5)
        ax_plot.xaxis.set_major_locator(majorLocator)
        ax_plot.xaxis.set_minor_locator(minorLocator)
        majorFormatter = FormatStrFormatter('%d')
        ax_plot.xaxis.set_major_formatter(majorFormatter)
        xlabels = ax_plot.get_xticklabels()
        plt.setp(xlabels, rotation=90, fontsize=10)
    elif (xDataVar == 'run_max' or xDataVar == 'run_min') and evenX:
        majorLocator = MultipleLocator(2)
        minorLocator = MultipleLocator(1)
        ax_plot.xaxis.set_major_locator(majorLocator)
        ax_plot.xaxis.set_minor_locator(minorLocator)
        xlabels = ax_plot.get_xticks().tolist()

        for i in range(2,len(xData),2):
            xlabels[i/2+1] = xData.tolist()[i-1]
        for i in range(len(xlabels)):
            if xlabels[i] < 200000: xlabels[i] = ''
        for index, row in iovs.iterrows():
            for j in range(0, len(xData)-1):
                if row['run'] >= xData[j] and row['run'] < xData[j+1]:
                    if 'Run' in row['info'] :
                        ax_plot.axvline(x=j+1, color='#EBAF3C',ls='--', zorder=5)
                    else :
                        ax_plot.axvline(x=j+1, color='#1072C3',ls='--', zorder=5)
        ax_plot.set_xticklabels(xlabels)
        xlabels = ax_plot.get_xticklabels()
        plt.setp(xlabels, rotation=90, fontsize=10)

    ax_plot.xaxis.grid(True, which="minor")
    ax_plot.yaxis.grid()
    ax_hist.xaxis.grid(True, which="minor")
    ax_hist.yaxis.grid()
    ax_hist.xaxis.set_ticks([])

    #Get and set the limits for the histogram
    for data_dataset,mc_dataset in zip(data_datasets,mc_datasets):
        if (len(mc_dataset) > 0):
            ymin = round(min(data_dataset.min(),mc_dataset.min())) - 1
            ymax = round(max(data_dataset.max(),mc_dataset.max())) + 1
        else:
            ymin = round(data_dataset.min()) - 1
            ymax = round(data_dataset.max()) + 1

    ax_plot.set_ylim((ymin,ymax))
    ax_hist.set_ylim((ymin,ymax))

    ax_plot.set_ylabel(label)

    nbin = 30

    hmaxes = []
    for data_dataset,colour in zip(data_datasets,colours):
        y,_,_ = ax_hist.hist(data_dataset, bins=nbin,orientation='horizontal', histtype='stepfilled', alpha=0.5,color=colour)
        hmaxes.append(y.max())
    hmax = max(hmaxes)
    ax_hist.set_xlim((0,hmax*1.1))

    ax_plot.set_title(region_labels[category] + "    " + label)

    #Annotate with mean and std dev
    step = 0.05
    for data_dataset,name in zip(data_datasets,names):
        npVals = np.asarray(data_dataset)
        name_gap = ' '*len(name)
        ax_hist.annotate('{:s}:\n$\mu$ = {:3.3f}, $\sigma$ = {:3.3f}'.format(name,np.mean(npVals),np.std(npVals)),(hmax/6,ymin-(ymax-ymin)*(0.1+step)),fontsize=11,annotation_clip=False,xycoords='data')
        step += 0.135
    
    #Add line for the MC 
    mc_dataset = mc_datasets[0]
    mc_errorset = mc_errorsets[0]
    if (len(mc_dataset) > 0):
        if evenX:
            if oldStyle:
                ax_plot.errorbar(xPlaceholder,mc_dataset,yerr=mc_errorset,capthick=0,marker='o',ms=4,ls='solid',c='Black')
            else:
                ax_plot.errorbar(xPlaceholder,mc_dataset,yerr=mc_errorset,capthick=0,marker='',ms=4,ls='dashed',c='Black')
        else:
            if oldStyle:
                ax_plot.errorbar(xData,mc_dataset,yerr=mc_errorset,capthick=0,marker='o',ms=4,ls='solid',c='Black')
            else:
                ax_plot.errorbar(xData,mc_dataset,yerr=mc_errorset,capthick=0,marker='',ms=4,ls='dashed',c='Black')

        if evenX:
            xNP = np.asarray(xPlaceholder)
        else:
            xNP = np.asarray(xData.tolist())

        mcNP    = np.asarray(mc_dataset.tolist())
        mcErrNP = np.asarray(mc_errorset.tolist())

        ax_plot.fill_between(xNP,mcNP-mcErrNP,mcNP+mcErrNP,alpha=0.5)

        if xVar == '':
            ax_hist.annotate('MC = {:3.3f} $\pm$ {:3.3f}'.format(mc_dataset[1],mc_errorset[1]),(hmax/6,ymin-(ymax-ymin)*(0.05+step)),fontsize=11,annotation_clip=False,xycoords='data')
        
    #Legend
    legend = ax_plot.legend(loc='lower left',numpoints=1,prop={'size':9})
    if (len(mc_datasets[0]) > 0):
        for i,name in enumerate(names):
            legend.get_texts()[i].set_text(name)
        legend.get_texts()[len(names)].set_text('MC')
    else:
        for i,name in enumerate(names):
            legend.get_texts()[i].set_text(name)
    legend.get_frame().set_alpha(0.5)
    legend.get_frame().set_linewidth(0.0)
    ax_hist.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
    ax_plot.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
    ax_plot.grid(which='minor', color='0.85', linestyle='--',dashes=(5,1),zorder=0)
    
    #Save
    if evenX:
        xDataVar = xDataVar + '_even'

    if not os.path.exists(path+'/'+xDataVar):
        os.makedirs(path+'/'+xDataVar)

    fileName = category + '_' + re.sub(r'[\\@#$/{}^]','', label)
    fileName = re.sub(r'[ ]','_',fileName)

    for fType in format_fig_output:
        print 'Saving plot: ' + path + xDataVar + '/' + fileName + '.' + fType
        plt.savefig( path + xDataVar + '/' + fileName+'.'+fType,
                     format=fType,orientation='landscape',
                     dpi=200,papertype='a4',pad_inches=0.1,
                     bbox_inches='tight')

    plt.close(fig)
