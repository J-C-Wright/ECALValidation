#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import sys, getopt
from stability  import ParseTable as pt 
from shutil     import copyfile
from optparse   import OptionParser
from pprint     import pprint

def get_options():
    parser = OptionParser()
    parser.add_option('-i', '--inv-mass',
                      dest='invMass', default='invMassSC',
                      help='''
                      Invariant mass variable.
                      ''')
    parser.add_option('-x','--xVar',
                      dest='xVar',default='',
                      help='''
                      Specify whether the plot should be over run range, or regions
                      ''')
    parser.add_option('-d','--outputName',
                      dest='outputName',default='',
                      help='''
                      Output directory for the plots
                      ''')
    parser.add_option('-c','--config',
                      dest='config',default='',
                      help='''
                      Config file containing runranges,tables,and names
                      ''')

    return parser.parse_args()
    
    
if __name__ == '__main__':

    (opt, args) =  get_options()

    assert opt.invMass     != '', 'No invariant mass name specified, use -i <invariant mass>'
    assert opt.outputName      != '', 'No output name specified'
    
    if opt.xVar != '':
        print 'x-axis variable is ', opt.xVar
    else:
        print 'x-axis variables are the times/runNumbers'
    
    if not os.path.exists('plot-stability/'):
        os.makedirs('plot-stability/')
    if not os.path.exists('data-stability/'):
        os.makedirs('data-stability/')
        
    
    config = opt.config.split('/')[-1]
    path = opt.config.split(config)[0]
    names,runRanges,tables = pt.get_tables_from_config(path = path,config = config)

    data_path = 'data-stability/' + opt.outputName + '/' + opt.invMass + '/' 
    plot_path = 'plot-stability/' + opt.outputName + '/' + opt.invMass + '/'

    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    if not os.path.exists(data_path):
        os.makedirs(data_path)


    #Get tables and runRanges from config + copy to data folder
    for table in tables:
        if table != data_path + table.split('/')[-1]:
            copyfile(table,data_path + table.split('/')[-1])

    if len(runRanges) != 0:
        for runRange in runRanges:
            if runRange != data_path + runRange.split('/')[-1]:
                copyfile(runRange,data_path + runRange.split('/')[-1])

    style = 'classic'

    #Reading regions from the table file
    regions = pt.read_regions_from_table(path=data_path,tableFile=tables[0].split('/')[-1],xVar=opt.xVar)
    print 'categories :: ', regions


    #Make plots
    print 'Starting plotmaking...'
    for region in regions:
        print 'Category: ',region

        #Prepare dataframe for each config entry
        dataFrames = []
        for table,runRange in zip(tables,runRanges):

            if opt.xVar != '':
                d=pt.parse_table_over_regions(path=data_path,tableFile=table,category=region,xVar=opt.xVar)
                dataFrames.append(d)
            else:
                #Get runrange and time info from the the runranges file
                d = pt.read_run_range(path=data_path,file=runRange.split('/')[-1])
                #Get variables information from the stability monitoring .tex file
                d = pt.append_variables(path=data_path,file=table.split('/')[-1],data=d,category=region)
                dataFrames.append(d)

        #Get variables to make plots of (data, not mc or err vars)
        variables = []
        if opt.xVar != '':
            xVars = [opt.xVar+'_min',opt.xVar+'_max',opt.xVar+'_mid']
        else:
            xVars = ['Nevents'     ,
                        'UnixTime'    ,
                        'run_number'  ,
                        'UnixTime_min',
                        'UnixTime_max',
                        'run_min'     ,
                        'run_max'     ,
                        'date_min'    ,
                        'date_max'    ,
                        'time'        ]

        for label in dataFrames[0].columns.values.tolist():
            if 'MC' not in label and label not in xVars and '_err' not in label:
                variables.append(label)

        #Loop over the vars
        for var in variables:
            #Get associated monte carlo info, or a placeholder
            varmc = var.replace('data','MC')

            mc_datasets = []
            mc_errorsets = []
            data_datasets = []
            data_errorsets = []

            for data in dataFrames:
                if 'MC' not in varmc:
                    print '[WARNING] MC counterpart not found for ', var
                    mc_datasets.append([])
                    mc_errorsets.append([])
                else:
                    mc_datasets.append(data[varmc])
                    mc_errorsets.append(data[varmc+'_err'])

                data_datasets.append(data[var])
                data_errorsets.append(data[var+'_err'])


            if opt.xVar == '':
                #Switches on whether the datapoints are evenly distributed along x
                evenXs = [False,True]
                #Plot as function of date or run numbers
                timevars = ['run_min','run_max','time']
                for timevar in timevars:
                    for evenX in evenXs:
                        pt.plot_stability( xData = dataFrames[0][timevar], data_datasets = data_datasets,
                                           data_errorsets = data_errorsets, mc_datasets = mc_datasets,
                                           mc_errorsets = mc_errorsets, label = pt.var_labels[var],
                                           category = region, path=plot_path, evenX = evenX,
                                           xVar=opt.xVar,names=names,style=style)
            else:
                xvars = [opt.xVar+'_min',opt.xVar+'_max',opt.xVar+'_mid']
                for xvar in xvars:
                    print 'xVar: ' + xvar
                    pt.plot_stability( xData = dataFrames[0][xvar], data_datasets = data_datasets,
                                       data_errorsets = data_errorsets, mc_datasets = mc_datasets,
                                       mc_errorsets = mc_errorsets, label = pt.var_labels[var],
                                       category = region, path=plot_path, evenX = False,
                                       xVar=opt.xVar,names=names,style=style)






