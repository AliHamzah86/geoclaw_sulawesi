
import os,sys

# create plots for all runs

driver_home = os.getcwd()
output_dir = '../geoclaw_output'
setplot_ffname = os.path.join(redclaw_home,'setplot.py')

def make_plot(arg, dirname, fnames):

    if 'run_' in dirname:
        cmd_list = ['cd ' + dirname
                    'cp ' + setplot_ffname + ' .',\
                    'make plots']
        cmd = '; '.join(cmd_list)
        print('running: ' + cmd)
        #os.system(cmd)


os.path.walk(output_dir,make_plot,[])

