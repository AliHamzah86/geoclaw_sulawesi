
import os
import shutil
import subprocess

# create plots for all runs

driver_home = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(driver_home, '..', 'geoclaw_output')
setplot_ffname = os.path.join(driver_home, 'setplot.py')

def make_plot(dirname):
    if not os.path.basename(dirname).startswith('run_'):
        return

    shutil.copy2(setplot_ffname, os.path.join(dirname, 'setplot.py'))
    print(f'Running make plots in {dirname}')
    subprocess.run(['make', 'plots'], cwd=dirname, check=True)


for dirname, _, _ in os.walk(output_dir):
    make_plot(dirname)
