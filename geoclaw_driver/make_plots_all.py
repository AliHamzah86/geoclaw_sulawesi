
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
    makefile_path = os.path.join(dirname, 'Makefile')
    if not os.path.isfile(makefile_path):
        # Skip aggregated directories (e.g. all_runs_npy_files) that lack a Makefile
        return
    output_dir = os.path.join(dirname, '_output')
    if not os.path.isdir(output_dir):
        return

    shutil.copy2(setplot_ffname, os.path.join(dirname, 'setplot.py'))
    print(f'Running make plots in {dirname}')
    try:
        subprocess.run(['make', 'plots'], cwd=dirname, check=True)
    except subprocess.CalledProcessError as err:
        print(f"[WARN] make plots failed in {dirname}: {err}")


for dirname, _, _ in os.walk(output_dir):
    make_plot(dirname)
