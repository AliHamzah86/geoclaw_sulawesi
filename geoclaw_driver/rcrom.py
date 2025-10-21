

import os
import shutil
import numpy as np
from clawpack.geoclaw import dtopotools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import pickle
import time, datetime


# ==============================================================================
#  Drom Class
# ==============================================================================

class Drom(object):
    r"""
    Delaunay-reversal Reduced Order Model

    

    """

    def __init__(self, path_common='../geoclaw_output/common', input_type='slip', drom_type=None):
        r"""
        initialize class

            * set common path (defaults to 'common')
            * initialize GeoClawInput class containing run info
            * set current run_id to zero
            * set pickle filename
            * set default initial time t0
            * node_id to run_id maps (and vice versa)
        
        r"""

        self.path_common = path_common      # set common path
        self.GeoClawInput = GeoClawInput()
        self.GeoClawInput.path_common = self.path_common 
        self.input_list = []
        
        self._continue = False              # continue previous run
        self._input_type = input_type       # input param type (TODO: remove)
        
        self._run_id = 0                    # current run_id
        self._run_id_list = []              # set list of run_id's
        self._runerror_id_list = []         # keep track of unfinished runs
       
        self._pickle_fname = None           # pickle filename
        
        self._t0 = None                     # set default t0

        self._delaunay = None               # delaunay tessellation
        self._run_id2node_id = None         # map run_id to node_id
        self._node_id2run_id = None         # map node_id to run_id



    def _check_run(self):
        r"""
        check previous run with rundir

            *if* the old run has all the proper count of outputs fort.q????
            then return True

            *if* old run has check-pointing then include it in 
            self._runerror_id_list then return True

            *else* return False
        
        r"""

        value = False     # default return value is False
        
        #rundir = 'run_' + str(run_id)
        rundir = self.GeoClawInput._rundir
        output_dir = os.path.join(rundir, '_output')

        # first check if rundir exists, if not, return False
        if not os.path.exists(output_dir):
            value = False
            return value
        
        
        # count the number of output files needed, depending on output_style
        output_style = self.GeoClawInput._rundata.clawdata.output_style
        if output_style == 1:
            noutput_times = \
                        self.GeoClawInput._rundata.clawdata.num_output_times
        elif output_style == 2:
            noutput_times = \
                        len(self.GeoClawInput._rundata.clawdata.output_times)
        elif output_style == 3:
            noutput_times = self.GeoClawInput._rundata.clawdata.total_steps
        
        fortt_fname_list = []
        fortq_fname_list = []

        # loop through all output files in the output directory
        output_fname_list = os.listdir(output_dir)
        for output_fname in output_fname_list:
            
            if 'fort.tck' in output_fname:
                # if checkpoint was found, return True and add run_id to
                # self._runerror_id_list
                print('(_check_run) this run was checkpointed: ' + str(run_id))
                if not run_id in self._runerror_id_list:
                    self._runerror_id_list.append(run_id)
                value = False
                return value
            
            else:
                if 'fort.t' in output_fname:
                    fortt_fname_list.append(output_fname)
                elif 'fort.q' in output_fname:
                    fortq_fname_list.append(output_fname)

        # sort filenames?
        #fortt_fname_list.sort()     
        #fortq_fname_list.sort()     

        # check fort.q* fort.t* count then return True if numbers are correct.
        nfortt = len(fortt_fname_list)
        nfortq = len(fortq_fname_list)
        
        value = ((nfortt == noutput_times) and (nfortq == noutput_times))
         
        return value



    def evaluate_hdm(self, save_drom=False, return_grb=False):
        """
        Run GeoClaw with current setup, writing input files,
        handling logs, and optionally returning gauge data.
        """
        import datetime, subprocess

        run_complete = False
        rundir = self.GeoClawInput._rundir
        run_id = self.GeoClawInput._run_id
        rundir_full = os.path.join(os.getcwd(), rundir)

        os.makedirs(rundir, exist_ok=True)

        if (self._check_run() and self._continue):
            print(f"(evaluate_hdm) old run exists, skipping run: {rundir_full}")
        else:
            # --- Write input files (includes dtopo, fgout, fgmax, plots) ---
            print("[DEBUG] Writing GeoClaw input files...")
            self.GeoClawInput.write_files()

            # --- Ensure xgeoclaw exists ---
            exe_source = os.path.join(self.path_common, "xgeoclaw")
            exe_makefile = os.path.join(self.path_common, "Makefile")
            if not os.path.exists(exe_source):
                os.makedirs(self.path_common, exist_ok=True)
                os.system(f"cp Makefile {self.path_common}")
                os.system(f"cd {self.path_common}; make new")

            # Symlink xgeoclaw into rundir
            exe_target = os.path.join(rundir, "xgeoclaw")
            if not os.path.exists(exe_target):
                os.system(f"ln -sf {os.path.relpath(exe_source, rundir)} {exe_target}")

            # --- Backup logs ---
            runlog_fname = os.path.join(rundir, "xgeoclaw.log")
            runerrlog_fname = os.path.join(rundir, "xgeoclaw.errlog")
            if os.path.exists(runlog_fname):
                os.rename(runlog_fname, runlog_fname + ".backup")
                if os.path.exists(runerrlog_fname):
                    os.rename(runerrlog_fname, runerrlog_fname + ".backup")

            # --- Timer ---
            timer_fname = os.path.join(rundir, "xgeoclaw.timer")
            tic = datetime.datetime.now()
            with open(timer_fname, "w") as timer_file:
                timer_file.write(f"started:\t{tic.strftime('%Y%m%d-%H%M%S')}\n")

            # --- Ensure dtopo.data is visible in _output/ where xgeoclaw runs ---
            output_dir = os.path.join(rundir, "_output")
            src_dtopo = os.path.join(rundir, "dtopo.data")
            dst_dtopo = os.path.join(output_dir, "dtopo.data")

            try:
                os.makedirs(output_dir, exist_ok=True)
                if os.path.exists(src_dtopo):
                    shutil.copy(src_dtopo, dst_dtopo)
                    print(f"[DEBUG] Copied dtopo.data → {dst_dtopo}")
                else:
                    print(f"[WARNING] Missing dtopo.data in {rundir}")
            except Exception as e:
                print(f"[ERROR] Could not copy dtopo.data: {e}")

            # --- Ensure Makefile is present in rundir ---
            makefile_src = os.path.join(self.path_common, "Makefile")
            makefile_dst = os.path.join(rundir, "Makefile")
            if not os.path.exists(makefile_dst):
                if os.path.exists(makefile_src):
                    shutil.copy(makefile_src, makefile_dst)
                    print(f"[DEBUG] Copied Makefile → {makefile_dst}")
                else:
                    print(f"[WARNING] Missing Makefile in {self.path_common}")

            # --- Run xgeoclaw ---
            print(f"[INFO] Starting xgeoclaw simulation in {rundir} >>>")

            import subprocess

            log_path = os.path.join(rundir, "xgeoclaw.log")
            err_path = os.path.join(rundir, "xgeoclaw.errlog")

            # Run GeoClaw cleanly, suppressing PETSc/MPI noise
            with open(log_path, "w") as out, open(err_path, "w") as err:
                ret = subprocess.run(
                    ["make", "output"],      # equivalent to "cd rundir; make output"
                    cwd=rundir,
                    stdout=out,
                    stderr=subprocess.DEVNULL   # hide PETSc cleanup messages
                ).returncode

            # --- Retry once if the run failed (e.g., dtopo.data missing) ---
            if ret != 0:
                output_dir = os.path.join(rundir, "_output")
                src_dtopo = os.path.join(rundir, "dtopo.data")
                dst_dtopo = os.path.join(output_dir, "dtopo.data")

                if os.path.exists(src_dtopo):
                    shutil.copy(src_dtopo, dst_dtopo)
                    print(f"[FIX] Re-copied dtopo.data → {dst_dtopo} after make cleaned _output/")
                    with open(log_path, "a") as out, open(err_path, "a") as err:
                        ret = subprocess.run(
                            ["make", "output"],
                            cwd=rundir,
                            stdout=out,
                            stderr=subprocess.DEVNULL
                        ).returncode
                else:
                    print(f"[WARNING] Could not retry: missing {src_dtopo}")

            # If Fortran failed due to missing dtopo.data, copy it again and retry once
            if ret != 0:
                output_dir = os.path.join(rundir, "_output")
                src_dtopo = os.path.join(rundir, "dtopo.data")
                dst_dtopo = os.path.join(output_dir, "dtopo.data")

                if os.path.exists(src_dtopo):
                    shutil.copy(src_dtopo, dst_dtopo)
                    print(f"[FIX] Re-copied dtopo.data → {dst_dtopo} after make cleaned _output/")
                    ret = os.system(f"cd {rundir}; make output 1> xgeoclaw.log 2> xgeoclaw.errlog")
                else:
                    print(f"[WARNING] Could not retry: missing {src_dtopo}")

            toc = datetime.datetime.now()
            with open(timer_fname, "a") as timer_file:
                      status = "OK" if ret == 0 else "FAILED"
                      timer_file.write(f"ended ({status}):\t{toc.strftime('%Y%m%d-%H%M%S')}\n")

            # --- Check completion ---
            run_complete = self._check_run()
            if run_complete:
                print(f"[INFO] Run completed successfully in {rundir}")
            else:
                print(f"[ERROR] xgeoclaw failed in {rundir}, check xgeoclaw.errlog")
                errlog = os.path.join(rundir, "xgeoclaw.errlog")
                if os.path.exists(errlog):
                    with open(errlog, "r") as f:
                        tail = f.readlines()[-10:]
                    print("[ERROR] Last 10 lines of xgeoclaw.errlog:")
                    for line in tail:
                        print(line.strip())
                return None  # early exit if failed
            
            # --- Check fort.amr for warnings/errors ---
            fort_amr = os.path.join(rundir, "_output", "fort.amr")
            if os.path.exists(fort_amr):
                with open(fort_amr, "r") as f:
                    tail = f.readlines()[-20:]
                print("[INFO] Tail of _output/fort.amr (last 20 lines):")
                for line in tail:
                    print(line.rstrip())
        
        # --- Update bookkeeping ---
        self._run_id_list.append(run_id)

        if save_drom:
            self.save_drom()

        if return_grb and run_complete:
            try:
                grb = self.read_gauges(run_id, rundir=rundir)  # pass correct path
                return grb
            except Exception as e:
                print(f"[ERROR] Could not read gauges: {e}")
                if run_id not in self._runerror_id_list:
                    self._runerror_id_list.append(run_id)
                return None

        def evaluate(self, m, response_type='gauge'):
            r"""
            evaluate reduced order model. 

            call self._delaunay.transform to find barycentric coordinates
            use self.read_gauge_response() and self._linear_comb()
            to compute gauge_response
                
            """

            if self._delaunay == None:
                # TODO: raise exception
                print('\t(evaluate) no delaunay tessellation found.')
                self.update_delaunay()

            # compute barycentric coordinates

            ndim = self.GeoClawInput.input_dim
            simplex_no = self._delaunay.find_simplex(m)
            nodes = self._delaunay.simplices[simplex_no]
            invT = self._delaunay.transform[simplex_no,:-1,:]

            r = self._delaunay.transform[simplex_no,-1,:]
            c = np.dot(invT , (np.array(m) - r.T))

            run_id_list = [self._node_id2run_id[j] for j in nodes]

            if response_type == 'gauge':
                grb_list = self.read_gauges_list(run_id_list)
                grb = self.linearc_grb(grb_list,c)

            elif response_type == 'eta_max':
                emb_list = self.read_eta_max_list(run_id_list)
                emb = self.linearc_emb(emb_list,c)

            return emb

    
    def linearc_grb(self,grb_list,c):
        r"""
        compute linear combination of grb's 

        TODO: if grb can have its own class binary operation
              this code will be simplified

        """
        grb = Gauge_response_bundle()
        grb.response_list = []

        n = len(grb_list)

        for j in range(n):
            m = len(grb_list[j].response_list)
            for k in range(m):
                t = grb_list[j].response_list[k].t
                self._update_t0(t)

        array_size = self._t0.shape[0]

        # compute redundant barycentric-coordinates
        c = np.concatenate((c,[1.- np.sum(c)])) 

        for i in range(m):

            gr = Gauge_response()

            gr.t = np.zeros(array_size)
            gr.h = np.zeros(array_size)
            gr.hu = np.zeros(array_size)
            gr.hv = np.zeros(array_size)
            gr.eta = np.zeros(array_size)

            gr.run_id = []

            for j in range(n):

                grb0 = grb_list[j]
                gr0 = grb0.response_list[i]
                gr1 = self.interp2t0(gr0)

                gr.gauge_id = gr1.gauge_id
                gr.x = gr1.x  
                gr.y = gr1.y  
                gr.t = gr1.t
                gr.h += gr1.h*c[j]
                gr.hu += gr1.hu*c[j]
                gr.hv += gr1.hv*c[j]
                gr.eta += gr1.eta*c[j]

                gr.run_id.append(gr1.run_id)

            grb.response_list.append(gr)

        return grb
    

    def linearc_emb(self,emb_list,c):
        r"""
        compute linear combination of grb's 

        TODO: if grb can have its own class binary operation
              this code will be simplified

        """
        emb = Eta_max_bundle()
        emb.response_list = []

        n = len(emb_list)
        m = len(emb_list[0].response_list)

        # compute redundant barycentric-coordinates
        c = np.concatenate((c,[1.- np.sum(c)])) 

        for i in range(m):

            em = Eta_max()

            em.t = None
            em.eta_max = 0.
            em.run_id = []

            for j,emb0 in enumerate(emb_list):

                em0 = emb0.response_list[i]

                em.gauge_id = em0.gauge_id
                em.x = em0.x  
                em.y = em0.y  
                em.eta_max += em0.eta_max*c[j]

                em.run_id.append(em0.run_id)

            emb.response_list.append(em)

        return emb
 
    def update_delaunay(self):
        r"""
        add point m to self._delaunay
        if tessellation isn't initialized, do so

        if self._delaunay is incremental, do .add_points,
        otherwise, recompute tessellation
        (wait - looks like Qhull has a bug. incremental mode for
        scipy.spatial is not to be trusted..)

        check to match with ._hash_id_list
        
        """
            
        if len(self.input_list) > self.GeoClawInput.input_dim:

            from scipy.spatial import Delaunay

            print('\t(update_delaunay) generating Delaunay')
            self._delaunay = Delaunay(self.input_list)
            # default indexing for _run_id2node_id
            self._run_id2node_id = range(len(self._run_id_list))
            self._node_id2run_id = range(len(self._run_id_list))
        else:
            # TODO: raise an exception
            print('\t(update_delaunay) error: need more eval pts')

    
    def update_delaunay_special(self, sub_run_id_list):
        r"""

        update delaunay with only subset of run points

        """

        # list() necessary?
        sub_input_list = [list(self.input_list[j]) for j in sub_run_id_list]

        if len(sub_input_list) > len(self.GeoClawInput.input_dim):

            from scipy.spatial import Delaunay

            print('\t(update_delaunay) generating Delaunay')
            self._delaunay = Delaunay(sub_input_list)

            # indexing for _run_id2node_id
            self._run_id2node_id = [-1] * len(self._run_id_list)
            self._node_id2run_id = []
            for j,run_id in enumerate(sub_run_id_list):
                self._run_id2node_id[run_id] = j 
                self._node_id2run_id.append(run_id)
        else:
            # TODO: raise an exception
            print('\t(update_delaunay) error: need more eval pts')
        

    def l2sqerror(self, m):
        r"""
        compute l2error difference between evaluate_hdm(m) and evaluate(m)
        then run self.update_delaunay(m)

        """

        grb1 = self.evaluate(m)
        grb2 = self.evaluate_hdm(m)

        # assume: grb1, grb2 already interp2t0
        l2sqerror_val = 0.
        
        for j in range(len(grb1.response_list)):
            l2sqerror_val += np.linalg.norm(grb1[j].eta - grb2[j].eta) ** 2

        return l2sqerror_val


    def interp2t0(self, gr, t0=None):
        r"""
        run self._update_t0(gauge_response)
        then interpolate all gauge values to time-tick t0
        
        INPUT gr in Gauge_response class
              t0 = ... set t0 manually

        """

        if t0 == None:
            t0 = self._t0

        igr = Gauge_response()

        igr.t = self._t0.copy()
        igr.x = gr.x 
        igr.y = gr.y 
        igr.t1 = gr.t1 
        igr.t2 = gr.t2
        igr.run_id = gr.run_id
        igr.gauge_id = gr.gauge_id

        igr.h = np.interp(t0, gr.t, gr.h)
        igr.hu = np.interp(t0, gr.t, gr.hu)
        igr.hv = np.interp(t0, gr.t, gr.hv)
        igr.eta = np.interp(t0, gr.t, gr.eta)

        return igr


    def _update_t0(self, t):
        r"""
        update self._t0 as necessary
        currently just sets up a uniform grid with smallest dt increment
        found when merging t and self._t0 (may be excessive?)


        INPUT t is the new time-grid
                should be a numpy array

        """
        

        if self._t0 == None:
            self._t0 = t.copy()

        elif len(self._t0) < 2:
            self._t0 = t.copy()

        else:
            tb = np.min([self._t0[0], t[0]])
            tf = np.max([self._t0[-1], t[-1]])
            tall = np.sort(np.unique(np.concatenate((self._t0, t))))
            self._t0 = tall


    def read_gauges_list(self, run_id_list):
        r"""
        read gauge response for a list of run_ids

        """
        grb_list = []

        for run_id in run_id_list:
            grb_list.append(self.read_gauges(run_id))

        return grb_list


    def read_gauges(self,run_id):
        r"""
        reads and returns gauge responses, organized in a dictionary 
        w/ gauge name as keys. each dict item (corresp. to a key) has a list of 
        gauge responses, as given by hash_id_list
    
        ready for multiple gauges...
    
        INPUT hash_id_list
                a list of integers, e.g. [0, 1, 2, 3]
    
        OUTPUT gauge_response
                a Gauge_response class
    
        TODO: pick which gauge
    
        """
    
        grb = Gauge_response_bundle()
        grb.response_list = []

        path_run = 'run_' + str(run_id)

        gauge_names_list = [gauge_info[0] for gauge_info in \
                                self.GeoClawInput._rundata.gaugedata.gauges]
        
        for j,gauge_name in enumerate(gauge_names_list):
            gauge_name_str = '{:05}'.format(gauge_name)
            gauge_fname = os.path.join(path_run,'_output','gauge'+ gauge_name_str + '.txt')
            gauge_data_array = np.loadtxt(gauge_fname,skiprows=2)

            gauge_response = Gauge_response()
            gauge_response.x = self.GeoClawInput._rundata.gaugedata.gauges[j][1]
            gauge_response.y = self.GeoClawInput._rundata.gaugedata.gauges[j][2]
            gauge_response.t = gauge_data_array[:,1]
            gauge_response.h = gauge_data_array[:,2]
            gauge_response.hu = gauge_data_array[:,3]
            gauge_response.hv = gauge_data_array[:,4]
            gauge_response.eta = gauge_data_array[:,5]
            gauge_response.run_id = run_id
            gauge_response.gauge_id = gauge_name
            grb.response_list.append(gauge_response)

        return grb
    

    def bubble_interp(self):
        r"""
        bubble interpolation
        
        """

    def manual_interp(self,m_list):
        r"""
        build piecewise linear model using points in m_list 
        
        """

        # save input list as path_common/input_list.txt
        # if old file exists, try to continue the old run
        m_array = np.array(m_list)
        input_list_fname = os.path.join(self.path_common, 'input_list.txt')
        
        if os.path.exists(input_list_fname):
            print('(manual_interp) already exists: ' + input_list_fname)
            m_array_old = np.loadtxt(input_list_fname)
            if (m_array == m_array_old).all():
                print('(manual_interp) old input list is identical.')
                print('(manual_interp) continuing old run...')
                self._continue = True
            else:
                print('(manual_interp) new input is different: ' \
                        + 'reconcile different inputs first.')
                raise StandardError 
        else: 
            if not os.path.exists(self.path_common):
                print('(manual_interp) creating path: ' + self.path_common)
                os.mkdir(self.path_common)
            np.savetxt(input_list_fname, m_array)
        
        
        for m in m_list:
            self.evaluate_hdm(m)


    def reversal1d(self, grb_list):

        n = len(grb_list)

        for j in range(n):
            m = len(grb_list[j].response_list)
            for k in range(m):
                self._update_t0(grb_list[j].response_list[k].t)
        
        # set pivot
        gr0 = self.interp2t0(grb_list[0].response_list[0])

        cb_list = [None]*n

        for j in range(n):
            m = len(grb_list[j].response_list)
            cb_list[j] = [None]*m
            for k in range(m):
                gr1 = self.interp2t0(grb_list[j].response_list[k])
                p = len(gr1.t)
                nu = 0
                min_l2diff = np.linalg.norm(gr0.eta - gr1.eta)
                for l in range(p):
                    ii = range(p)
                    ii = ii[l:] + ii[:l] 
                    l2diff = np.linalg.norm(gr0.eta - gr1.eta[ii])


                    if l2diff < min_l2diff:
                        min_l2diff = l2diff
                        nu = l

                cb_list[j][k] = nu


        return cb_list


    def read_eta_max(self,run_id):
        r"""
        read_eta_max

        """
        import pickle

        run_id_dir = 'run_' + str(run_id)
        pkl_fname = os.path.join(run_id_dir, 'eta_max_bundle.pkl')

        if os.path.exists(pkl_fname):

            #print('\t(read_eta_max) the pickle exists. ')
            #print('\t(read_eta_max) reading...' + pkl_fname)

            with open(pkl_fname, mode='r') as infile:
                emb = pickle.load(infile)

        else:
            print('\t(read_eta_max) no pickle, processing gauges... ')
            grb = self.read_gauges(run_id)
            emb = Eta_max_bundle()

            emb.response_list = []

            for i, gr in enumerate(grb.response_list):
                em = Eta_max()

                em.run_id = gr.run_id
                em.x = gr.x
                em.y = gr.y

                imax = np.argmax(gr.eta)
                em.t = gr.t[imax]
                em.eta_max = gr.eta[imax] - gr.eta[0]

                emb.response_list.append(em)

            with open(pkl_fname, mode='w') as outfile:
                print('\t(read_eta_max) saving to file: ' + pkl_fname)
                pickle.dump(emb,outfile)

        return emb



    def compute_total_eta_max(self):
        r"""
        compute the sum of eta_max over all runs over the bundles

        """

        total_eta_max_list = []
        n = len(self.input_list)

        for j in range(n):
            emb = self.read_eta_max(j)
            val = 0.
            for k,em in enumerate(emb.response_list):
                val += em.eta_max

            total_eta_max_list.append(val)

        return total_eta_max_list



    def save_eta_max_txt(self):

        n = len(self.input_list)

        for j in range(n):
            emb = self.read_eta_max(j)

            rundir = 'run_' + str(j)
            xgrid_fname = os.path.join(rundir,'x_grid.txt')
            ygrid_fname = os.path.join(rundir,'y_grid.txt')
            eta_max_fname = os.path.join(rundir,'eta_max_bundle'+str(j) +'.txt')

            x = []
            y = []
            eta_max = []

            for k,em in enumerate(emb.response_list):
                x.append(em.x)
                y.append(em.y)
                eta_max.append(em.eta_max)
                
            x = np.array(x)
            y = np.array(y)
            eta_max = np.array(eta_max)
            np.savetxt(xgrid_fname,x)
            np.savetxt(ygrid_fname,y)
            np.savetxt(eta_max_fname,eta_max)



    def read_eta_max_list(self,run_id_list):

        emb_list = []

        for run_id in run_id_list:
            emb_list.append(self.read_eta_max(run_id))

        return emb_list



    def save_drom(self, fname='drom.pkl'):
        r"""
        pickle object and save to common/drom.pkl

        apparently, simple pickling does not work.

        """

        import pickle

        fullfname = os.path.join(self.path_common,fname)
        self._pickle_fname = fname

        if os.path.exists(fullfname):
            print('\t(save_drom) pickle exists: will be backed up.')
            backupname = fullfname + '.backup'
            os.system('mv ' + fullfname + ' ' + backupname)

        with open(fullfname, mode='w') as outfile:
            print('\t(save_drom) saving to file: ' + fullfname)
            pickle.dump(self,outfile)



    def cleanup_pickle(self):
        r"""
        remove pickle

        """
        fullfname = os.path.join(self.path_common,self._pickle_fname)
        if os.path.exists(fullfname):
            print('\t(cleanup_pickle) removing file: ' + fullfname)
            os.system('rm ' + fname)
        else:
            print('\t(cleanup_pickle) missing file.')


    def cleanup_runs(self):
        r"""

        remove all subdirectorys with name run_#

        """

        for j in self._run_id_list:
            runpath = 'run_' + str(j)
            if os.path.exists(runpath):
                print('removing dir: ' + runpath)
                os.system('rm -rf ' + runpath)
            else:
                print('dir: ' + runpath + ' : not found.')

        self.GeoClawInput = GeoClawInput()
        self.GeoClawInput.path_common = self.path_common # this looks bad
        self.input_list = []
        self._delaunay = None
        self._run_id_list = []
        self._run_id2node_id = None
        self._node_id2run_id = None
        self._run_id = 0
        self._t0 = None


# ==============================================================================
#  GeoClawInput Class
# ==============================================================================

class GeoClawInput(object):
    r"""
    GeoClawInput object containing input info for a geoclaw run

    TODO (2017/03/12)

        * make class iterable
        * user supplies the class GeoClawInput_iter
          or defaults to one of the built-in iterables
        - set_input() should be given a more specific name, like 
          set_KL_slips()
        * set up a built-in iterable which receives a list of KL-inputs
          -> simplify the CC_CSZ_South run.
        * handle different executables / different common files
        * write "clean_up" function that will remove all directories 
          will iterate over all GeoClawInput and delete directories
    
    """


    #@property
    #def input_dim(self):
    #    if not self.KL:
    #        self._slip_dim = len(self.fault.subfaults)
    #    return self._slip_dim

    # TODO: remove of this
    @property
    def input_lims(self):
        if not self.KL:
            self._input_lims = [0.,100.]    # beyond 100m slip is unlikely
        if self.KL:
            self._input_lims = [-4.,4.]     # covers > 99.80 % of std normal
        return self._input_lims

    @property
    def KL_Mw_desired(self):
        return self._KL_Mw_desired

    @property
    def KL_Mw_desired(self):
        return self._KL_Mw_desired

    @KL_Mw_desired.setter
    def KL_Mw_desired(self, value):
        if value != self._KL_Mw_desired:
            self._KL_recalculate = True
        self._KL_Mw_desired = value

    
    def __init__(self, iter_fun=None, run_home=None):
        r"""
        GeoClawInput class initialization

            * class containing input values to geoclaw

        """

        self.path_common = '.'
        self._KL_input = None               # current input
        self._slip_dim = None               # slip dimension (varies dep on KL)
        self._input_lims = None             # limits
        
        self.fault = dtopotools.Fault()     # dtopotools.Fault()
        self.dtopo = None                   # dtopotools.DTopography()
        
        self.KL = False                     # use KL expansion
        self._KL_Mw_desired = 8.0           # desired Mw
        self._KL_taper = False              # use taper
        self._KL_nterms = None              # number of terms in KL expansion
        self._KL_tau = None                 # supply taper
        self._KL_distribution = None        # distribution for KL
        self._KL_Lstrike = None             # corr length strike direction
        self._KL_Ldip = None                # corr length dip direction
        self._KL_H = None                   # Hurst exponent
        self._KL_ACF = None                 # Autocorrelation
        self._KL_alpha = None               # alpha parameter
        self._KL_mean_slip = None           # mean slip 
        self._KL_modes = None               # KL evals and evectors
        self._KL_Mo_desired = None          # desired Mw
        
        self._KL_recalculate = False        # recompute KL expansion
        #self.input_type = input_type       # fault parameter is main example
        
        self._rundata = None
        self._claw_pkg = None
        self._num_dim = None                # number of spatial dimensions
        self._dtopo_fname = None

        self._run_id = 0                    # given iteration id
        self._rundir = None                 # designated rundir
        self._next = iter_fun               # iterator 
        self._input_info = None

        # set run directory
        if run_home == None:
            self._run_home = os.getcwd()
        else:
            self._run_home = run_home


    def __iter__(self):
        return self

    def __next__(self):
        r"""
        
            users use user-supplied iterator
        
        """
        if self._next == None:
            raise StopIteration()
        return self._next(self)

    def set_iter(self, iter_fun):
        r"""
        
            set iteration function
        
        """
        self._next = iter_fun

    def set_rundata(self, setrun=None, setgeo=None):
        r"""
        sets data.ClawRunData and saves it in self._rundata
        
            * user supplies setrun() and setgeo() routines as kwargs
             (when setrun is supplied, setgeo is also assumed to be supplied)
            * if setrun is left as None, it will use a default setrun
        
        """
        if (setrun == None):
            # if setrun is not supplied, setup an *empty* rundata
            setrun = self._setrun
            self._rundata = setrun()
        elif ((setrun != None) and (setgeo != None)):
            self._rundata = setrun(setgeo)
        else:
            # TODO raise error here
            print('setrun and setgeo both must be supplied!')


    def set_KL_slip(self, m):
        r"""
        
        sets input parameter m. 
        dimension of array m should match self._slip_dim
        
        if self.KL = True then compute the slips in KL basis then 
        set the slips for self.fault

        if self.KL = False then directly set slips for self.fault 

        the ordering of m: for KL its obvious, for non-KL it is 
        the list ordering of fault.subfaults

        """

        self._KL_input = m     # save current input
        num_slips = len(self.fault.subfaults)
        
        # compute slip for each input
        if self.KL:
            slip_array = self._KL_mean_slip.copy()
            eigenvals,V = self._KL_modes
            for k in range(self._KL_nterms):
                slip_array += m[k] * np.sqrt(eigenvals[k]) * V[:,k]
        else:
            slip_array = np.array(m)

        if self._KL_distribution == 'Lognormal':
            slip_array = np.exp(slip_array)

        # set slips for each subfault
        if self._KL_taper:
            for j,s in enumerate(self.fault.subfaults):
                s.slip = slip_array[j] * self._KL_tau(s.depth)

            # rescale to have desired magnitude
            Mo = self.fault.Mo()
            slip_array *= self._KL_Mo_desired / Mo

            for j,s in enumerate(self.fault.subfaults):
                s.slip = slip_array[j] * self._KL_tau(s.depth)

    def _KL_reexpand(self):
        r"""
            re-expand the KL if flagged
            re-expansion is done if
            
                * self.KL is true (KL is to be used at all)
                * self._KL_recalculate (some parameters have been changed)
            
            the input parameters for the expansion are used as they
            are stored in the variables:
            
                * nterms = self._KL_nterms
                * Lstrike = self._KL_Lstrike
                * Ldip = self._KL_Ldip
                * H = self._KL_H
                * ACF = self._KL_ACF
                * distribution = self._KL_distribution
                * alpha = self._KL_alpha
                * tau = self._KL_tau
                * KL_Mw_desired = self._KL_Mw_desired
        
        """
        if (self.KL and self._KL_recalculate):
            self.KL_expand(\
                nterms = self._KL_nterms, \
                Lstrike = self._KL_Lstrike, \
                Ldip = self._KL_Ldip, \
                H = self._KL_H, \
                ACF = self._KL_ACF, \
                distribution = self._KL_distribution, \
                alpha = self._KL_alpha, \
                tau = self._KL_tau, \
                KL_Mw_desired = self._KL_Mw_desired)
            self.set_KL_slip(self._KL_input)
            self._KL_recalculate = False


    def write_files(self, rundir=None):
        """
        Write all GeoClaw input files, fgmax/fgout data, and generate dtopo + plots.
        Ensures compatibility with evaluate_hdm.
        """

        import pylab as pl
        import shutil

        if rundir is None:
            rundir = self._rundir
        if not os.path.exists(rundir):
            os.makedirs(rundir)
            print(f"[DEBUG] Created run directory {rundir}")

        # --- 1. Write standard .data files ---
        self._rundata.write(out_dir=rundir)
        print(f"[DEBUG] Wrote standard .data files into {rundir}")

        # --- 2. Verify required .data files ---
        for fname in [
            "claw.data", "amr.data", "geoclaw.data", "refinement.data",
            "regions.data", "gauges.data", "friction.data", "qinit.data",
            "multilayer.data", "topo.data", "surge.data"
        ]:
            fpath = os.path.join(rundir, fname)
            if os.path.exists(fpath):
                print(f"[DEBUG] Verified {fname}")
            else:
                print(f"[WARNING] {fname} not found in {rundir}")

        # --- 3. FGout and FGmax info ---
        if hasattr(self._rundata, "fgout_data"):
            print(f"[DEBUG] FGout grids in rundata: {self._rundata.fgout_data.fgout_grids}")
        if hasattr(self._rundata, "fgmax_data"):
            print(f"[DEBUG] FGmax grids in rundata: {self._rundata.fgmax_data.fgmax_grids}")

        # --- 4. Ensure xgeoclaw binary is available ---
        exe = os.path.join(rundir, "xgeoclaw")
        driver_home = os.getcwd()
        if not os.path.exists(exe):
            src_exe = os.path.join(driver_home, "xgeoclaw")
            if os.path.exists(src_exe):
                shutil.copy(src_exe, exe)
                print(f"[DEBUG] Copied xgeoclaw → {exe}")
            else:
                print(f"[WARNING] xgeoclaw not found in {driver_home}")
        else:
            print(f"[DEBUG] xgeoclaw already present in {rundir}")

        # --- 5. Generate dtopo and plots ---
        dtopo_fname = os.path.join(rundir, "dtopo.tt3")
        dtopo_plot_fname = os.path.join(rundir, "dtopo_plot.png")
        subfaults_plot_fname = os.path.join(rundir, "subfaults_plot.png")

        # Create dtopography
        print("[DEBUG] Creating dtopo...")
        x, y = self.fault.create_dtopo_xy()
        self.dtopo = self.fault.create_dtopography(x, y, times=[1.0])
        self.dtopo.write(dtopo_fname, dtopo_type=3)
        print(f"[DEBUG] Wrote dtopo file {dtopo_fname}")

        # save plot of subfaults
        fig0,ax0 = pl.subplots()
        self.fault.plot_subfaults(axes=ax0, slip_color=True)
        fig0.savefig(subfaults_plot_fname,bbox_inches='tight')
        pl.close(fig0)

        # save plot of dtopo
        fig1,ax1 = pl.subplots()
        self.dtopo.plot_dZ_colors(0.,axes=ax1)
        fig1.savefig(dtopo_plot_fname,bbox_inches='tight')
        pl.close(fig1)

        # Recalculate KL expansion if needed
        if hasattr(self, "_KL_reexpand"):
            print("[DEBUG] Recomputing KL expansion...")
            self._KL_reexpand()

        # --- 6. Write dtopo.data properly for GeoClaw ---
        # Configure BEFORE writing to ensure dt_max_dtopo is saved
        dtopo_type = 3  # Okada model
        self._rundata.dtopo_data.dtopofiles = [[dtopo_type, 'dtopo.tt3']]
        # self._rundata.dtopo_data.dtopofiles.append([dtopo_type, 'dtopo.tt3'])        
        # self._rundata.dtopo_data.dt_max_dtopo = 0.2       # in setrun we set 1e+99

        print(f"[DEBUG] Writing GeoClaw input files → {rundir}")
        self._rundata.write(out_dir=rundir)
        print(f"[DEBUG] Updated dtopo.data written with dt_max_dtopo = 0.2")

        # --- 7. Ensure all .data files exist in _output for xgeoclaw ---
        output_dir = os.path.join(rundir, "_output")
        os.makedirs(output_dir, exist_ok=True)

        # List of data files GeoClaw expects
        required_data_files = [
            "claw.data", "amr.data", "geoclaw.data", "refinement.data",
            "regions.data", "gauges.data", "friction.data", "qinit.data",
            "multilayer.data", "topo.data", "surge.data", "fgmax_grids.data",
            "fgout_grids.data", "dtopo.data"
        ]

        for fname in required_data_files:
            src = os.path.join(rundir, fname)
            dst = os.path.join(output_dir, fname)
            if os.path.exists(src):
                shutil.copy(src, dst)
            elif os.path.exists(dst):
                # already there, fine
                continue
            else:
                print(f"[WARNING] Missing {fname} in {rundir} and _output")

        print(f"[DEBUG] Copied all .data files into {output_dir}")

       
          
    def compute_subfault_distances(self,fault):
        """
        Estimate the distance between subfaults i and j for every pair in the list
        fault.subfaults.
    
        :Inputs:
          -  *fault* of class dtopotools.Fault or some subclass,
        
        :Outputs:
          - *D* array of Euclidean distances based on longitudes, latitudes, and depths
          - *Dstrike* array of estimated distances along strike direction
          - *Ddip* array of estimated distances along dip direction
        with D**2 = Dstrike**2 + Ddip**2 to within roundoff.
    
        For each array, the [i,j] entry is distance from subfault i to j when
        ordered in the order the subfaults appear in the list fault.subfaults.
    
        Distance in dip direction based on differences in depth.  
    
        Code by R.J. LeVeque
        
        """
    
        import numpy
        from numpy import pi,sqrt,sin,cos,tan
    
        rad = pi/180.       # conversion factor from degrees to radians
        rr = 6.378e6        # radius of earth
        lat2meter = rr*rad  # conversion factor from degrees latitude to meters
        
        nsubfaults = len(fault.subfaults)
        D = numpy.zeros((nsubfaults,nsubfaults))
        Dstrike = numpy.zeros((nsubfaults,nsubfaults))
        Ddip = numpy.zeros((nsubfaults,nsubfaults))
        for i,si in enumerate(fault.subfaults):
            xi = si.longitude
            yi = si.latitude
            zi = si.depth
            for j,sj in enumerate(fault.subfaults):
                xj = sj.longitude
                yj = sj.latitude
                zj = sj.depth
                dx = abs(xi-xj)*cos(0.5*(yi+yj)*pi/180.) * lat2meter
                dy = abs(yi-yj) * lat2meter
                dz = abs(zi-zj)
    
                # Euclidean distance:
                D[i,j] = sqrt(dx**2 + dy**2 + dz**2)
                
                # estimate distance down-dip based on depths:
                dip = 0.5*(si.dip + sj.dip)
                ddip1 = dz / sin(dip*pi/180.)
                Ddip[i,j] = ddip1 
                if Ddip[i,j] > D[i,j]:
                    # should not happen...
                    if 0:
                        print ("i,j,dx,dy,dz: ",i,j,dx,dy,dz)
                        print ("*** Ddip = %s, D = %s" % (Ddip[i,j], D[i,j]))
    
                # compute distance in strike direction to sum up properly:
                dstrike2 = max(D[i,j]**2 - Ddip[i,j]**2, 0.)
                Dstrike[i,j] = sqrt(dstrike2)
                    
        return D,Dstrike,Ddip
        
        
    def KL_expand(self, nterms=2, Lstrike=400e3, Ldip=40e3, H=0.3,\
          ACF='Exponential', distribution='Gaussian', alpha=0.5, tau=None,\
          KL_Mw_desired = None):
        r"""
        performs KL expansion using spatial correlation designated by
        Lstrike and Ldip


        set self._KL_nterms = nterms after the expansion.

        H: Hurst exponent

        REQUIRES slip_tools

        TODO: more reasonable defaults for Lstrike and Ldip

            original code by R.J. LeVeque

        TODO: 
            * add in slip_tools (remove slip_tools)
            * add a template function for tau fctn
            
        """

        if tau != None:
            self._KL_taper = True
            self._KL_tau = tau

        self._KL_Lstrike = Lstrike 
        self._KL_Ldip = Ldip 
        self._KL_H = H 
        self._KL_ACF = ACF 
        self._KL_alpha = alpha 

        n = len(self.fault.subfaults)

        D, Dstrike, Ddip = self.compute_subfault_distances(self.fault)

        #print "Correlation lengths: Lstrike = %g, Ldip = %g" % (Lstrike,Ldip)
        #print "Hurst exponent = %4.2f" % H

        # exponential correlation
        r = np.sqrt((Dstrike/Lstrike)**2 + (Ddip/Ldip)**2)
        C = np.exp(-r)

        # TODO: allow various ACF: von Karman, Gaussian
        #C = np.exp(-((Dstrike/Lstrike)**2 + (Ddip/Ldip)**2))
        
        lengths = np.array([s.length for s in self.fault.subfaults])
        widths = np.array([s.width for s in self.fault.subfaults])
        areas = lengths * widths
        total_area = sum(areas)

        self._KL_Mw_desired = KL_Mw_desired
        
        Mo_desired = 10.**(1.5*KL_Mw_desired + 9.05)
        self._KL_Mo_desired = Mo_desired
        mean_slip = Mo_desired / (self.fault.subfaults[0].mu * total_area)
        print('mean_slip ' + str(mean_slip) + ' meters required ' + \
              'for Mw ' + str(KL_Mw_desired))
        
        # Turn this into a constant vector:
        mean_slip = mean_slip * np.ones(n)
        sigma_slip = alpha * mean_slip

        if distribution == 'Gaussian':
            
            Cov = sigma_slip * (C*sigma_slip).T
            self._KL_distribution = 'Gaussian'

        elif distribution == 'Lognormal':

            Cov = np.log((sigma_slip/mean_slip)*(C*(sigma_slip/mean_slip)).T\
                          + 1.)
            mean_slip = np.log(mean_slip) - np.diag(Cov)/2.
            self._KL_distribution = 'Lognormal'

        # Find eigenvalues, and eigenvector matrix.
        # Columns V[:,k] are eigenvectors.

        print ("Finding eigenmodes from %s by %s matrix C" % (n,n))
        eigenvals, V = np.linalg.eig(Cov)

        eigenvals = np.real(eigenvals)  # imag parts should be at rounding level
        V = np.real(V)

        # Sort eigenvalues:
        i = list(np.argsort(eigenvals))
        i.reverse()
        eigenvals = eigenvals[i]
        V = V[:,i]

        lam = eigenvals

        if not self.KL:
            self.KL = True
        
        self._KL_nterms = nterms
        self._KL_mean_slip = mean_slip.copy()
        self._KL_modes = (eigenvals,V)


    def KL_save(self, fname):
        r"""
        saves KL expansion in numpy.array binary format to a file 

        """

    def _setrun(self,claw_pkg='geoclaw',num_dim=2):
        r"""
        return an *empty* data.ClawRunData object
        
        """

        from clawpack.clawutil import data

        self._claw_pkg = claw_pkg
        self._num_dim = num_dim
        
        rundata = data.ClawRunData(claw_pkg, num_dim)
        rundata = self._setgeo(rundata)

        return rundata


    def _setgeo(self,rundata):
        r"""
        default setgeo routine
        """
        
        return rundata


# ==============================================================================
#  Gauge_response_bundle Class
# ==============================================================================

class Gauge_response_bundle(object):
    r"""
    
    bundle of gauge responses

    """

    @property
    def xy_list(self):
        if self._xy_list is None:
            self._get_bundle_info()
        return self._xy_list
    
    @property
    def gauge_id_list(self):
        if self._gauge_id_list is None:
            self._get_bundle_info()
        return self._gauge_id_list
    
    @property
    def run_id(self):
        if self._run_id is None:
            self._get_bundle_info()
        return self._run_id

    def __init__(self):

        #super(gauge_response_bundle, self).__init__()
        
        self.response_list = None
        self._xy_list = None
        self._gauge_id_list = None
        self._run_id = None
        
    def _get_bundle_info(self):
        r"""
        loop through gauge responses and gather information
        
        set the variables
        
        self.response_list 
        self.xy_list 
        self.gauge_id_list 
        self.run_id 
        
        """
        N = len(self.response_list)
        if  N > 0:
            self._run_id = self.response_list[0].run_id
            self._gauge_id_list = []
            self._xy_list = []
            for j in range(N):
                xy = (self.response_list[j].x, self.response_list[j].y)
                self._xy_list.append(xy)
                self._gauge_id_list.append(self.response_list[j].gauge_id)



    
# ==============================================================================
#  Gauge_response Class
# ==============================================================================

class Gauge_response(object):
    r"""
    a class containing a gauge output, from a single gauge

    TODO: the output variables should be allowed to be arbitrary.
    """


    def __init__(self):

        self.x = None
        self.y = None
        self.t1 = None
        self.t2 = None
        self.t = None
        self.h = None
        self.hu = None
        self.hv = None
        self.eta = None
        self.run_id = None
        self.gauge_id = None


# ==============================================================================
#  Eta_max_bundle Class
# ==============================================================================

class Eta_max_bundle(object):
    r"""
    
    bundle of gauge responses

    """

    @property
    def xy_list(self):
        if self._xy_list is None:
            self._get_bundle_info()
        return self._xy_list
    
    @property
    def gauge_id_list(self):
        if self._gauge_id_list is None:
            self._get_bundle_info()
        return self._gauge_id_list
    
    @property
    def run_id(self):
        if self._run_id is None:
            self._get_bundle_info()
        return self._run_id

    def __init__(self):

        #super(gauge_response_bundle, self).__init__()
        
        self.response_list = None
        self._xy_list = None
        self._run_id = None
        
    def _get_bundle_info(self):
        r"""
        loop through gauge responses and gather information
        
        set the variables
        
        self.response_list 
        self.xy_list 
        self.gauge_id_list 
        self.run_id 
        
        """
        N = len(self.response_list)
        if  N > 0:
            self._run_id = self.response_list[0].run_id
            self._gauge_id_list = []
            self._xy_list = []
            for j in range(N):
                xy = (self.response_list[j].x, self.response_list[j].y)
                self._xy_list.append(xy)
                self._gauge_id_list.append(self.response_list[j].gauge_id)
    



# ==============================================================================
#  Eta_max Class
# ==============================================================================

class Eta_max(object):
    r"""
    a class containing eta_max (maximum elevation), from a single gauge

    """

    def __init__(self):

        self.x = None
        self.y = None
        self.t = None
        self.eta_max = None
        self.run_id = None
        self.gauge_id = None


