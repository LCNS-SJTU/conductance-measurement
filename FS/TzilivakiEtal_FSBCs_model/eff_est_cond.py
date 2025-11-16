from neuron import h
import matplotlib.colors as mcolors
from neuron.units import ms, mV
from matplotlib import pyplot, cm
import numpy as np
import matplotlib.pyplot as plt
from neuron import clear_gui_callback
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from numpy.linalg import pinv, eig
import numpy as np
from multiprocessing import Pool
font2={'family':'Times New Roman',
'weight':'bold',
'size': 10}
# 自定义刻度标签显示格式


h.load_file('nrngui.hoc')
h.load_file("import3d.hoc")

morphologyFilename = "/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/FS/TzilivakiEtal_FSBCs_model/Multicompartmental_Biophysical_models/experiment/model_reconstructions/PFC/Mar11.hoc"
biophysicalModelTemplateFilename = "/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/FS/TzilivakiEtal_FSBCs_model/Multicompartmental_Biophysical_models/experiment/PFCtemplate.hoc"

h.load_file(biophysicalModelTemplateFilename)

def error_eff_cond_est_two(run_dt, fE, fI, C=0.00015, GL = 1e-05, locE=58, locI=37,
                         plot =False):
    
    neuron = h.FScell(morphologyFilename)
        
    h.dt = run_dt
    h.tstop = 600.

    h.v_init = -65.46
    Vclamp_init = -65.46
    recv = h.Vector().record(neuron.soma[0](0.5)._ref_v)
    
    # Initialize an excitatory and inhibitory synaptic input
    # make a new stimulator
    stimE = h.NetStim()
    stimI = h.NetStim()

    synE = h.Exp2Syn(neuron.dend[int(locE)](0.5))
    synI = h.Exp2Syn(neuron.dend[int(locI)](0.5))
    synE.tau1 = 1.  # ms
    synE.tau2 = 5.  # ms
    synI.tau1 = 1.  # ms
    synI.tau2 = 5.  # ms

    synE.e = 4.54  # mV for E input and -80 for I input 
    synI.e = -75.46
        
    stimE.number = 1
    stimE.start = h.tstop - 100.  # ms
    stimI.number = 1
    stimI.start = h.tstop - 100.  # ms
    ncstimE = h.NetCon(stimE, synE)
    ncstimE.delay = 0.
    ncstimE.weight[0] = 0e-2  # uS
    ncstimI = h.NetCon(stimI, synI)
    ncstimI.delay = 0.
    ncstimI.weight[0] = 0e-2  # uS


    # Initialize voltage clamp
    vclamp = h.SEClamp(neuron.soma[0](0.5))
    vclamp.amp1 = Vclamp_init
    vclamp.dur1 = h.tstop
    vclamp.dur2 = h.tstop
    vclamp.amp2 = Vclamp_init - 20.

    irec = h.Vector()
    irec.record(vclamp._ref_i)

    # Function to initialize the simulation
    def initialize():
        h.t = 0
        h.finitialize(h.v_init)
        h.fcurrent()

    # Function to integrate the simulation
    def integrate():
        while h.t < h.tstop:
            h.fadvance()

    # Function to run the simulation
    def go():
        initialize()
        integrate()

    # Main simulation loop
    def main():

        Vclamp = np.arange(-110.,-60.,10)
        tstart = int(stimI.start/run_dt+1)
        
        #  pure E for the estimation of true GE_est
        vclamp.dur1 = 0.
        vclamp.dur2 = 0.
        vclamp.dur3 = 0.
        ncstimE.weight[0] = fE
        ncstimI.weight[0] = 0.
        go()
        recvE0 = np.array(recv.to_python())
        V_rest = recvE0[tstart]
        recvE0_norm = recvE0[tstart:tstart + 1000] - V_rest
        epsilonE = 4.54 - V_rest
        epsilonI = -75.46 - V_rest
        epsilonI1 = -95.46 - V_rest
        dVdtE = np.diff(recvE0_norm)/run_dt
        dVdtE = np.append(dVdtE, 0)
        GE_est = (C*dVdtE + GL*recvE0_norm)/(epsilonE - recvE0_norm)
        
        #  pure I input for the estimation of GI_est
        vclamp.dur1 = 0.
        vclamp.dur2 = 0.
        vclamp.dur3 = 0.
        ncstimE.weight[0] = 0.
        ncstimI.weight[0] = fI
        go()
        recvI0 = np.array(recv.to_python())
        recvI0_norm = recvI0[tstart:tstart+1000] - V_rest
        dVdtI = np.diff(recvI0_norm)/run_dt
        dVdtI = np.append(dVdtI, 0)
        GI_est = (C*dVdtI + GL*recvI0_norm)/(epsilonI - recvI0_norm)
        if plot:
            plt.plot(recvE0_norm, 'r')
            plt.plot(recvI0_norm, 'b')
        
        
    
        Isyn_inj_t = np.zeros((len(Vclamp),1000))
        Isyn_epsilonI1_inj_t = np.zeros((len(Vclamp),1000))
        relative_error = np.zeros(4) 

        

        for i in range(len(Vclamp)):
            ncstimE.weight[0] = 0.
            ncstimI.weight[0] = 0.
            vclamp.amp1 = Vclamp[i]
            vclamp.dur1 = h.tstop
            h.v_init = V_rest
            go()
            Iinj0 = np.array(irec.to_python())  # Assuming irec is a NEURON Vector
            
            synE.e = 4.54
            synI.e = -75.46
            ncstimE.weight[0] = fE
            ncstimI.weight[0] = fI
            vclamp.amp1 = Vclamp[i]
            vclamp.dur1 = h.tstop
            h.v_init = V_rest
            go()
            Iinj1 = np.array(irec.to_python())  # Assuming irec is a NEURON Vector
            Isyn = (Iinj1[tstart:tstart+1000] - Iinj0[tstart:tstart+1000])
            Isyn_inj_t[i,:] = Isyn
            
            synE.e = 4.54
            synI.e = -95.46
            ncstimE.weight[0] = fE
            ncstimI.weight[0] = fI
            vclamp.amp1 = Vclamp[i]
            vclamp.dur1 = h.tstop
            h.v_init = V_rest
            go()
            Iinj1_epsilonI1 = np.array(irec.to_python())  # Assuming irec is a NEURON Vector
            Isyn_epsilonI1 =  (Iinj1_epsilonI1[tstart:tstart+1000] - Iinj0[tstart:tstart+1000])
            Isyn_epsilonI1_inj_t[i,:] = Isyn_epsilonI1
            
                                                         
        traditional_effective_conductance = np.zeros((2,1000))
        intercept_effective_conductance = np.zeros((2,1000))
        A = np.ones((len(Vclamp),2))
        A[:,0] = Vclamp - V_rest
        epsilonE = 4.54 - V_rest
        epsilonI = -75.46 - V_rest
        epsilonI1 = -95.46 - V_rest
        for i in range(1000):
            x, residuals, rank, singular_values = np.linalg.lstsq(A, Isyn_inj_t[:,i], 
                                                                  rcond=None)
            x_epsilonI1, residuals, rank, singular_values = np.linalg.lstsq(A, 
                                Isyn_epsilonI1_inj_t[:,i], rcond=None)
            
            traditional_effective_conductance[1, i] = (x[0]*epsilonE + x[1])/(epsilonE - epsilonI)
            traditional_effective_conductance[0, i] = (x[0]*epsilonI + x[1])/(epsilonI - epsilonE)
           
            intercept_effective_conductance[1,i] = (x[1] - x_epsilonI1[1])/(epsilonI1 - epsilonI)
            intercept_effective_conductance[0,i] = (-intercept_effective_conductance[1,i]*epsilonI - x[1])/epsilonE
        
        relative_error[0] = np.linalg.norm(GE_est*1e6-intercept_effective_conductance[0,:]*1e3)/np.linalg.norm(GE_est*1e6)
        relative_error[1] = np.linalg.norm(GE_est*1e6-traditional_effective_conductance[0,:]*1e3)/np.linalg.norm(GE_est*1e6)
        relative_error[2] = np.linalg.norm(GI_est*1e6-intercept_effective_conductance[1,:]*1e3)/np.linalg.norm(GI_est*1e6)
        relative_error[3] = np.linalg.norm(GI_est*1e6-traditional_effective_conductance[1,:]*1e3)/np.linalg.norm(GI_est*1e6)
        
        GE_est_norm = GE_est*1e6
        GI_est_norm = GI_est*1e6
        GE_est_inter = intercept_effective_conductance[0,:]*1e3  
        GE_est_trad = traditional_effective_conductance[0,:]*1e3
        GI_est_inter = intercept_effective_conductance[1,:]*1e3 
        GI_est_trad = traditional_effective_conductance[1,:]*1e3 
        
        return relative_error, GE_est_norm, GI_est_norm, GE_est_inter, GE_est_trad, GI_est_inter, GI_est_trad
    
    
    relative_error, GE_est_norm, GI_est_norm, GE_est_inter, GE_est_trad, GI_est_inter, GI_est_trad= main()
    N = len(relative_error) + len(GE_est_norm) + len(GI_est_norm) + len(GE_est_inter) + len(GE_est_trad) + len(GI_est_inter) + len(GI_est_trad)
    N1 = len(relative_error)
    N2 = N1 + len(GE_est_norm)
    N3 = N2 + len(GI_est_norm)
    N4 = N3 + len(GE_est_inter)
    N5 = N4 + len(GE_est_trad)
    N6 = N5 + len(GI_est_inter)
    N7 = N6 + len(GI_est_trad)
    results = np.zeros(N)
    results[:N1] = relative_error
    results[N1:N2] = GE_est_norm
    results[N2:N3] = GI_est_norm
    results[N3:N4] = GE_est_inter
    results[N4:N5] = GE_est_trad
    results[N5:N6] = GI_est_inter
    results[N6:N7] = GI_est_trad
    return results


# Your existing functions and code
def error_eff_cond_est_two_parallel(args,):
    stimu0, = args
    results = error_eff_cond_est_two(0.1, 7e-5, 7e-5, C=0.0001470913112090827, GL=1.0017936585289593e-05,
            plot = False, locE = stimu0[0], locI = stimu0[1])
     
    return results
loc = np.arange(0,68,1)
# 生成网格
X, Y = np.meshgrid(loc, loc)

# 将网格重塑为9x9的数组
# 首先将X和Y平铺为一维数组，然后将这两个数组堆叠成二维数组，每一行是一个(x,y)对
stimu = np.vstack([X.ravel(), Y.ravel()]).T



# Create a Pool with desired number of processes (adjust as needed)
num_processes = 55
with Pool(num_processes) as pool:
    # Prepare arguments for parallel processing
    args_list = [(stimu0,) for stimu0 in stimu]

    # Apply the function in parallel
    results = pool.map(error_eff_cond_est_two_parallel, args_list)

results_total = np.zeros((len(stimu), 6004))

for i, result in enumerate(results):
    results_total[i, :] = result
np.savetxt("pair_eff_cond_E7e5_I7e5_whole_relative_reversal_potential.txt", results_total)