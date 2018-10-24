from modsim import *

init = State(Vin=0,
             Vm=0,
             Vout=0)

params = Params(C1=10e-7, #Capacitor value 1
         C2=10e-9, #Capacitor value 2
         R1=300, #Resistor value 1
         R2=15000) #resistor value 2


setsystem = System(f1=1, #10^f1 freq lower bound
                    f2=5, #10^f2 freq upper bound
                    A=0.5, #Amplitude of Vin wave
                    init=init,
                    t0=0,
                    freqs = 100, #number of frequencies to sweep through
                    stepres = 200,
                    numwavels = 4) #num of wavelengths simulated per freq)

def make_system(params):
    system = System(params)
    return system

def slope_func(init, t, system):
    unpack(system)
    vin, vm, vout = init

    dvin = 2 * np.pi * A * f * np.cos(2*np.pi*f*t)
    dvm = dvin - (vm/R1 + vout/R2) / C1
    dvout = dvm - vout / (R2*C2)

    return dvin, dvm, dvout

def run_bode(system):
    unpack(system)

    farray = np.logspace(f1, f2, freqs)
    Re = TimeSeries()

    for f in farray:
        system.set(f=f, t_end = numwavels / f)
        max_step = (system.t_end - t0) / (stepres)
        results, details = run_ode_solver(system, slope_func, max_step = max_step)
        tail = int(details.nfev/(2*np.pi*numwavels))
        amplitudeM = results.Vout.tail(tail).ptp()
        Re[f] = amplitudeM
    return Re

def run_calc(system):
    unpack(system)

    farray = np.logspace(f1, f2, freqs)
    C = TimeSeries()

    for f in farray:
        w = 2*np.pi*f
        rcw1 = R1*C1*w
        rcw2 = R2*C2*w
        amplitudeC = (rcw1*rcw2) / (np.sqrt(1 + rcw1**2) * np.sqrt(1 + rcw2**2))
        C[f] = amplitudeC

    return C

def error_func(params):
    system = make_system(params)

    results = run_bode(system)
    data = run_calc(system)

    errors = results - data
    return errors

error_func(params)

#best_params, fit
"""
system = make_system(params)
results = run_bode()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
lns1 = ax.plot(results, label = 'simulated')
lns2 = ax.plot(run_calc(), label = 'calculation')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'best')
savefig('graphs\BodePlotLowRes.png')
"""
