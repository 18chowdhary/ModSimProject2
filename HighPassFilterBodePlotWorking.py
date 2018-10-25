from modsim import *

init = State(Vin=0,
             Vm=0,
             Vout=0)

params = Params(R1=150, #Resistor value 1
         R2=150) #resistor value 2


setsystem = System(f1=1, #10^f1 freq lower bound
                    f2=5, #10^f2 freq upper bound
                    tau1 = 1.5e-5,
                    tau2 = 1.5e-5,
                    A=0.5, #Amplitude of Vin wave
                    init=init,
                    t0=0,
                    freqs = 10, #number of frequencies to sweep through
                    stepres = 200,
                    numwavels = 4) #num of wavelengths simulated per freq)

def make_system(params, setsystem):
    setsystem.set(params = params)
    system = setsystem
    return system

def slope_func(init, t, system):
    unpack(system)
    R1, R2 = system.params
    C1 = tau1/R1
    C2 = tau2/R2

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
    print('bode run')
    return Re

def run_calc(system):
    unpack(system)
    R1, R2 = system.params
    C1 = tau1/R1
    C2 = tau2/R2

    farray = np.logspace(f1, f2, freqs)
    C = TimeSeries()

    for f in farray:
        w = 2*np.pi*f
        rcw1 = R1*C1*w
        rcw2 = R2*C2*w
        amplitudeC = (rcw1*rcw2) / (np.sqrt(1 + rcw1**2) * np.sqrt(1 + rcw2**2))
        C[f] = amplitudeC

    return C

def error_func(params, setsystem):
    system = make_system(params, setsystem)

    results = run_bode(system)
    data = run_calc(system)

    errors = results - data

    return errors

best_params, fit_details = fit_leastsq(error_func, params, setsystem, maxfev=1)

print(best_params, setsystem.tau1/best_params.R1, setsystem.tau2/best_params.R2)


system = make_system(best_params, setsystem)
results = run_bode(system)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
lns1 = ax.plot(results, label = 'simulated')
lns2 = ax.plot(run_calc(system), label = 'calculation')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'best')
savefig('graphs\BodePlotLowRes2.png')

