from modsim import *

init = State(Vm=0,
             Vout=0)

params = Params(R1=10000, #Resistor value 1
         R2=10000) #resistor value 2


setsystem = System(tau = 1.5e-7,
                   A=0.5, #Amplitude of Vin wave
                   init=init,
                   t0=0,
                   freqs = 50, #number of frequencies to sweep through
                   stepres = 200,
                   numwavels = 4) #num of wavelengths simulated per freq)

def make_system(params, setsystem):
	#@param params: the configuration for the circuit (R1, R2)

	#reates a system object representing the circuit with the correct configuration.
    R1, R2 = params

    setsystem.set(params = params)
    setsystem.set(C1 = setsystem.tau/R1, C2 = setsystem.tau/R2)
    fc = 1/(2*np.pi*setsystem.tau)
    flow = int(np.log10(fc))-2
    fhigh = int(np.log10(fc))+2
    setsystem.set(f1 = flow, f2 = fhigh)
    system = setsystem
    return system

def slope_func(init, t, system):

	#@param init: the initial state of the system
	#@param t: the time at which this function will be evaluated
	#@param system: the system object

	#Returns the changes in Vin, Vm, and Vout
    unpack(system)
    R1, R2 = system.params

    vm, vout = init

    # ODEs for Vin, Vm, Vout
    dvin = 2 * np.pi * A * f * np.cos(2*np.pi*f*t)
    dvm = dvin - (vm/R1 + vout/R2) / C1
    dvout = dvm - vout / (R2*C2)

    return dvm, dvout

def run_bode(system):

	#@param system: the system object

	#Generates a bode plot for a series of frequencies

    unpack(system)

    # Creates a set of frequencies on a log scale
    farray = np.logspace(f1, f2, freqs)

    Re = TimeSeries()

    for f in farray:
        system.set(f=f, t_end = numwavels / f)

        # Determines the maximum step length, to force the bode plot to run faster
        max_step = (system.t_end - t0) / (stepres)
        results, details = run_ode_solver(system, slope_func, max_step = max_step)

        # Uses nfev for the # of steps and to select out the tail
        tail = int(details.nfev/(2*np.pi*numwavels))
        amplitudeM = results.Vout.tail(tail).ptp()
        Re[f] = amplitudeM
    print('bode run')
    return Re

def run_calc(system):
    unpack(system)
    R1, R2 = system.params

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
    print(params)
    system = make_system(params, setsystem)

    results = run_bode(system)
    data = run_calc(system)

    errors = results - data

    return errors

best_params, fit_details = fit_leastsq(error_func, params, setsystem, maxfev=7)

print(best_params, setsystem.tau/best_params.R1, setsystem.tau/best_params.R2)


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
