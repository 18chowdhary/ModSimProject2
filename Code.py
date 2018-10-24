from modsim import *


init = State(Vin=0, Vm=0, Vout=0)
system = System(C1=10e-7, C2=10e-7, R1=1580, R2=1580, f1=0, f2=4, A=0.5, init=init, t0=0, freqs = 10, stepresf = 5, numwavels = 3)


def slope_func(init, t, system):
    unpack(system)
    vin, vm, vout = init

    dvin = 2 * np.pi * A * f * np.cos(2*np.pi*f*t)
    dvm = dvin - (vm/R1 + vout/R2) / C1
    dvout = dvm - vout / (R2*C2)

    return dvin, dvm, dvout

def run_bode():
    unpack(system)

    farray = np.logspace(f1, f2, freqs)

    Re = TimeSeries()

    for f in farray:
        system.set(f=f, t_end = numwavels / f)
        max_step = (system.t_end - t0) / (stepresf * f)
        results, details = run_ode_solver(system, slope_func, max_step = max_step)
        amplitudeM = results.Vout.max() - results.Vout.min()
        Re[f] = amplitudeM

    return Re

def run_calc():
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


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
lns1 = ax.plot(run_bode(), label = 'simulated')
lns2 = ax.plot(run_calc(), label = 'calculation')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'best')

savefig('graphs\BodePlot.png')
