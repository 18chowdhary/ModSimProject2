from modsim import *

init = State(Vin=0,
             Vm=0,
             Vout=0)

system = System(C1=10e-7, #Capacitor value 1
                C2=10e-7, #Capacitor value 2
                R1=15.8, #Resistor value 1
                R2=15.8, #resistor value 2
                f=100000,
                A=0.5, #Amplitude of Vin wave
                init=init,
                t0=0,
                stepres = 200,
                numwavels = 4)


def slope_func(init, t, system):
    unpack(system)
    vin, vm, vout = init

    dvin = 2 * np.pi * A * f * np.cos(2*np.pi*f*t)
    dvm = dvin - (vm/R1 + vout/R2) / C1
    dvout = dvm - vout / (R2*C2)

    return dvin, dvm, dvout

def run_bode():
    unpack(system)
    system.set(f=f, t_end = numwavels / f)
    max_step = (system.t_end - t0) / (stepres)
    results, details = run_ode_solver(system, slope_func, max_step = max_step)
    tail = int(details.nfev/(2*np.pi*numwavels))
    amplitudeM = results.Vout.tail(tail).ptp()
    print(amplitudeM)
    print(details.nfev)
    print(tail)
    print('maxx')
    print(results.Vout.tail(tail).idxmax())
    print(results.Vout.tail(tail).max())
    print('minx')
    print(results.Vout.tail(tail).idxmin())
    print(results.Vout.tail(tail).min())
    print('tailx')
    print(results.Vout.tail(tail).first_valid_index())
    print(system.t_end-1/f)

    return results

fig = plt.figure()
ax = fig.add_subplot(111)
lns = ax.plot(run_bode().Vout, label = 'simulated')
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'best')

savefig('graphs\SingleFrequency.png')
