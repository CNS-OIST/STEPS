
import os.path as path
from . import rallpack1_membres as rallpack1

def setup_module():
    global C

    # defaults
    C={ 'meshdir': 'validation_efield/meshes',
        'mesh': 'axon_cube_L1000um_D866nm_1978tets',
        'datadir': 'validation_efield/data/rallpack1_correct',
        'v0data': 'v0',
        'v1data': 'vx',
        'seed': 7 }



def test_rallpack1():
    params = rallpack1.sim_parameters

    meshfile = path.join(C['meshdir'],C['mesh'])
    v0data = path.join(C['datadir'],C['v0data'])
    v1data = path.join(C['datadir'],C['v1data'])
    seed = C['seed']

    simdata, rms_err_0um, rms_err_1000um = rallpack1.run_comparison(seed, meshfile, v0data, v1data)

    print("MEMBRES: rms error at 0um = " + str(rms_err_0um*1e3)+ " mV")
    print("MEMBRES: rms error at 1000um = " + str(rms_err_1000um*1e3) + "mV")
    
    import matplotlib.pyplot as plt
    plt.subplot(211)
    plt.plot(simdata[0,:], simdata[2,:], 'k-' ,label = 'Correct, 0um', linewidth=3)
    plt.plot(simdata[0,:], simdata[1,:], 'r--', label = 'STEPS, 0um', linewidth=3)
    plt.legend(loc='best')
    plt.ylabel('Potential (mV)')
    plt.subplot(212)
    plt.plot(simdata[0,:], simdata[4,:], 'k-' ,label = 'Correct, 1000um', linewidth=3)
    plt.plot(simdata[0,:], simdata[3,:], 'r--', label = 'STEPS, 1000um', linewidth=3)
    plt.legend(loc='best')
    plt.ylabel('Potential (mV)')
    plt.xlabel('Time (ms)')
    plt.show(block=True)

    max_rms_err = 1.e-3
    assert(rms_err_0um < max_rms_err)
    assert(rms_err_1000um < max_rms_err)

