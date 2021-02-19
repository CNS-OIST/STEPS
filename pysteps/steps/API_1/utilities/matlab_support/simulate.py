from __future__ import division

import sys
import os
import json
import numpy
import datetime
import argparse
import time

if(sys.version_info >= (3,0)):
    import configparser as cfg_parser
else:
    import ConfigParser as cfg_parser

import steps
import steps.API_1.model as smod
import steps.API_1.geom as swm
import steps.API_1.rng as srng
import steps.API_1.solver as ssolver

# get the directory of this file
mod_path = os.path.dirname(__file__)
steps_path = os.path.dirname(steps.__file__) + '/..'

# read the config
config = cfg_parser.ConfigParser()
config.read(mod_path + '/simulate.cfg')
log_dest = config.get('logging', 'destination')
log_level = config.get('logging', 'level')
log_timestamp = config.getboolean('logging', 'use_timestamp') 
print('log level: ' + log_level)

# logging configuration and setup
import logging

# in any other case use screen as fallback
logname = ''
if(log_dest == 'file'):
    if(not os.path.exists(steps_path + '/.logs')):
        os.makedirs(steps_path + '/.logs')

    if(log_timestamp):
       # get timestamp 
       logname = steps_path + '/.logs/STEPS_for_Matlab_' + str(datetime.datetime.now()) + '.log'
    else:
       logname = steps_path + '/.logs/STEPS_for_Matlab.log'

# get the log level
llevel = logging.getLevelName(log_level.upper())

logging.basicConfig(filename=logname,
    format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%d-%m-%Y:%H:%M:%S',
    level=llevel)

logger = logging.getLogger()

# parse a Matlab mathematical expression and convert it to Python
from steps.API_1.utilities.matlab_support.event_parser import ExprParser

def run(model_string, initial_condition_string, stop_time, event_string='', seed=23412, solver='wmdirect', dt=1, rk4dt=1e-5, realisations=1, file=False ):
    """ 
       Run a Matlab Simbiology model in STEPS
       using the Wmdriect or Wmrk4 solver.

       This is the entry point called by Matlab via the main code.

       Input parameters:
       - STEPS model as a string,
       - initial condition part of the STEPS model as a string,
       - stop time
       - pseudo random number generator seed,
       - solver: wmdirect or wmrk4,
       - dt: time increment betwen sampling points
       - realisations: number of times the simulation is to be run
       - file: flag to indicate if the generated model should be 
               written to file. The file name carries a timestamp.
    """
    if(not solver.upper() == 'WMDIRECT' and not solver.upper() == 'WMRK4'):
        logger.error('run unknown solver string: ' + solver)
        raise Exception('run unknown solver string: ' + solver)

    if(file):
        if(not os.path.exists(steps_path + '/.logs')):
            os.makedirs(steps_path + '/.logs')

        outfile = open(steps_path + '/.logs/matlab_steps_model_' + str(datetime.datetime.now()), 'w')
        outfile.write(model_string)
        outfile.write(initial_condition_string)
        outfile.close()

    exec(model_string.encode('ascii').decode('unicode-escape'), globals())

    # setup the solver, only two solvers are supported, 
    # so a simple if construct is sufficient.
    global sim
    if(solver.upper() == 'WMDIRECT'):
        r = srng.create('mt19937', 256)
        r.initialize(seed)
        sim = ssolver.Wmdirect(model, geometry, r)
    elif(solver.upper() == 'WMRK4'):
        sim = ssolver.Wmrk4(model, geometry)
        sim.setRk4DT(rk4dt)
    else:
        logger.error('run unknown solver string: ' + solver)
        raise Exception('run unknown solver string: ' + solver)

    # steup the return data structures
    nSpecies = len(model.getAllSpecs())
    list_of_species = model.getAllSpecs()
    specie_names = []
    for s in list_of_species:
        specie_names.append(s.getID())
    tvector = numpy.arange(0, stop_time, dt)

    # get a list of the compartment names
    compartment_names = []
    compartments = geometry.getAllComps()
    for comp in compartments:
        compartment_names.append(str(comp.getID()))

    res = numpy.zeros([realisations, len(compartments), len(tvector), nSpecies])
    
    # parse the events and
    # store them in a dictionary
    event_string = json.loads(event_string.replace('][', ', '))
    logger.debug('run got events: ' + str(event_string))
    
    # iterate the names of the reaction rates
    rate_names = event_string[1]

    # iterate the events
    event_string = event_string[0]
    
    # number of events, each event is a dict with trigger and event entry
    nevents = len(event_string)    

    event_dict = {}
    # iterate all the dictionaries
    for i in range(nevents):
        # extract the trigger and event
        event_dict.update(event_string[i])

    event_times = {}
    nevents = len(event_dict)
    
    if(nevents % 2 != 0):
        logger.error('run internal error, number of triggers and events do not match: ' + str(event_dict))
        raise RuntimeError('run internal error, number of triggers and events do not match: ' + str(event_dict))
  
    for i in range(nevents // 2):
        trigger = event_dict['trigger_' + str(i)]
        
        if(trigger.find('>=') >= 0):
            tmp = trigger.split('>=')
        elif(trigger.find('>') >= 0):
            tmp = trigger.split('>')
        else:
            logger.error('run invalid event trigger: ' + e)
            raise RuntimeError('run invalid event trigger: ' + e)
        
        #if(event_times.has_key(float(tmp[1]))):
        if(float(tmp[1]) in event_times):
            event_times[float(tmp[1])].append(event_dict['event_' + str(i)])  
        else:
            event_times[float(tmp[1])] = [ event_dict['event_' + str(i)] ]

    # build the event queue - sort the events
    event_queue = sorted(event_times.keys())

    if(not event_queue):
        # run without event queue
        for i in range(realisations): # this is an index
            sim.reset()
            exec(initial_condition_string.encode('ascii').decode('unicode-escape'), globals())
    
            for t in range(len(tvector)): # this is an index
                for c in range(len(compartment_names)): # this is an index
                    for s in range(len(specie_names)):
                        res[i, c, t, s] = sim.getCompCount(compartment_names[c], specie_names[s])
                sim.run(tvector[t])
    else: 
        # run with event queue
        logger.debug('run event times: ' + str(event_times))
        logger.debug('run event queue ' + str(event_queue))
        logger.debug('run event dict: ' + str(event_dict))
        for i in range(realisations):
            # copy the event queue for this realisation
            # FIXME: replace copying the queue with a counter
            eq = list(event_queue)
            t_event = float(eq.pop(0))

            # initial conditions
            sim.reset()
            exec(initial_condition_string.encode('ascii').decode('unicode-escape'), globals())
            
            for t in range(len(tvector)):
                if tvector[t] >= float(t_event):
                    logger.debug('run time of next event: ' + str((float(t_event))))
                    sim.run(float(t_event))
                    
                    # handle all events event queued for this t_event
                    for e in event_times[t_event]:
                        tmp =  handleEvent(e, sim, model, spec_names, comp_names, rate_names, reac_to_comp, sreac_to_patch, kcst_si_factor)
                        exec(tmp)
                    
                    # get next event
                    if eq:
                        t_event = eq.pop(0)
                    else:
                        t_event = float('inf')
                # end if
                
                sim.run(tvector[t])
                
                # get molecule counts
                for c in range(len(compartment_names)): # this is an index
                    #sim.run(tvector[t])
                    for s in range(len(specie_names)):
                        res[i, c, t, s] = sim.getCompCount(compartment_names[c], specie_names[s])
   
    # return the data json encoded
    data = {'species':specie_names, 'compartments':compartment_names, 'data':res.tolist()}
    sys.stdout.write(json.dumps(data))
    sys.stdout.flush()
    sys.exit(0)
    # done run

# ------------------------------
# internal functions to handle
# the incoming model
# -----------------------------

def handleEvent(s_string, sim, model, spec_names, comp_names, rate_names, reac_to_comp, sreac_to_patch, kcst_si_factor):
    """ Parse and convert to executable Python code 
    of a string definig a Matlab Simbiology event. 

    We assume the event string is very close to executable Python code.
    This seems reasonable since we put it together in the Matlab part.
    """
    p = ExprParser(spec_names, comp_names, rate_names, reac_to_comp, sreac_to_patch, kcst_si_factor)

    if not type(s_string) == type([]):
        s_string = [s_string]

    e_code = ''
    for s in s_string:
        logger.debug('handleEvent event string: ' + s)
        e_code += p.toPython(str(s), generate=True)
        logger.debug('handleEvent executable python code: ' + e_code)
    return(e_code)
    # done handleEvent

def validSpecie(spec, model):
    """Returns true is spec is a
       valid species in the model,
       false otherwise.
    """
    species_list = model.getAllSpecs()
    
    # check if species exists in the STEPS model
    for s in species_list:
        if spec == s.getID():
            return(True)
 
    return(False)
    # done validSpecie


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Call STEPS from within Matlab.')
    parser.add_argument('model_string', help='STEPS model as a string.')
    parser.add_argument('initial_condition_string', help='The initial conditions as string.')
    parser.add_argument('stop_time', type=float, help='The final time until when the simulation should be run.')
    parser.add_argument('--events', help='JSON encoded distionary cantaining triggers and events', default='')
    parser.add_argument('--seed', type=int, help='PRNG seed', default=23412)
    parser.add_argument('--solver', help='solver name', default='wmdirect')
    parser.add_argument('--dt', type=float, default=1)
    parser.add_argument('--rk4dt', type=float, default=1e-5)
    parser.add_argument('--realisations', type=int, default=1)
    parser.add_argument('--file', action='store_true', default=False) 

    args = parser.parse_args()
    run(args.model_string, args.initial_condition_string, args.stop_time, args.events, args.seed, args.solver, args.dt, args.rk4dt, args.realisations, args.file)
    
