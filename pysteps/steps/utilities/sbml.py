####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

try:
    import libsbml
except:
    print "Unable to import SBML"
import sys
import os
import math
import random

import steps.model as smodel
import steps.geom as sgeom

####################################################################################################
#                                           Constants                                              #
####################################################################################################

# Avogadro's number. Source: http://physics.nist.gov/cgi-bin/cuu/Value?na 23/03/2011
AVOGADRO = 6.02214179e23

####################################################################################################
#                                       Auxilliary functions                                       #
####################################################################################################

def is_num(n):
    """
    Returns whether argument is a number or not
    """
    return isinstance(n, (float, int, long, complex))

def _float_approx_equal(x, y, tol=1e-32, rel=1e-7):
    """
    Return whether two floats are approximately equal
    """
    
    # NOTE: tolerance of 1e-16 was allowing some clearly different rates
    # to be equated as equal, due to the small random numbers used now. 
    # A very low absolute tolerance is required
    
    if tol is rel is None:
        raise TypeError('cannot specify both absolute and relative errors as None')
    tests = []
    if tol is not None: tests.append(tol)
    if rel is not None: tests.append(rel*abs(x))
    assert tests
    return abs(x - y) <= max(tests)
    
####################################################################################################

def MLfunc(level, refs):
    """
    Returns the result of a function imported from mathML -> python list
    
    Arguments:
        level: the function in list form
        refs: a dictionary of any names -> [value (concs, time , species, params), factor]
            NOTE::  values will be given in SBML model units, but they should be converted to given
                    units (which may be anything) as they are supposedto be if these differ.
                    SBML value / factor gives the quantity in given units. 
    NOTE: Return value should be in given units for whatever it is. 
    """
    
    # This is why we need to do the units conversion, the number should be in given units
    if (is_num(level)): return level 
        
    if (isinstance(level, (str))):
        if (level in refs): 
            # Converting to given units 
            if refs[level][1]: return refs[level][0]/refs[level][1]
            else: return refs[level][0]
        else:
            raise NotImplementedError("Unknown parameter '%s' in maths."%level)       
    
    if (not level or (len(level) != 2) or (len(level) != 2)):
        raise NotImplementedError("Symbol in maths not a number, a known variable, or an expression.")
        
    chrct = level[0]
    
    if(chrct == '+'):
        return MLfunc(level[1][0], refs) + MLfunc(level[1][1], refs)
    elif(chrct == '*'):
        return MLfunc(level[1][0], refs) * MLfunc(level[1][1], refs)
    elif(chrct == '-'):
        return MLfunc(level[1][0], refs) - MLfunc(level[1][1], refs)
    elif(chrct == '/'):
        if  MLfunc(level[1][1], refs) == 0.0: return 0.0
        else: return MLfunc(level[1][0], refs) / MLfunc(level[1][1], refs)    
    elif (chrct == 'power'):
        return math.pow(MLfunc(level[1][0], refs), MLfunc(level[1][1], refs))
    elif (chrct == 'root'):
        return math.pow(MLfunc(level[1][0], refs), (1.0/MLfunc(level[1][1], refs)))
    elif (chrct == 'lt'):
        return (MLfunc(level[1][0], refs) < MLfunc(level[1][1], refs))
    elif (chrct == 'leq'):
        return (MLfunc(level[1][0], refs) <= MLfunc(level[1][1], refs))
    elif (chrct == 'gt'):
        return (MLfunc(level[1][0], refs) > MLfunc(level[1][1], refs))
    elif (chrct == 'geq'):
        return ((MLfunc(level[1][0], refs)) >= (MLfunc(level[1][1], refs)))
    elif (chrct == 'eq'):
        return (MLfunc(level[1][0], refs) == MLfunc(level[1][1], refs))
    elif (chrct == 'neq'):
        return (MLfunc(level[1][0], refs) != MLfunc(level[1][1], refs))
    elif (chrct == 'exp'):
        return math.exp(MLfunc(level[1][0], refs))
    elif (chrct == 'abs'):
        return abs(MLfunc(level[1][0], refs))
    elif (chrct == 'ln'):
        return math.log(MLfunc(level[1][0], refs))
    elif (chrct == 'log'):
        return math.log10(MLfunc(level[1][0], refs))
    elif (chrct == 'floor'):
        return math.floor(MLfunc(level[1][0], refs))
    elif (chrct == 'ceiling'):
        return math.ceil(MLfunc(level[1][0], refs))
    elif (chrct == 'piecewise'):
        # The length of the list /2 will give an integer - the number of piecwise tests
        for i in range(len(level[1])/2):
            if (MLfunc(level[1][(i*2)+1], refs)): 
                return MLfunc(level[1][i*2], refs)
        return MLfunc(level[1][-1], refs)
    elif (chrct == 'and'):
        return (MLfunc(level[1][0], refs) and MLfunc(level[1][1], refs))
    elif (chrct == 'or'):
        return (MLfunc(level[1][0], refs) or MLfunc(level[1][1], refs))
    elif (chrct == 'not'):
        return (not MLfunc(level[1][0], refs))
    elif (chrct == 'sin'):
        return math.sin(MLfunc(level[1][0], refs))
    elif (chrct == 'cos'):
        return math.cos(MLfunc(level[1][0], refs))
    elif (chrct == 'tan'):
        return math.tan(MLfunc(level[1][0], refs))
    elif (chrct == 'sinh'):
        return math.sinh(MLfunc(level[1][0], refs))
    elif (chrct == 'cosh'):
        return math.cosh(MLfunc(level[1][0], refs))
    elif (chrct == 'tanh'):
        return math.tanh(MLfunc(level[1][0], refs))
    elif (chrct == 'arcsin'):
        return math.asin(MLfunc(level[1][0], refs))
    elif (chrct == 'arccos'):
        return math.acos(MLfunc(level[1][0], refs))
    elif (chrct == 'arctan'):
        return math.atan(MLfunc(level[1][0], refs))
    else: 
        raise NotImplementedError("Unknown character in maths.")

####################################################################################################

def rate_func(level, refs, params_refs, return_list, bad_params):
    """
    Returns the result of a function imported from mathML -> python list
    Note::  Seperate function needed specifically for rate functions, 
            which should have a specific form. No 'piecwise' for example.
    
    Arguments:
        level: the function in list form
        refs: a dictionary of any names -> [value (concs, time , species, params), factor]
        return_list: to return a list of used parameters
        bad_params: a list of things that are known not to be reaction parameters.
    """
    
    # This is why we need to do the units conversion, the number is in SBML units
    if (is_num(level)): return level
    
    if (isinstance(level, (str))):
        
        # Dividing by factor to keep in SBML units
        if (level in refs): 
            if (level in params_refs): 
                return_list.append(level)
            # Converting to given units 
            if refs[level][1]: return refs[level][0]/refs[level][1]
            else: return refs[level][0]
        else:
            if (level in bad_params):
                raise NotImplementedError("Parameter '%s' in reaction maths has unexpected units."%level)
            else:
                raise NotImplementedError("Unknown or unsupported parameter '%s' in reaction maths."%level)
    
    if (not level or (len(level) != 2) or (len(level) != 2)):
        raise NotImplementedError("Symbol in reaction maths not a number, a known variable, or an expression.")
    
    chrct = level[0]
    
    if(chrct == '+'):
        return rate_func(level[1][0], refs, params_refs, return_list, bad_params) \
            + rate_func(level[1][1], refs, params_refs, return_list, bad_params)
    elif(chrct == '*'):
        return rate_func(level[1][0], refs, params_refs, return_list, bad_params) \
            * rate_func(level[1][1], refs, params_refs, return_list, bad_params)
    elif(chrct == '-'):
        return rate_func(level[1][0], refs, params_refs, return_list, bad_params) - \
            rate_func(level[1][1], refs, params_refs, return_list, bad_params)
    elif(chrct == '/'):
        return rate_func(level[1][0], refs, params_refs, return_list, bad_params) / \
            rate_func(level[1][1], refs, params_refs, return_list, bad_params)    
    elif (chrct == 'power'):
        return math.pow(rate_func(level[1][0], refs, params_refs, return_list, bad_params), \
            rate_func(level[1][1], refs, params_refs, return_list, bad_params))
    else: 
        raise NotImplementedError("Unsupported form of reaction rate maths.")

####################################################################################################

# Expects a list of id strings and will return a list of the multiples 
# Horrible complication that the species may be amount or conc. If amount convert to a conc
# Spec_comp should be id of compartment

def make_rate_list(varibs, spec_subs_units, spec_comp):
    """
    Creates a reaction maths list of the form expected by STEPS. 
    If this doesn't match the reaction maths in SBML the import will fail. 
    """
    if (len(varibs) == 0): return 1
    var = varibs.pop()
    if (var in spec_subs_units):
        if spec_subs_units[var] == True:
            var_list = ['/', [var,spec_comp]]
        else: var_list = var
    else: var_list = var
    
    if (len(varibs) == 1): 
        if (varibs[0] in spec_subs_units):
            if (spec_subs_units[varibs[0]] == True): 
                var_list0 = ['/', [varibs[0],spec_comp]]
                return ['*', [var_list0, var_list]]
            else: return ['*', [varibs[0], var_list]]
        else: return ['*', [varibs[0], var_list]]
    else: 
        return ['*', [var_list, make_rate_list(varibs, spec_subs_units, spec_comp)]]

####################################################################################################

def MLtoList(math_obj, function_defs = {}, replace_dict={}):
    """
    Converts a mathMl object to a more readable Python layered list. Easiest to explain with examples:
    
    1 + 2  ->  ['+', [1,2]]
    
    1 + (2 * 3)  ->  ['+', [1, ['*', [2,3] ] ] ]
    
    etc etc
    
    Replace_dict is an optional dictionary of function arguments to real arguments, i.e model variables
    This is useful for lambda functions
    """
    
    if (math_obj.isReal()): 
        return math_obj.getReal()
    elif (math_obj.isInteger()): 
        return math_obj.getInteger() 
    elif (math_obj.isName()): 
        mo_name = math_obj.getName()
        if (mo_name in replace_dict): 
            return replace_dict[mo_name]
        else: 
            return mo_name
    elif (math_obj.isOperator()):
        if (math_obj.getNumChildren() == 2):
            return [math_obj.getCharacter(), [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict), \
                MLtoList(math_obj.getRightChild(), function_defs, replace_dict)]]
        elif (math_obj.getNumChildren() == 1):
            return [math_obj.getCharacter(), [0.0, MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        else: assert(False)
    elif (math_obj.isRelational()):
        assert(math_obj.getNumChildren() == 2)
        return [math_obj.getName(), [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict), \
            MLtoList(math_obj.getRightChild(), function_defs, replace_dict)]]
    elif (math_obj.isFunction()):
        fname = math_obj.getName()
        # Add some built-in functions here. BUILD ON THIS
        if (fname == 'power'):
            # Here the left child is the expression and the right child is the power, may be integer or other
            return ['power', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict), \
                MLtoList(math_obj.getRightChild(), function_defs, replace_dict)]]   
        elif (fname == 'root'):
            # Here the left child is the expression and the right child is the power, may be integer or other
            return ['root', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict), \
                MLtoList(math_obj.getRightChild(), function_defs, replace_dict)]] 
        elif (fname == 'exp'):
            # Here the left child is the variable
            return ['exp', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]] 
        elif (fname == 'abs'):
            # Here the left child is the variable
            return ['abs', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'ln'):
            # Here the left child is the variable
            return ['ln', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'log'):
            # Here the left child is the variable
            return ['log', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'floor'):
            # Here the left child is the variable
            return ['floor', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]        
        elif (fname == 'ceiling'):
            # Here the left child is the variable
            return ['ceiling', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]] 
        elif (fname == 'piecewise'):
            nchildren = math_obj.getNumChildren()
            if (nchildren %2 != 1):
                raise NotImplementedError("Piecewise function has unexpected number of arguments (%d)."%nchildren)
            templist = []
            for i in range(nchildren):
                templist.append(MLtoList(math_obj.getChild(i), function_defs, replace_dict))
            return ['piecewise', templist]
        elif (fname == 'sin'):
            return ['sin', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'cos'):
            return ['cos', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'tan'):
            return ['tan', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'sinh'):
            return ['sinh', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'cosh'):
            return ['cosh', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'tanh'):
            return ['tanh', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'arcsin'):
            return ['arcsin', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'arccos'):
            return ['arccos', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (fname == 'arctan'):
            return ['arctan', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        elif (math_obj.getName() not in function_defs):
            raise NotImplementedError("Function '%s' not known."%math_obj.getName())
        arglist = []
        n_children  = math_obj.getNumChildren()
        for n in range(n_children):
            # Arguments might contain some maths
            arglist.append(MLtoList(math_obj.getChild(n), function_defs, replace_dict))
        # Just so we're clear we expect a lambda object here
        lambda_obj = function_defs[math_obj.getName()]
        if (lambda_obj.isLambda() == False):
            raise NotImplementedError("Expected lambda function")
        nchildren = lambda_obj.getNumChildren()
        # The last child should be the math as far as I can tell
        lambda_math = lambda_obj.getChild(nchildren-1)
        # We expect the variables list and number of lambda arguments to be equal 
        # Create an object to get from labmda arg to real sim variable
        larg_to_list = {}
        for i in range(nchildren-1):
            larg_to_list[lambda_obj.getChild(i).getName()] = arglist[i]
        return MLtoList(lambda_math, function_defs, larg_to_list)
    elif (math_obj.isLogical()):
        fname = math_obj.getName()
        if (fname == 'and'):
            return ['and', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict), \
                MLtoList(math_obj.getRightChild(), function_defs, replace_dict)]]  
        elif (fname == 'or'):
            return ['or', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict), \
                MLtoList(math_obj.getRightChild(), function_defs, replace_dict)]] 
        elif (fname == 'not'):
            return ['not', [MLtoList(math_obj.getLeftChild(), function_defs, replace_dict)]]
        else:
            raise NotImplementedError("Math object unknown logical operator: '%s'"%math_obj.getName())
    elif (math_obj.isConstant()):
        # TODO: Add to this, what other constants are possible?
        cname = math_obj.getName()
        if (cname == 'pi'):
            return math.pi
        elif (cname == 'false'):
            return False
        elif (cname == 'true'):
            return True
        elif (cname == 'exponentiale'):
            return math.e
        else:
            raise NotImplementedError("Math object unknown constant: '%s'"%cname)
    else: 
        raise NotImplementedError("Math object unknown type: '%s'"%math_obj.getName())            
  

####################################################################################################

# Convert units to base SI units used by STEPS. Notable exception to the rule in STEPS that concentration are molar units
# i.e. bsaed on litres not cubic metres
def unitsDef_to_STEPS(mult, sca, exp):
    """
    Convert units to base SI units used by STEPS. Notable exception to the rule in STEPS that 
    concentration are molar units, i.e. based on litres not cubic metres, but this appears to be 
    the assumption in SBML too.
    """
    return (math.pow((mult *math.pow(10, sca)), exp))


####################################################################################################
####################################################################################################

class Interface(object):
    """
    The interface to the SBML importer.
    """
    
    def __init__(self, filename, timeunits_def = 1.0, volunits_def = 1.0e-3, subsunits_def = 1.0, volume_def = False, area_def = False, strict_mode = False):
        """
        Construction::
        
            iSbml = steps.utilities.sbml.Interface(sbmlFile, defvolunits_litre = True)
            
        Construct an SBML interface object, and optionally declare whether default volume units are
        litres rather than cubic meters. 
        
        Arguments (required):
            * string sbmlFile
            
        Arguments (optional)
            * float timeunits_def (default = 1.0)
            * float volunits_def (default = 1.0e-3)
            * float subsunits_def (default = 1.0)
            * float volume_def (default = False)
            * float area_def (default = False)
            * bool strict_mode (default = False)
        
        The 'units' arguments are a means to set the default model units and unit itself shoud relate to s.i units.
        For example: to set default model units to ms, timeunits_def = 1.0e-3. To set default model substance units to micromol, subsunits_def=1.0e-6.
        NOTE: These default units will only be used when units are not explicitly declared in the SBML file.
        
        The volume_def and area_def arguments are a way to set volume of compartments or area of patches with size 1.0 and no units,
        which are very common in SBML files. This may allow for a model to be run stochastically without modifying the file. 
        Values should be given in s.i. units, based on metres.
        For example: to set default volume to 100 femtolitre:  volume_def = 100*1.0e-15*1.0e-3
        
        strict_mode is a bool. False mean that all reactions will be allowed, True means that only fundamental reactions that can be represented in 
        the SSA without modification will be allowed, the import failing if any other types of reactions are present. The vast majority of SBML
        models include some non-fundamental reactions so setting this to True will mean that most SBML imports will fail. 
                
        """
        self.__reader = libsbml.SBMLReader()
        self.__document = self.__reader.readSBML(filename)
        self.__model = self.__document.getModel()
        if self.__document.getNumErrors() > 0 :
            for err in range(self.__document.getNumErrors()):
                pm = self.__document.getError(err)
                raise IOError(pm.getMessage()) 
        
        # Need to save for reset() function
        self.__timeunits_def = timeunits_def
        self.__volunits_def = volunits_def
        self.__subsunits_def = subsunits_def
        
        self.__defvolume = volume_def
        self.__defarea = area_def
        
        self.__strict_mode = strict_mode
        
        # Get the time structure, with units, and the volume units from the model in one fell swoop
        # NOTE: I am assuming the SBML documentation is wrong here and that volume units has a bearing
        # on reactions other than order 1
        # The 'units' variables are for conversion to STEPS (s.i.) units: 
        #   value in SBML units * model units => value in STEPS units
        self.__time, self.__time_units, self.__volume_units, self.__substance_units, self.__area_units, self.__length_units = self._parseModelUnits() 
        
        # Sets 
        # __globalParamters: a dictionary with elements list len<2>; {'parameter_id': [value, factor]
        #       Value is, as will be standard, in SBML units
        #       Value / factor => value in given units 
        # __glob_params_order: a dictionary mapping: {'parameter_id': order}
        #       Order may be 'unknown', if units are not specified or parameter is not a reaction constant
        #
        self.__globalParameters, self.__glob_params_order = self._parseGlobalParameters()
        
        # Basic function parsing. May be used in lots of places in the model.
        self.__function_defs = self._parseFunctionDefs()
        
        # Returns a dictionaray with elements list len<2>; {'compr_id': [value, factor]}
        # Value is, as will be standard, in SBML units 
        # SBML value / factor -> value in given units
        # ie. if the comp volume is specified as 1 liter, but SBML volume units are 
        # m^3, comp value will be 0.001(m^3) and factor will be 0.001 too.   
        self.__comps, self.__patches = self._parseComps(volume_def, area_def)
        
        # Create the geometry and model container objects in STEPS
        self.__geom = sgeom.Geom()
        self.__mdl = smodel.Model()
        
        # Create a volume system for each compartment, and store dictionary- {comp_id:volsys_id}
        # volsys_id will just be comp_id + 'volsys'
        self.__comp_volsys = {}
        
        # Similarly store surface systems for patch
        self.__patch_surfsys={}
        
        # Simply creates the compartments in STEPS and adds 'volsys' volume system.
        self._setupComps()
        
        # NOTE: Unfortunately we don't have any connectivity information in SBML, 
        # so can't set up the patches yet
        
        # A lot of information about the species in the model, so stored in separate objects:
        #
        # self.__species is a dictionary mapping: {'species_id': [value, factor]}
        #   however the value can be concentration or amount: 
        #       If concentration; SBML value (substance units/volume or area units) / factor => value in given units
        #       If amount; SBML value (substance units) / factor => value in given units
        # 
        # self.__species_amount_flag: Used to be important but now kind of defunct:
        #   is True if quanitity is Amount, False if Concentration
        #
        # self.__species_comps: maps species id to the id of its compartment in SBML. 
        #   NOTE: species can only belong to one compartment in SBML, not the case in STEPS
        #
        # self.__species_const_bc: maps id to tuple of two booleans: (constant, boundary_condition)
        #
        # NOTE: The values, which may be concentrations or amounts, are Molar for volume, amount/m^2 for area
        # (species may belong to 3D or 2D compartment)
        self.__species, self.__species_amount_flag, self.__species_subst_units,\
            self.__species_comps, self.__species_const_bc, = self._parseSpecies()
        
        # Create the chemical species in STEPS
        # NOTE: This must come after _setupComps  because it needs acces to the volume systems and surface systems
        self._setupModel1()
        
        # Set Intial Assignments where possible. 
        #   NOTE: This needs to come before _setupReactions() because Initial Assignments can change reaction constants
        #   NOTEL: Create dictionary for species references- available for _parseInitialAssignments function 
        self.__spec_refs = {}
        self._parseInitialAssignments()
 
        # Create assignment and rate rule structures 
        #   NOTE: This must come before _setupReactions is called for obvious reasons - not necessary now because of 
        #       the new method of updating concentrations manually and not trying to create reactions
        #   NOTE: This must also come before _parse_Reactions, because Assignment rules can contain species reference
        #
        # __rules_ass; dictionary with list elements len<2>: {assignmentrule_id: [var_id, assignment_math]}
        # __rules_rate; dictionary with list elements len<2>: {raterule_id: [var_id, rate_math]}
        #
        self.__rules_ass, self.__rules_rate = self._parseRules()
        
        # Create local reaction objects of class Reaction
        #   reaction constant from SBML units to STEPS units
        self.__reactions, self.__math_reactions, self.__surface_reactions, self.__surface_math_reactions = self._parseReactions()
        
        # Now create the patches now that we have the connectivity information
        self._setupPatches()
        
        # Create the chemical species in STEPS
        # NOTE: This must come after _setupPatches because it needs acces to the surface systems
        self._setupModel2()

        # Create the steps.model.Reac objects from local Reaction objects
        self._setupReactions()
        
        # Read the events and store the information in separate objects:
        #
        # __evnts_lt; dictionary with list elements len<2>: {event_id: [left_math, right_math]}
        #   These events will fire when the value of the lhs is LESS than the right hand side
        #
        # __evnts_gt; dictionary with list elements len<2>: {event_id: [left_math, right_math]}
        #   These events will fire when the value of the lhs is GREATER THAN the right hand side      
        # 
        # __evnts_ass; The event assignemnts
        # __evnts_dl; The event delay maths and time: [math, time]. The time will be calcualted from the math at time of trigger
        # __evnts_flip; TO store if event has gone from 'True' to 'False'  
        
        self.__evnts_trig, self.__evnts_ass, self.__evnts_dl, self.__evnts_flip, self.__evnts_trigvals = self._parseevnts()
        
        # Object that will temporarily store the event to fire, as a means to allow for delays
        self.__evnts_fire = {}
        # And an object to store values if they are to be evaluated at trigger time
        self.__evnts_vals = {}
        
        # Store hard-value of parameter conversions to save recalculating on the fly everytime:
        self.__param_converter_vol={}
        self.__param_converter_area={}
        # Go from zero to 5
        for order in range(0, 6):
            self.__param_converter_vol[order] = (math.pow(self.__volume_units, order-1)*(math.pow(1000, order-1)))/(self.__time_units*math.pow(self.__substance_units, order-1))
            self.__param_converter_area[order] = math.pow(self.__area_units, order-1)/(self.__time_units*math.pow(self.__substance_units, order-1))
            
    ################################################################################################
    
    def reset(self):
        """
        Resets the model to initailly imported state.
        This function basically repeats the import as that is the simplest and safest way.
        """
        self.__time, self.__time_units, self.__volume_units, self.__substance_units, self.__area_units, self.__length_units = self._parseModelUnits() 
        self.__globalParameters, self.__glob_params_order = self._parseGlobalParameters()
        self.__function_defs = self._parseFunctionDefs()
        self.__comps, self.__patches = self._parseComps(self.__defvolume, self.__defarea)
        self.__geom = sgeom.Geom()
        self.__mdl = smodel.Model()
        self.__comp_volsys = {}
        self.__patch_surfsys={}
        self._setupComps()
        self.__species, self.__species_amount_flag, self.__species_subst_units,\
            self.__species_comps, self.__species_const_bc, = self._parseSpecies()
        self._setupModel1()
        self.__spec_refs = {}
        self._parseInitialAssignments()
        self.__rules_ass, self.__rules_rate = self._parseRules()
        self.__reactions, self.__math_reactions, self.__surface_reactions, self.__surface_math_reactions = self._parseReactions()
        self._setupPatches()
        self._setupModel2()
        self._setupReactions() 
        self.__evnts_trig, self.__evnts_ass, self.__evnts_dl, self.__evnts_flip, self.__evnts_trigvals = self._parseevnts()
        self.__evnts_fire = {}
        self.__evnts_vals = {}

    ################################################################################################
        
    def _parseFunctionDefs(self):
        """
        Import the function definitions. 
        """
        fdefs = {}
        for funcdef in self.__model.getListOfFunctionDefinitions():
            fdefs[funcdef.getId()] = funcdef.getMath()
        return fdefs
    
    ################################################################################################
        
    def _parseModelUnits(self):
        """
        Attempt to read model units for time, volume and substance, though often absent. 
        """
        
        ret_time = {}
        ret_time_units = 1.0
        unitsTime = self.__model.getTimeUnits()
        if not unitsTime: unitsTime = 'time'
        unitdef = self.__model.getUnitDefinition(unitsTime)
        if unitdef: 
            units = unitdef.getListOfUnits()
            if (len(units) == 1):
                unit=units[0]
                if (unit.isSecond()):
                    expo = unit.getExponent()
                    scale = unit.getScale()
                    multiplier = unit.getMultiplier()
                    factor = unitsDef_to_STEPS(multiplier, scale, expo)
                    # Add all the different symbols I found in semantic tests
                    ret_time['time'] = [0.0, 1.0]
                    ret_time['t'] = [0.0, 1.0] 
                    ret_time['s'] = [0.0, 1.0] 
                    ret_time['Time'] = [0.0, 1.0]
                    ret_time_units = factor
                    
                else:
                    raise NotImplementedError("Time base units are not seconds!")
            else:
                print "WARNING: Failed to read model unit for time. Model time unit set to default (%s seconds)"%str(self.__timeunits_def)
                ret_time['time'] = [0.0, 1.0]
                ret_time['t'] = [0.0, 1.0]
                ret_time['s'] = [0.0, 1.0]
                ret_time['Time'] = [0.0, 1.0]
                ret_time_units = self.__timeunits_def
                
        else:
            print "WARNING: Failed to read model unit for time. Model time unit set to default (%s seconds)"%str(self.__timeunits_def)
            ret_time['time'] = [0.0, 1.0]
            ret_time['t'] = [0.0, 1.0]
            ret_time['s'] = [0.0, 1.0]
            ret_time['Time'] = [0.0, 1.0]
            ret_time_units = self.__timeunits_def
        
        ret_vol_units = 0
        unitsVol = self.__model.getVolumeUnits()
        if not unitsVol: unitsVol = 'volume'
        unitdef = self.__model.getUnitDefinition(unitsVol)
        if unitdef: 
            units = unitdef.getListOfUnits()
            if (len(units) == 1):
                unit = units[0]
                expo = unit.getExponent()
                scale = unit.getScale()
                multiplier = unit.getMultiplier()
                if (unit.isLitre()):
                    if (expo != 1): raise NotImplementedError("Volume units in model are not volume units!")
                    ret_vol_units = unitsDef_to_STEPS(multiplier, scale, expo)*0.001
                elif (unit.isMetre()):
                    if (expo != 3): raise NotImplementedError("Volume units in model are not volume units!")
                    ret_vol_units = unitsDef_to_STEPS(multiplier, scale, expo)
                else: raise NotImplementedError("Volume units in model are not volume units!")
            else:
                print "WARNING: Failed to read model unit for volume. Model volume unit set to default (%s m^3)"%str(self.__volunits_def)
                ret_vol_units =  self.__volunits_def
        else:
            print "WARNING: Failed to read model unit for volume. Model volume unit set to default (%s m^3)"%str(self.__volunits_def)
            ret_vol_units =  self.__volunits_def
        
        ret_area_units = 0
        unitsArea = self.__model.getAreaUnits()
        if not unitsArea: unitsArea = 'area'
        unitdef = self.__model.getUnitDefinition(unitsArea)
        if unitdef: 
            units = unitdef.getListOfUnits()
            if (len(units) == 1):
                unit = units[0]
                expo = unit.getExponent()
                scale = unit.getScale()
                multiplier = unit.getMultiplier()
                if (unit.isMetre()):
                    if (expo != 2): raise NotImplementedError("Area units in model are not area units!")
                    ret_area_units = unitsDef_to_STEPS(multiplier, scale, expo)
                else: raise NotImplementedError("Area units in model are not area units!")
            else:
                # On second thoughts, don't think a warning is necessary here
                #print "WARNING: Failed to read model unit for area. Default unit will be based on model volume units."
                ret_area_units =  math.pow(ret_vol_units, 2.0/3.0)
        else:
            # On second thoughts, don't think a warning is necessary here
            # print "WARNING: Failed to read model unit for area. Default unit will be based on model volume units."
            ret_area_units =  math.pow(ret_vol_units, 2.0/3.0)

        ret_length_units = 0
        unitsLength = self.__model.getLengthUnits()
        if not unitsLength: unitsLength = 'length'
        unitdef = self.__model.getUnitDefinition(unitsLength)
        if unitdef: 
            units = unitdef.getListOfUnits()
            if (len(units) == 1):
                unit = units[0]
                expo = unit.getExponent()
                scale = unit.getScale()
                multiplier = unit.getMultiplier()
                if (unit.isMetre()):
                    if (expo != 1): raise NotImplementedError("Length units in model are not length units!")
                    ret_length_units = unitsDef_to_STEPS(multiplier, scale, expo)
                else: raise NotImplementedError("Length units in model are not length units!")
            else:
                # On second thoughts, don't think a warning is necessary here
                #print "WARNING: Failed to read model unit for length. Default unit will be based on model volume units."
                ret_length_units =  math.pow(ret_vol_units, 1.0/3.0)                    
        else:
            # On second thoughts, don't think a warning is necessary here
            #print "WARNING: Failed to read model unit for length. Default unit will be based on model volume units."
            ret_length_units =  math.pow(ret_vol_units, 1.0/3.0)  
        
        ret_subs_units = 0
        unitsSubs = self.__model.getSubstanceUnits()
        if not unitsSubs: unitsSubs = 'substance'
        unitdef = self.__model.getUnitDefinition(unitsSubs)
        if unitdef:
            units = unitdef.getListOfUnits()
            if (len(units) == 1):
                unit = units[0]
                expo = unit.getExponent()
                scale = unit.getScale()
                multiplier = unit.getMultiplier()
                if (unit.isMole()):
                    ret_subs_units = unitsDef_to_STEPS(multiplier, scale, expo)
                elif (unit.isDimensionless() or unit.isItem()):
                    ret_subs_units = unitsDef_to_STEPS(multiplier, scale, expo)/AVOGADRO
                else: raise NotImplementedError("Substance units in model are not supported units.")
            else:
                print "WARNING: Failed to read model unit for substance. Model substance unit set to default (%s mole)"%str(self.__subsunits_def)
                ret_subs_units =  self.__subsunits_def
        else:
            print "WARNING: Failed to read model unit for substance. Model substance unit set to default (%s mole)"%str(self.__subsunits_def)
            ret_subs_units =  self.__subsunits_def
        
        return ret_time, ret_time_units, ret_vol_units, ret_subs_units, ret_area_units, ret_length_units
    
    ################################################################################################
    
    def _parseComps(self, defvolume, defarea):
        """
        Import the compartments.
        """
        
        # If user specified a default volume for compartments use that, otherwise just one model volume unit
        # Same for area
        if defvolume and is_num(defvolume): volume_default = defvolume/self.__volume_units
        else: volume_default = 1.0
        
        if defarea and is_num(defarea): area_default = defarea/self.__area_units
        else: area_default = 1.0
        
        ListOfComp = self.__model.getListOfCompartments()
        
        comps = {}
        patches = {}
        for comp in ListOfComp:
        
            idComp = str(comp.getId())
            
            dim = -1
            # It appears now in level 3 that a user doesn't even need to set the dimensions of a compartment. 
            # How I wish they would consider developers and not just modellers!!
            # This is basically useless now because it doesn't have to be specified so probably just ignore this test
            if (comp.isSetSpatialDimensions()):
                dim = comp.getSpatialDimensions()
                if (dim != 3 and dim != 2):
                    raise NotImplementedError("Compartments must be 3D or 2D ('patches').", \
                    "Compartment'%s' in this model %d dimensions." %(idComp, comp.getSpatialDimensions()))
            # What choice do we have? Try to find out from the units, or if not then assume 3D 
            
            # Need a bool to keep track if we have already found sufficient volume for comp
            added = False
            
            sizeComp = comp.getSize()
            
            # Discovered that volume may be nan
            if (sizeComp!=sizeComp):
                raise NotImplementedError("Compartment size is not-a-number", \
                                        " for compartment '%s'." %idComp)   
                # NOTE: Possibility here that compartment volume should be initialised by 
                # some assignment rule or something, but unfortunately in STEPS
                # we can't start with zero volume so I'm going to keep the condition that
                # the compartment has to be initialised with some volume
            elif (sizeComp):
                # Find the units:
                unitsComp = comp.getUnits()
                # Try to get from the model
                # We're going to assume 3D if not set, no other choice really
                if not unitsComp: 
                    unitsComp = self.__model.getVolumeUnits()
                if not unitsComp: 
                    if dim ==3: unitsComp = 'volume'
                    elif dim == 2: unitsComp = 'area'
                    else:
                        # What choice? Assume 3D for now? 
                        unitsComp = 'volume'
                if (unitsComp):
                    # Using a for loop that can be exited with a break, 
                    # for want of a better alternative
                    for j in [0]:
                        unitdef = self.__model.getUnitDefinition(unitsComp)
                        if not (unitdef):
                            break
                        if (unitdef.isVariantOfVolume()):
                            dim = 3
                            units = unitdef.getListOfUnits()
                            if (len(units) != 1):
                                print "WARNING: Failed to read units for compartment '%s'. " \
                                "Model volume units will be assumed."%idComp
                                break
                            unit = units[0]
                            if (unit.isMetre()):
                                if (unit.getExponent() != 3):
                                    raise NotImplementedError("Compartment '%s' dimensions not 3D."%idComp)
                                scale = unit.getScale()
                                multiplier = unit.getMultiplier()
                                # Convert to STEPS (s.i) units and then to SBML
                                factor = unitsDef_to_STEPS(multiplier, scale, 3)/self.__volume_units
                                sizeComp *= factor
                                comps[idComp] = [sizeComp, factor]
                                added = True
                                break
                            elif (unit.isLitre()):
                                if (unit.getExponent() != 1):
                                    raise NotImplementedError("Compartment '%s' dimensions not 3D."%idComp)
                                scale = unit.getScale()
                                multiplier = unit.getMultiplier()
                                # Convert to STEPS (s.i) units and then to SBML
                                factor = (unitsDef_to_STEPS(multiplier, scale, 1)*1.0e-3)/self.__volume_units
                                sizeComp *= factor
                                comps[idComp] = [sizeComp, factor]
                                added = True
                                break
                            else: 
                                print "WARNING: Failed to read units for compartment '%s'. " \
                                "Model volume units will be assumed."%idComp
                                break
                            break
                        elif(unitdef.isVariantOfArea()):
                            dim = 2
                            units = unitdef.getListOfUnits()
                            unit = units[0]
                            if (unit.isMetre()):
                                if (unit.getExponent() != 2):
                                    raise NotImplementedError("Patch '%s' dimensions not 2D."%idComp)
                                scale = unit.getScale()
                                multiplier = unit.getMultiplier()
                                # Convert to STEPS (s.i) units and then to SBML
                                factor = unitsDef_to_STEPS(multiplier, scale, 2)/self.__area_units
                                sizeComp *= factor
                                patches[idComp] = [sizeComp, factor]
                                added = True
                                break
                            else: 
                                print "WARNING: Failed to read units for patch '%s'. " \
                                "Model area units will be assumed."%idComp
                                break
                            break
                        else:
                            raise NotImplementedError("Compartment '%s' dimensions not 2D or 3D."%idComp)
                                            
                    # If we got here we breaked. 
            else:
                if dim == 3:
                    print "WARNING: No volume specified for compartment '%s'. " \
                        "Volume set to default value (%sm^3)."%(idComp, str(volume_default*self.__volume_units))
                    sizeComp = volume_default
                    comps[idComp] = [sizeComp, 1.0]
                    added = True
                    continue
                elif dim == 2:
                    print "WARNING: No area specified for patch '%s'. " \
                            "Area set to default value (%m^2s)."%(idComp, str(area_default*self.__area_units))
                    sizeComp = area_default
                    patches[idComp] = [sizeComp, 1.0]
                    added = True
                    continue
                else:
                    raise NotImplementedError("Failed to add compartment", \
                                        "'%s'. Could not find dimensions." %idComp)                      
            
            # This will make the default a volume compartment
            if dim != 2:        
                # Special case if size == 1 without units, a common occurunce in SBML
                if (comp.getSize() == 1.0 and not comp.getUnits()):
                    print "WARNING: Compartment '%s' size = 1.0 wih no units. " \
                        "Volume set to default value (%sm^3)."%(idComp, str(volume_default*self.__volume_units))
                    sizeComp = volume_default
                    
                    comps[idComp] = [sizeComp , 1.0]
                    added = True
                    continue  
                
                # Special case if size == 1 without default units (litre), also a common occurance
                if (comp.getSize() == 1.0 and idComp not in comps):
                    print "WARNING: Compartment '%s' size = 1.0 wih default units. " \
                        "Volume set to default value (%sm^3)."%(idComp, str(volume_default*self.__volume_units))
                    sizeComp = volume_default
                    comps[idComp] = [sizeComp , 1.0]
                    added = True
                    continue   
                
                if not added:
                    # If we got here we have a volume, which is not equal to 1.0, 
                    # but couldn't read the units. Assume model units
                    comps[idComp] = [sizeComp, 1.0]
                    added = True
                    continue
                    
            else:    
                # Special case if size == 1 without units, a common occurunce in SBML
                if (comp.getSize() == 1.0 and not comp.getUnits()):
                    print "WARNING: Patch '%s' size = 1.0 wih no units. " \
                            "Area set to default value (%m^2s)."%(idComp, str(area_default*self.__area_units))
                    sizeComp = area_default                    
                    patches[idComp] = [sizeComp , 1.0]
                    added = True
                    continue  
                
                # Special case if size == 1 without default units (litre), also a common occurance
                if (comp.getSize() == 1.0 and idComp not in patches):
                    print "WARNING: Patch '%s' size = 1.0 wih default units. " \
                            "Area set to default value (%m^2s)."%(idComp, str(area_default*self.__area_units))
                    sizeComp = area_default                    
                    patches[idComp] = [sizeComp , 1.0]
                    added = True
                    continue   
                
                if not added:
                    # If we got here we have an area, which is not equal to 1.0, 
                    # but couldn't read the units. Assume model units
                    patches[idComp] = [sizeComp, 1.0]
                    added = True
                    continue
        
        return comps, patches
  
    ################################################################################################

    def _setupComps(self):
        """ 
        Create the compartments in STEPS.
        """
        
        for c in self.__comps:
            # Values are now stored in SBML units, the conversion to STEPS is the model units of volume
            comp = sgeom.Comp(c, self.__geom, vol = self.__comps[c][0]*self.__volume_units)
            vsysid = c+'_volsys'
            comp.addVolsys(vsysid)
            # Doing the volsys here just to keep it in one place
            vsys = smodel.Volsys(vsysid, self.__mdl)
            self.__comp_volsys[c] = vsysid

    ################################################################################################
    
    
    def _setupPatches(self):
        """ 
        Create the patches in STEPS, as well as the surface systems.
        """
        
        # OK, this is the complicated part. 
        
        # Can't create the patches yet because there is no connectivity information in SBML. 
        # First create a dictionary of known patches to 'inner' and 'outer' comp, found from the surface reactions
        # If 2 reactions differ in their 'inner' and 'outer' from this description, one of them will have to be flipped. 
        # If they differ so much that an entirely new compartment is involved, that's a STEPS error
        
        # Dictionary: {patch_ids: [innercomp_id, outercomp_id]}
        patch_connectivity={}
        
        for p in self.__patches:
            patch_connectivity[p]=['','']
            
            # All 'math' reactions are included in surface_reaction, no need for an extra loop
            
            # Store a list of ones to flip, that is the inner comp matches the outer comp of the patch and vice-versa
            sr_flip = []
            for sr in self.__surface_reactions:
                sr_surf = sr.getSLhs()
                if not sr_surf: sr_surf = sr.getSRhs()
                if sr_surf: 
                    if self.__species_comps[sr_surf[0].getID()] == p:
                        sr_ilhs = sr.getILhs()
                        sr_olhs = sr.getOLhs()
                        sr_irhs = sr.getIRhs()
                        sr_orhs = sr.getORhs()
                        
                        sr_icomp = ''
                        sr_ocomp = ''
                                                
                        if sr_ilhs: sr_icomp = self.__species_comps[sr_ilhs[0].getID()]
                        if sr_olhs: sr_ocomp = self.__species_comps[sr_olhs[0].getID()]
                        # OK, we might end up getting them twice but so what?
                        if sr_irhs: sr_icomp = self.__species_comps[sr_irhs[0].getID()]
                        if sr_orhs: sr_ocomp = self.__species_comps[sr_orhs[0].getID()]
                        
                        # Lots of possibilities here. Both the surface reaction and the patch may 
                        # or may not have an inner and outer comp. If they are there and 
                        # match that is fine, but also if they are reversed (the inner comp of the SR
                        # mathces the outer comp of the patch etc).
                        
                        # Marker to see if we can use this surface reaction.
                        ok = False
                        
                        patch_icomp = patch_connectivity[p][0]
                        patch_ocomp = patch_connectivity[p][1]
                        if sr_icomp:
                            if sr_ocomp:
                                # Got both inner and outer comp for SR:
                                if patch_icomp:
                                    if patch_ocomp:
                                        # Got both inner and outer comp for patch
                                        if (sr_icomp == patch_icomp and sr_ocomp == patch_ocomp):
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            continue
                                        elif (sr_icomp == patch_ocomp and sr_ocomp == patch_icomp):
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            sr_flip.append(sr)
                                            continue
                                        else:
                                            raise NotImplementedError("Compartment connectivity cannot be represented", \
                                            "In STEPS for 2D compartment '%s'." %p)
                                    else:
                                        # Got inner comp for patch, but not outer
                                        if (sr_icomp == patch_icomp):
                                            patch_connectivity[p][1] = sr_ocomp
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            continue
                                        elif (sr_ocomp == patch_icomp):
                                            patch_connectivity[p][1] = sr_icomp
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            sr_flip.append(sr)
                                            continue
                                        else:
                                            raise NotImplementedError("Compartment connectivity cannot be represented", \
                                            "In STEPS for 2D compartment '%s'." %p)                                      
                                elif patch_ocomp:
                                        # Got outer comp for patch, but not inner
                                        if (sr_ocomp == patch_ocomp):
                                            patch_connectivity[p][0] = sr_icomp
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            continue
                                        elif (sr_icomp == patch_ocomp):
                                            patch_connectivity[p][0] = sr_ocomp
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            sr_flip.append(sr)
                                            continue
                                        else:
                                            raise NotImplementedError("Compartment connectivity cannot be represented", \
                                            "In STEPS for 2D compartment '%s'." %p)                                                                                          
                                else:
                                    # Got no comps for patch, easy
                                    patch_connectivity[p][0] = sr_icomp
                                    patch_connectivity[p][1] = sr_ocomp
                                    sr.setSsys(p+'_ssys')
                                    sr.setPatchID(p)
                                    continue
                            else:
                                # Got inner, but no outer comp for SR
                                if patch_icomp:
                                    if patch_ocomp:
                                        # Got both inner and outer comp for patch
                                        if (sr_icomp == patch_ocomp):
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            sr_flip.append(sr)
                                            continue
                                        elif (sr_icomp == patch_icomp):
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            continue
                                        else:
                                            raise NotImplementedError("Compartment connectivity cannot be represented", \
                                            "In STEPS for 2D compartment '%s'." %p)                                             
                                    else:
                                        # Got inner comp for patch, but not outer
                                        if (sr_icomp == patch_icomp):
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            continue
                                        else:
                                            sr.setSsys(p+'_ssys')
                                            sr.setPatchID(p)
                                            patch_connectivity[p][1] = sr_icomp
                                            sr_flip.append(sr)
                                            continue                                            
                                elif patch_ocomp:
                                    # Got outer comp for patch, but no inner
                                    if (sr_icomp == patch_ocomp):
                                        sr.setSsys(p+'_ssys')
                                        sr.setPatchID(p)
                                        sr_flip.append(sr)
                                        continue                                        
                                    else:
                                        patch_connectivity[p][0] = sr_icomp
                                        sr.setSsys(p+'_ssys')
                                        sr.setPatchID(p)
                                        continue
                                else:
                                    # Got no comps for patch, easy
                                    patch_connectivity[p][0] = sr_icomp
                                    sr.setSsys(p+'_ssys')
                                    sr.setPatchID(p)
                                    continue
                        elif sr_ocomp:
                            # Got only outer comp for SR
                            if patch_icomp:
                                if patch_ocomp:
                                    # Got both inner and outer comp for patch
                                    if (sr_ocomp == patch_ocomp):
                                        sr.setSsys(p+'_ssys')
                                        sr.setPatchID(p)
                                        continue    
                                    elif (sr_ocomp == patch_icomp):
                                        sr.setSsys(p+'_ssys')
                                        sr.setPatchID(p)
                                        sr_flip.append(sr)
                                        continue   
                                    else:
                                        raise NotImplementedError("Compartment connectivity cannot be represented", \
                                        "In STEPS for 2D compartment '%s'." %p)                                          
                                else:
                                    # Got inner comp for patch, but not outer
                                    if (sr_ocomp == patch_icomp):
                                        sr.setSsys(p+'_ssys')
                                        sr.setPatchID(p)
                                        sr_flip.append(sr)
                                        continue      
                                    else:
                                        sr.setSsys(p+'_ssys')
                                        sr.setPatchID(p)
                                        patch_connectivity[p][1] = sr_ocomp
                                        continue                                          
                            elif patch_ocomp:
                                # Got outer comp for patch, but no inner
                                if (sr_ocomp == patch_ocomp):
                                    sr.setSsys(p+'_ssys')
                                    sr.setPatchID(p)
                                    continue      
                                else:
                                    patch_connectivity[p][0] = sr_ocomp
                                    sr.setSsys(p+'_ssys')
                                    sr.setPatchID(p)
                                    sr_flip.append(sr)
                                    continue                                     
                            else:
                                # Got no comps for patch, easy  
                                # Has to be the 'inner' compartment to ensure inner is set
                                patch_connectivity[p][0] = sr_ocomp
                                sr_flip.append(sr)
                                sr.setSsys(p+'_ssys')
                                sr.setPatchID(p)
                                continue                                 
                        else:
                            # No compartments for surface reaction, don't care about connectivity
                            sr.setSsys(p+'_ssys')
                            sr.setPatchID(p)
                            continue     
                else: 
                    # Here, no surface species for surface reaction - there may or may not be a patch created in 
                    # SBML, so create a virtual one if not. Anyway, this is not the place to do it
                    continue
            
            # Now we are here with the 'flip-list'. Invoke method.
            for sr in sr_flip: sr.flip()
                
        # Now the volume to other volume reactions, where a virtual patch will have to be created if there isn't one there
        # A compartment can have as many patches as it likes, so create a virtual patch for all volume-volume reactions. 
                
        for sr in self.__surface_reactions:
            # If it has a surface system it's already been dealt with
            if sr.getSsys(): continue
            
            # Double-check
            sr_surf = sr.getSLhs()
            if not sr_surf: sr_surf = sr.getSRhs()
            assert(not sr_surf)
            
            sr_ilhs = sr.getILhs()
            sr_olhs = sr.getOLhs()
            sr_irhs = sr.getIRhs()
            sr_orhs = sr.getORhs()
            
            sr_icomp = ''
            sr_ocomp = ''
                                                
            if sr_ilhs: sr_icomp = self.__species_comps[sr_ilhs[0].getID()]
            if sr_olhs: 
                    sr_ocomp = self.__species_comps[sr_olhs[0].getID()]
            # OK, we might end up getting them twice but so what?
            if sr_irhs: sr_icomp = self.__species_comps[sr_irhs[0].getID()]
            if sr_orhs: 
                sr_ocomp = self.__species_comps[sr_orhs[0].getID()]    
            
            # Depending on boundary species trick volume-volume reaction may have had a species removed from
            # products, so may actually only now include inner comp. In that case use the virtaul comp  
            #assert (sr_icomp and sr_ocomp);
            assert (sr_icomp or sr_ocomp)
            
            # Dimensions completely make no difference
            self.__patches['virtual_patch_'+sr.getName()] = [1,1]
            sr.setSsys('virtual_patch_'+sr.getName()+'_ssys')
            sr.setPatchID('virtual_patch_'+sr.getName())
            
            patch_connectivity['virtual_patch_'+sr.getName()] = [sr_icomp, sr_ocomp]

        # Finally create the patches!
        # patch may have no connectivity if it involves purely surface surface reactions
        # Create a 'virtual comp' for this case
        virt_comp = sgeom.Comp('virtual_comp', self.__geom, vol = 1)
        for p in patch_connectivity:
            if not patch_connectivity[p][0]:
                assert (not patch_connectivity[p][1])
                patch_connectivity[p][0] = 'virtual_comp'
        
        for p in self.__patches:
            if patch_connectivity[p][1]: 
                # Area is now saved in SBML units
                patch = sgeom.Patch(p, self.__geom, icomp = self.__geom.getComp(patch_connectivity[p][0]), \
                    ocomp = self.__geom.getComp(patch_connectivity[p][1]),  area = self.__patches[p][0]*self.__area_units)       
            elif patch_connectivity[p][0]:
                patch = sgeom.Patch(p, self.__geom, icomp = self.__geom.getComp(patch_connectivity[p][0]), \
                    area = self.__patches[p][0]*self.__area_units)        
            else: assert(False)
            
            ssys_id = p + '_ssys'
            patch.addSurfsys(ssys_id)
            # Doing the volsys here just to keep it in one place
            ssys = smodel.Surfsys(ssys_id, self.__mdl)
            self.__patch_surfsys[p] = ssys_id

    ################################################################################################    
    
    def _parseSpecies(self):
        """
        Import the chemical species.
        """
        
        # Getting the various species
        ListOfSpecies = self.__model.getListOfSpecies()
        species = {}
        species_amount_flag = {}
        species_subst_units = {}
        species_comps = {}
        species_patches = {}
        species_const_bcs = {}
        
        factor = 1
        
        for specie in ListOfSpecies:
            # Any volume units should match the compartment units
            compId = specie.getCompartment()
            
            species_comps[str(specie.getId())] = compId
            
            vol_species = False
            surf_species = False
            if compId in self.__comps:
                vol_species = True
            elif compId in self.__patches:
                surf_species = True
            else:
                raise NotImplementedError("Species '%s' belongs to unknown compartment."%specie.getId())
            
            species_const_bcs[str(specie.getId())] = (specie.getConstant(), specie.getBoundaryCondition())
            
            if specie.isSetInitialAmount():
                value = specie.getInitialAmount()
                species_amount_flag[str(specie.getId())] = True
                isAmount = True
            elif specie.isSetInitialConcentration():
                value = specie.getInitialConcentration()
                species_amount_flag[str(specie.getId())] = False
                isAmount = False
            # Undefined so set to True - doesn't matter anyway
            else: 
                value = specie.getInitialAmount()
                species_amount_flag[str(specie.getId())] = True
                isAmount = True
            
            # THis is going to be the important flag. Internally, after this, the species
            # will follow these units regardless of what was injected initially
            if (specie.isSetHasOnlySubstanceUnits()): 
                hasAmountUnits = specie.getHasOnlySubstanceUnits()
            else: 
                hasAmountUnits = isAmount
            species_subst_units[specie.getId()] = hasAmountUnits
            
            # Get the units for the amount
            unitsSpec = specie.getSubstanceUnits()
            # Try to get from the model if not explicitly set
            if not unitsSpec: unitsSpec = self.__model.getSubstanceUnits()
            # Lets just try setting manually to substance then shall we?
            if not unitsSpec: unitsSpec = 'substance'
            # To convert any given species value, in rates, assignments etc, to the correct unit in STEPS
            factor = 1
            
            if (unitsSpec): 
                # Using a for loop that can be exited with a break, 
                # for want of a better alternative
                for j in [0]:
                    unitdef = self.__model.getUnitDefinition(unitsSpec) 
                    if not (unitdef):
                        #print "WARNING: Failed to read units of substance " \
                        #"for Species '%s'. Default units (mole) will be assumed."%specie.getId()
                        break
                    # Note: All units here should be of substance, i.e. dimensionless
                    # Units for the volume (if concentration) are taken from 
                    # the units for volume of the compartment
                    units = unitdef.getListOfUnits()
                    if (len(units) != 1):
                        print "WARNING: Failed to read units of substance " \
                        "for Species '%s'. Model substance unit will be assumed."%specie.getId()
                        break
                    unit = units[0]
                    if not (unit.isItem() or unit.isDimensionless() or unit.isMole()):
                        raise NotImplementedError("Species '%s' substance units are not supported units."%specie.getId())
                    if (unit.getExponent() != 1):
                        print "WARNING: Failed to read units of substance " \
                        "for Species '%s'. Model substance unit will be assumed."%specie.getId()
                        break                                
                    multiplier = unit.getMultiplier()
                    scale = unit.getScale()
                    # Convert to s.i. units, then SBML ones:
                    factor = unitsDef_to_STEPS(multiplier, scale, 1)/self.__substance_units
                    value *= factor
                    
                    # Convert to mole if we had something dimensionless
                    if (unit.isMole() != True): 
                        value /= AVOGADRO
                        factor /= AVOGADRO
                    break
            else: 
                print "WARNING: Failed to read units of substance " \
                    "for Species '%s'. Model substance unit will be assumed."%specie.getId()
            
            # Confusing, but value could still be an amount or conc, but the amount part is converted to 
            # model substance units. If a conc we still need to convert the volume part. 
            # Another complication is that we are going to set to conc or amount by the hasSubstanceUnits flag. 
            if (isAmount):
                if (hasAmountUnits): 
                    # Easy, we have an amount already converted to SBML units and we want to store an amount.
                    species[str(specie.getId())] = [value, factor] 
                else: 
                    # We have an amount in mole, but we want to store the concentration. Simply 
                    # divide by the comp volume (converted from^3 to litre) to get the conc in mole/litre.
                    # The units will be factor for the amount divided by the comp units, because it SBML conc -> STEPS conc
                    if (vol_species): 
                        # Volume units for conc are assumed to be compartment volume units
                        conc = value/(self.__comps[compId][0])
                        # The factor must convert to SBML volume units, which comaprtment is stored in:
                        species[str(specie.getId())] = [conc, factor/(self.__comps[compId][1])]
                    else:
                        conc = value/(self.__patches[compId][0])
                        species[str(specie.getId())] = [conc, factor/(self.__patches[compId][1])]
            else:
                if (hasAmountUnits):
                    # We've got a conc (the volume part still in compartment units), but have to convert to amount for expressions
                    # We have to multiply by the volume of the compartment converted to SBML units
                    if vol_species:
                        amount = value*(self.__comps[compId][0])
                        species[str(specie.getId())] = [amount, factor]
                    else:
                        amount = value*(self.__patches[compId][0])
                        # Units are just going to be factor
                        species[str(specie.getId())] = [amount, factor]                        
                else:
                    # We've got a conc and want to store a conc. Just need to convert the units for comp.
                    # Concentration units for compartment are stored in SBML units
                    if vol_species:
                        conc = value
                        # Still need to store the conversion from this conc to STEPS one and factor hs only had the substance part converted 
                        species[str(specie.getId())] = [conc, factor/(self.__comps[compId][1])]                        
                    else:
                        conc = value
                        species[str(specie.getId())] = [conc, factor/(self.__patches[compId][1])]                        
            
            # Crazily, sometimes the default amount for a species that is not declared in SBML is nan
            if (species[str(specie.getId())][0] != species[str(specie.getId())][0]): species[str(specie.getId())][0] = 0.0
            if (species[str(specie.getId())][1] != species[str(specie.getId())][1]): species[str(specie.getId())][1] = 0.0
        
        return species, species_amount_flag, species_subst_units, species_comps, species_const_bcs
    
    ################################################################################################

    def _getSpecies(self):
        return self.__species
        
    def _getSpeciesComps(self):
        return self.__species_comps
        
    ################################################################################################
    
    def _setupModel1(self):
        """ 
        Setup the model in STEPS.
        """
        
        for specie in self.__species:
            mol = smodel.Spec(specie, self.__mdl)

    ################################################################################################

    def _setupModel2(self):
        """ 
        Add the nothing reactions to comps and patches to make sure they are defined in them for rules.
        """
        
        for specie in self.__species:
            # This is a bit of a hack, but a few occasions when a species does not appear in
            # reactions, or stoichiometry changed becase of BCs so doesn't appear in reactions, 
            # and therefore not added to comp by STEPS.
            comp_id= self.__species_comps[specie]
            if comp_id in self.__comps:
                reac = smodel.Reac(specie+'reac', self.__mdl.getVolsys(self.__comp_volsys[comp_id]), lhs = [self.__mdl.getSpec(specie)], \
                    rhs = [self.__mdl.getSpec(specie)], kcst = 0.0)
            elif comp_id in self.__patches:
                sreac = smodel.SReac(specie+'sreac', self.__mdl.getSurfsys(self.__patch_surfsys[comp_id]), slhs = [self.__mdl.getSpec(specie)], \
                    srhs = [self.__mdl.getSpec(specie)], kcst = 0.0)  
    
    ################################################################################################
            
    def _parseReactions(self):
        """
        Import the reactions.
        """
       
        strict = self.__strict_mode
        
        listOfReactions = self.__model.getListOfReactions()
        
        reactions = []
        # The reactions that will be treated separately- different maths from expected
        math_reactions=[]
        
        # Reactions involving a 2D compartment, or two 3D compartments
        surface_reactions = []
        
        surface_math_reactions = []
        
        
        for reac in listOfReactions:
            
            reversible = reac.getReversible()
            
            # Reactants: list of ids of the reactants, separated for reversible reactions
            reacts_f = []
            if (reversible): reacts_b = []
            
            # Products: list of ids of the products, separated for reversible reactions  
            prods_f = []
            if (reversible): prods_b = []
            
            products = reac.getListOfProducts()            
            reactants = reac.getListOfReactants()
            
            totsto = 0
            for reactant in reactants:
                sto = reactant.getStoichiometry()
                stomath = reactant.getStoichiometryMath()
                if (stomath):
                    stonode = stomath.getMath()
                    # The stoichiometry could be set by complicated maths involving parameters
                    # elsewhere defined. Ours not to reason why...
                    slist = MLtoList(stonode, self.__function_defs)
                    sto = MLfunc(slist, self.__globalParameters)
                # Sto may be None or nan (thanks SBML)
                if (not sto or sto!=sto): 
                    # Hmm. Probably a species reference. 
                    # Bizarrely, .getId returns the species reference 
                    spec_ref = reactant.getId()
                    if spec_ref: sto = self.__spec_refs[spec_ref]
                    else:
                        # What choice? Assume 1
                        sto=1
                # stoichiometry must be a whole number (allowing a small error)
                if (sto % 1 > 0.001) : 
                    raise NotImplementedError("Partial stoichiometry (%f) in reaction '%s'."%(sto, reac.getId()))
                sto = int(round(sto))
                totsto+=sto
                if (totsto > 4): raise NotImplementedError("Stoichiometry (at least %f) in reaction '%s' is more than STEPS maximum of 4."%(totsto, reac.getId()))
                # More than one molecule of each reactant species may be present in the reaction
                for j in range(sto):
                    if(reactant.getSpecies() != "Empty"):
                        reacts_f.append(reactant.getSpecies())
                        if (reversible): prods_b.append(reactant.getSpecies())
            
            for product in products:        
                sto = product.getStoichiometry()
                stomath = product.getStoichiometryMath()
                if (stomath):
                    stonode = stomath.getMath()
                    # The stoichiometry could be set by complicated maths involving parameters
                    # elsewhere defined.
                    slist = MLtoList(stonode, self.__function_defs)
                    sto = MLfunc(slist, self.__globalParameters)
                # In versions 2, sto may not be specified
                if (not sto or sto!=sto): 
                    spec_ref = product.getId()
                    sto = self.__spec_refs[spec_ref]
                if (sto % 1 ) : 
                    raise NotImplementedError("Partial stoichiometry (%f) in reaction '%s'."%(sto, reac.getId()))
                sto = int(round(sto))
                for j in xrange(sto):
                    if(product.getSpecies() != "Empty"):
                        prods_f.append(product.getSpecies())
                        if (reversible): reacts_b.append(product.getSpecies())
            
            """
            # Modifiers REMOVED BECAUSE THEY DON'T SEEM TO DO ANYTHING
            mods = reac.getListOfModifiers()
            """
            
            # A CHECK AS TO WHETHER SPECIES ARE IN DIFFERENT COMPS OR NOT, and sort into reacs and surface-reacs
            
            # These are gonna be dictionaries 'compId':[list of species ids]
            reac_comp1 = {}
            reac_comp2 = {}
            reac_surf = {}

            reac_comp1_id = ''
            reac_comp2_id = ''
            reac_surf_id = ''
            
            # Store the important multiplying factor for the rates, the volume or 
            # area in which the reaction aoccurs.
            reactant_space_id = ''
            product_space_id = ''
            
            # This is going to be useful for math reactions.
            # A rate of a math reaction may be expressed in terms of some 
            # units other than the model units (if for example separate units 
            # are declared for parameters). The conversion will tell how that rate should be
            # converted to SBML units. For example some 'rate' may be micromol/second
            # and model substance units are mol wheras the math law assumes the rate is 
            # essentially mol/second  
            # One unsurmountable obstacle is that different species may have different 
            # factors (different units) but we can only assume they are the same at least within a reaction
            conv_factor_subs  = 0.0
            conv_factor_space = 0.0
            
            for r in reacts_f:
                if not conv_factor_subs: conv_factor_subs = self.__species[r][1]
                if not conv_factor_space : conv_factor_space
                # self.__species_comps may actually store a 'comp' or 'patch' string
                reac_comp = self.__species_comps[r]
                if reac_comp in self.__comps:
                    reactant_space_id = reac_comp
                    if not reac_comp1: 
                        if not conv_factor_space : conv_factor_space = self.__comps[reac_comp][1]
                        reac_comp1[reac_comp] = [r]
                        reac_comp1_id = reac_comp
                    elif reac_comp1.has_key(reac_comp): 
                        if not conv_factor_space : conv_factor_space = self.__comps[reac_comp][1]
                        reac_comp1[reac_comp].append(r)
                    else: raise NotImplementedError("Reaction '%s': reactants appear in different volumes."%reac.getId())
                elif reac_comp in self.__patches :
                    if not reactant_space_id: reactant_space_id = reac_comp
                    if not reac_surf: 
                        if not conv_factor_space : conv_factor_space = self.__patches[reac_comp][1]
                        reac_surf[reac_comp] = [r]
                        reac_surf_id = reac_comp
                    elif reac_surf.has_key(reac_comp): 
                        if not conv_factor_space : conv_factor_space = self.__patches[reac_comp][1]
                        reac_surf[reac_comp].append(r)
                    else: raise NotImplementedError("Reaction '%s': reactants and products appear in more than 1 surface."%reac.getId())
                else: 
                    raise NotImplementedError("Reaction '%s': reactants appears in unknown compartment."%reac.getId())
            for r in prods_f:
                if not conv_factor_subs: conv_factor_subs = self.__species[r][1]
                
                # self.__species_comps may actually store a 'comp' or 'patch' string
                reac_comp = self.__species_comps[r]
                if reac_comp in self.__comps:
                    product_space_id = reac_comp
                    if not reac_comp1: 
                        if not conv_factor_space : conv_factor_space = self.__comps[reac_comp][1]
                        reac_comp1[reac_comp] = [r]
                        reac_comp1_id = reac_comp
                    elif reac_comp1.has_key(reac_comp):
                        if not conv_factor_space : conv_factor_space = self.__comps[reac_comp][1]
                        reac_comp1[reac_comp].append(r)
                    elif not reac_comp2: 
                        if not conv_factor_space : conv_factor_space = self.__comps[reac_comp][1]
                        reac_comp2[reac_comp] = [r] 
                        reac_comp2_id = reac_comp
                    elif reac_comp2.has_key(reac_comp): 
                        if not conv_factor_space : conv_factor_space = self.__comps[reac_comp][1]
                        reac_comp2[reac_comp].append(r)
                    else: raise NotImplementedError("Reaction '%s': reactants and products appear in more than 2 volumes."%reac.getId())
                elif reac_comp in self.__patches :
                    if not product_space_id: product_space_id = reac_comp 
                    if not reac_surf: 
                        if not conv_factor_space : conv_factor_space = self.__patches[reac_comp][1]
                        reac_surf[reac_comp] = [r]
                        reac_surf_id = reac_comp
                    elif reac_surf.has_key(reac_comp): 
                        if not conv_factor_space : conv_factor_space = self.__patches[reac_comp][1]
                        reac_surf[reac_comp].append(r)
                    else: raise NotImplementedError("Reaction '%s': reactants and products appear in more than 1 surface."%reac.getId())
                else: 
                    raise NotImplementedError("Reaction '%s': reactants appears in unknown compartment."%reac.getId())
            
            if reac_comp1 and not (reac_comp2 or reac_surf): reac_type = 'volume reaction'
            elif reac_comp2 and not reac_surf: reac_type = 'volume volume reaction'
            elif reac_surf and not reac_comp1: reac_type = '2D surface reaction'
            elif reac_comp1 and reac_surf: reac_type = 'surface reaction'
            else: assert(false)
            
                                    
            if reac_type == 'volume reaction':
                r_f = Reaction()
                #reversible = reac.getReversible()
                if (reversible): 
                    r_f.setName(str(reac.getId())+"_f") 
                    r_f.setCompID(reac_comp1_id)
                    r_b = Reaction()
                    r_b.setName(str(reac.getId())+"_b")
                    r_b.setCompID(reac_comp1_id)
                else : 
                    r_f.setName(str(reac.getId())) 
                    r_f.setCompID(reac_comp1_id)
                
                # Tell the reaction which volsys it belongs to: compid+'_volsys'
                r_f.setVsys(reac_comp1_id+'_volsys')
                if (reversible): 
                    r_b.setVsys(reac_comp1_id+'_volsys')
                
                # Tricky situation where the species is a boundary spec, but not constant. This means the 
                # quantity can change in rules and assignments but NOT reactions. No direct support for 
                # that in STEPS really (a spec is either clamped or not) but can do a little trick here 
                # where I equate the rhs stoichiometry to whatever's on the lhs.
                # Extra difficulty is that the reaction may be reversible - treat separately.
                
                # First create a temporary structure to map lhs spec id -> stoichiometry
                lhs_stoich = {}
                for l in reacts_f:
                    if l in lhs_stoich:
                        lhs_stoich[l] += 1
                    else:
                        lhs_stoich[l] = 1   
                
                # Now do the same to map rhs spec id -> stoichiometry
                rhs_stoich = {}
                for r in prods_f:
                    if r in rhs_stoich:
                        rhs_stoich[r] += 1
                    else:
                        rhs_stoich[r] = 1  
                
                for lsid in lhs_stoich: 
                    if (self.__species_const_bc[lsid][1] == True):
                        # First find the number on the lhs
                        n_lhs = lhs_stoich[lsid]
                        # And rhs
                        if lsid in rhs_stoich: n_rhs = rhs_stoich[lsid]
                        else: n_rhs = 0
                        # Now calculate the difference
                        lhs_diff_rhs = int(n_lhs - n_rhs)
                        
                        # If we have more on the lhs than the right we need to add to the rhs
                        if (lhs_diff_rhs > 0):
                            while (lhs_diff_rhs > 0):
                                prods_f.append(lsid)
                                lhs_diff_rhs -= 1
                        elif (lhs_diff_rhs < 0):
                            while (lhs_diff_rhs < 0):
                                prods_f.remove(lsid)
                                lhs_diff_rhs += 1
                
                # Now recalculate the stoichiometries to check the rhs
                lhs_stoich = {}
                for l in reacts_f:
                    if l in lhs_stoich:
                        lhs_stoich[l] += 1
                    else:
                        lhs_stoich[l] = 1 
                
                rhs_stoich = {}
                for r in prods_f:
                    if r in rhs_stoich:
                        rhs_stoich[r] += 1
                    else:
                        rhs_stoich[r] = 1                            
                                
                for rsid in rhs_stoich:
                    if (self.__species_const_bc[rsid][1] == True):
                        n_rhs = rhs_stoich[rsid]
                        if (rsid in lhs_stoich): n_lhs = lhs_stoich[rsid]
                        else: n_lhs = 0
                        rhs_diff_lhs = int(n_rhs-n_lhs)
                        
                        # If we have more on the rhs than the left we need to remove from the rhs
                        if (rhs_diff_lhs > 0):
                            while (rhs_diff_lhs > 0):
                                prods_f.remove(rsid)
                                rhs_diff_lhs -= 1
                        # If less on the rhs we have to add some
                        elif (rhs_diff_lhs < 0):
                            while (rhs_diff_lhs < 0):
                                prods_f.append(rsid)
                                rhs_diff_lhs += 1
                
                if (reversible):
                    # First create a temporary structure to map lhs spec id -> stoichiometry
                    lhs_stoich = {}
                    for l in reacts_b:
                        if l in lhs_stoich:
                            lhs_stoich[l] += 1
                        else:
                            lhs_stoich[l] = 1   
                    
                    # Now do the same to map rhs spec id -> stoichiometry
                    rhs_stoich = {}
                    for r in prods_b:
                        if r in rhs_stoich:
                            rhs_stoich[r] += 1
                        else:
                            rhs_stoich[r] = 1  
                    
                    for lsid in lhs_stoich: 
                        if (self.__species_const_bc[lsid][1] == True):
                            # First find the number on the lhs
                            n_lhs = lhs_stoich[lsid]
                            # And rhs
                            if lsid in rhs_stoich: n_rhs = rhs_stoich[lsid]
                            else: n_rhs = 0
                            # Now calculate the difference
                            lhs_diff_rhs = int(n_lhs - n_rhs)
                            
                            # If we have more on the lhs than the right we need to add to the rhs
                            if (lhs_diff_rhs > 0):
                                while (lhs_diff_rhs > 0):
                                    prods_b.append(lsid)
                                    lhs_diff_rhs -= 1
                            elif (lhs_diff_rhs < 0):
                                while (lhs_diff_rhs < 0):
                                    prods_b.remove(lsid)
                                    lhs_diff_rhs += 1
                    
                    # Now recalculate the stoichiometries to check the rhs
                    lhs_stoich = {}
                    for l in reacts_b:
                        if l in lhs_stoich:
                            lhs_stoich[l] += 1
                        else:
                            lhs_stoich[l] = 1 
                    
                    rhs_stoich = {}
                    for r in prods_b:
                        if r in rhs_stoich:
                            rhs_stoich[r] += 1
                        else:
                            rhs_stoich[r] = 1                            
                                    
                    for rsid in rhs_stoich:
                        if (self.__species_const_bc[rsid][1] == True):
                            n_rhs = rhs_stoich[rsid]
                            if (rsid in lhs_stoich): n_lhs = lhs_stoich[rsid]
                            else: n_lhs = 0
                            rhs_diff_lhs = int(n_rhs-n_lhs)
                            
                            # If we have more on the rhs than the left we need to remove from the rhs
                            if (rhs_diff_lhs > 0):
                                while (rhs_diff_lhs > 0):
                                    prods_b.remove(rsid)
                                    rhs_diff_lhs -= 1
                            # If less on the rhs we have to add some
                            elif (rhs_diff_lhs < 0):
                                while (rhs_diff_lhs < 0):
                                    prods_b.append(rsid)
                                    rhs_diff_lhs += 1                
                
                lhsList_f = []
                rhsList_f = []
                for rct_f in reacts_f:
                    lhsList_f.append(self.__mdl.getSpec(rct_f))
                for pro_f in prods_f:
                    rhsList_f.append(self.__mdl.getSpec(pro_f))
                  
                if (reversible):
                    lhsList_b = []
                    rhsList_b = []
                    for rct_b in reacts_b:
                        lhsList_b.append(self.__mdl.getSpec(rct_b))
                    for pro_b in prods_b:
                        rhsList_b.append(self.__mdl.getSpec(pro_b))
                
                r_f.setReacts(reacts_f)
                r_f.setLhs(lhsList_f)
                order_f = len(lhsList_f)
                
                if (reversible): 
                    r_b.setProds(prods_b)
                    r_b.setRhs(rhsList_b)
                                        
                r_f.setProds(prods_f)
                r_f.setRhs(rhsList_f)
                if (reversible):
                    r_b.setReacts(reacts_b)
                    r_b.setLhs(lhsList_b)
                    order_b = len(lhsList_b)
                
                #######################################
                
                # Kinetic constants
                params = {}
                kLaw = reac.getKineticLaw()
                if not kLaw: 
                    print "WARNING: Reaction '%s' has undefined kinetic law " \
                    "and will be ignored."%reac.getId()
                    continue
                parameters = kLaw.getListOfParameters()
                
                # Attempt to match the parameters to the right reaction (if reversible, 
                # i.e. forward or back - very difficult from SBML). This 
                # involves looking at the parameters' units and trying to match
                # to the reaction by order, if orders are different.
                # Even if reaction is not reversible, a number of parameters may exist
                # possibly representing things we don't support right now, current,
                # temperature etc 
                
                # Start with a list of params by [id, order, value(STEPS units)]
                params = []
                # Store a list of bad parameters for error checking
                bad_params = []
                
                for p in parameters:
                    p_id = str(p.getId())
                    p_value = p.getValue()
                    p_units = p.getUnits()
                    if (not p_value): 
                        if p_id in self.__globalParameters:
                            p_value, p_factor  = self.__globalParameters[p_id]     
                            p_order = self.__glob_params_order[p_id]
                            if (p_order == 'not reac'): 
                                if (strict): raise NotImplementedError("Reaction '%s' parameter '%s' not the right units for a reaction parameter."%(reac.getId()), p_id)
                            params.append([p_id, p_order, p_value, p_factor])
                            # PARAMS ARE IN MODEL UNITS!
                            # Will ckeck later if this param is not a reaction param
                            continue  
                    else: pass
                    
                    # For local parameters we are looking specifically for things involved in a reaction. If 
                    # they are other things (checked in the units if specified) we simply don't add them.
                    if (p_units): 
                        
                        # Using a for loop that can be exited with a break, 
                        # for want of a better alternative
                        for j in [0]:
                            badunit = False
                        
                            unitdef = self.__model.getUnitDefinition(p_units) 
                            if not (unitdef):
                                if strict: 
                                    print "WARNING: Local parameter '%s' has unknown unit definition " \
                                    "and model units will be assumed."%p_id
                                    break
                                else:
                                    # This seems to happen with dimensionless units now
                                    badunit = True
                                    bad_params.append([p_id, -1, p_value, p_factor])
                                    break
                            units = unitdef.getListOfUnits()  
                            
                            # Units for a reaction should be (conc)^-(n-1).(time)^-1
                            # Where n gives the order of the reaction
                            # NOTE: If a unit is bad it is simply not added to the list of parameters for the search.
                            
                            # PARAM MUST NOW BE KEPT IN SBML UNITS
                            gotTime = False
                            gotLitre = False
                            gotMetre = False
                            gotMole = False
                            gotSecond = False
                            order = -1
                            p_factor = 1
                            for u in units:
                                mult = u.getMultiplier()
                                sca = u.getScale()
                                exp = u.getExponent()
                                # pfactor will convert a number in base units to SBML
                                p_factor *= unitsDef_to_STEPS(mult, sca, exp)
                                p_value *= unitsDef_to_STEPS(mult, sca, exp)
                                if (u.isDimensionless() or u.isItem()):
                                    # OK, do nothing
                                    # Don't think we'll get here with new libsbml, bug?
                                    continue
                                elif (u.isLitre()):
                                    if (gotLitre or gotMetre): badunit = True
                                    elif (exp%1 != 0): badunit = True
                                    elif (order != exp+1) and (order != -1): badunit = True 
                                    order = exp+1
                                    if (order < 0): badunit = True
                                    # Convert to meters^3
                                    p_factor *= math.pow(0.001, exp)
                                    p_value *= math.pow(0.001, exp)
                                    
                                    p_factor/=math.pow(self.__volume_units, exp)
                                    p_value/=math.pow(self.__volume_units, exp)

                                    gotLitre = True
                                elif (u.isMetre()): 
                                    if (gotLitre or gotMetre): badunit = True
                                    elif (exp%3 != 0): badunit = True
                                    elif (order != (exp/3)+1) and (order != -1): badunit = True
                                    order  = (exp/3)+1
                                    if (order < 0): badunit = True
                                    gotMetre = True
                                    
                                    p_factor/=math.pow(self.__length_units, exp)
                                    p_value/=math.pow(self.__length_units, exp)
                                    
                                elif (u.isMole()):
                                    if (gotMole): badunit = True
                                    elif (exp%1 != 0): badunit  = True
                                    elif (order != -(exp-1)) and (order != -1): badunit = True
                                    order = -(exp-1)
                                    if (order < 0): badunit = True
                                    gotMole = True
                                    
                                    p_factor/=math.pow(self.__substance_units, exp)
                                    p_value/=math.pow(self.__substance_units, exp) 
                                                                       
                                elif (u.isSecond()):
                                    if (gotSecond): badunit = True
                                    elif (exp != -1):  badunit = True
                                    gotSecond = True
                                    
                                    p_factor/=math.pow(self.__time_units, exp)
                                    p_value/=math.pow(self.__time_units, exp) 
                                                                        
                                else: 
                                    raise NotImplementedError("Reaction '%s' parameter '%s' contains unsupported unit."%(reac.getId(), p_id))
                                    badunit = True
                                # Breaking here means the local parameter is not a reaction parameter and will not be added to params
                                if(badunit): 
                                    if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                    "parameter units and will cause error if in rate expression."%p_id
                                    bad_params.append([p_id, order, p_value, p_factor])    
                                    break
                            # Set first order if so
                            if (gotSecond and order == -1): order = 1    
                            # Some sanity checks: break if true and don't add param to params (print warning??)
                            if (badunit or order < 0 or gotSecond == False or p_value < 0): 
                                if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                "parameter units and will cause error if in rate expression."%p_id
                                bad_params.append([p_id, order, p_value, p_factor])
                                break
                            if order != 1: 
                                if not (gotLitre or gotMetre): 
                                    if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                    "parameter units and will cause error if in rate expression."%p_id
                                    bad_params.append([p_id, order, p_value, p_factor])                            
                                    break
                            # We may still have say (for second order reaction): m^3/s but substance units are molar. Need to convert in this case:
                            if (gotMole): params.append([p_id, order, p_value, p_factor])
                            else:
                                # First convert to moles:
                                p_value*=math.pow(AVOGADRO, order-1)
                                p_factor*=math.pow(AVOGADRO, order-1)
                                # And into SBML units, remmeber here we don't have the (usually negative) exp:
                                # NOTE: If we had molar units we already converted with a (usually) negative exp
                                p_factor*=math.pow(self.__substance_units, order-1)
                                p_value*=math.pow(self.__substance_units, order-1)                             
                                params.append([p_id, order, p_value, p_factor])
                            break
                        
                    else: 
                        # No units, assume default, which are MODEL UNITS OF VOLUME, SUBSTANCE, TIME 
                        order = 'unknown'
                        # Also we don't know the STEPS conversion (meters, moles, seconds) because we 
                        # don't know the order:
                        params.append([p_id, order, p_value, 1.0]) 
                
                kMath = kLaw.getMath()
                # Attempt to check the form of the rate in kMath to what is expected
                # Parameters are kind of a free thing - we don't know which parameter is 
                # which (if it is reversible)
                # First provide a new function to read the list and return parameters (from the 
                # global or local parameters) involved with a return value. Then construct two lists-
                # one for each combination of parameters (for a reversible reaction) of the 
                # form that is expected. If the return from one of the lists is the same as the original return 
                # we have a good reaction. 
                # One complication, as usual, is that species quantities in the expression may be 
                # amount or concentration. The concentration form is quite straight forward, but
                # if amounts the compartment volume will take a power depending on the order.
                rate_list = MLtoList(kMath, self.__function_defs)
                
                # Supply a dictionary of math parameters, but give random values to work out rates. 
                # This is because any parameter could be initialized with zero giving zero rates which would compare as equal
                # Temporary structures must be used otherwise global values will get overwritten
                
                all_variables = {}
                params_variables = {}
                
                # First the local parameters 
                paramsloc = {}
                paramsloc_temp = {}
                for p in params:
                    paramsloc[p[0]] = [p[2], p[3]]
                    paramsloc_temp[p[0]] = [random.random(), p[3]]                
                # Note: paramsloc_temp should come after global parameters because it may be a shadower
                
                params_temp = {}
                for p in self.__globalParameters:
                    params_temp[p] = [random.random(), self.__globalParameters[p][1] ]
                    
                species_temp = {}
                for s in self.__species:
                    species_temp[s] = [random.random(), self.__species[s][1] ]
                
                comps_temp = {}
                for c in self.__comps:
                    comps_temp[c] = [random.random(), self.__comps[c][1] ]
                
                all_variables.update(params_temp)
                all_variables.update(paramsloc_temp)
                all_variables.update(species_temp)
                all_variables.update(comps_temp)
                
                params_variables.update(self.__globalParameters)
                params_variables.update(paramsloc_temp)
                
                
                mathreac = False
                params_id = []
                
                if strict:
                    rate = rate_func(rate_list, all_variables, params_variables, params_id, bad_params)
                    # Separate out the parameters
                    if not reversible:
                        if len(params_id) != 1:
                            if (strict): raise NotImplementedError("Reaction rate maths in reaction '%s' is not of supported form."%reac.getId())
                    else:
                        if len(params_id) != 2:
                            raise NotImplementedError("Reaction rate maths in reaction '%s' is not of supported form."%reac.getId())
                else:
                    try: rate = rate_func(rate_list, all_variables, params_variables, params_id, bad_params)
                    except: 
                        mathreac = True
                    if not mathreac:
                        if not reversible:
                            if len(params_id) != 1: mathreac = True
                        else:
                            if len(params_id) != 2: mathreac = True
                
                # Update the params list
                for p_id in params_id:
                    if (p_id in bad_params):
                        if strict: raise NotImplementedError("Reaction rate maths in reaction '%s' includes" \
                            "parameter '%s' with unexpected units."%(reac.getId(), p_id))
                    if (p_id in self.__globalParameters):
                        gotP = False
                        for p in params:
                            if (p[0] == p_id): gotP = True 
                        if not gotP:
                            p_value, p_factor  = self.__globalParameters[p_id]    
                            p_order = self.__glob_params_order[p_id] 
                            if (p_order == 'not reac'): 
                                if strict: raise NotImplementedError("Reaction '%s' parameter '%s' not the right units for a reaction parameter."%(reac.getId()), p_id)
                            elif (type(p_order) in (int, float)):
                                pass
                            params.append([p_id, p_order, p_value, p_factor])
                
                # OK, now we are here we may have some parameters that have not specified units, the flag is p_order = 'unknown'
                # otherwise they have been converted already
                
                # Right, hopefully that worked (error will have been thrown if not). 
                # Now the tricky part of constructing the expected list
                
                # Lets treat forward and backward separately (for reversible reactions), but create a holder
                # for a reverse reaction regardless (initialised to 0, ignored for non-reversible reactions)
                
                # We may already know that the raection has bad parameters etc and is a 'math reac'
                if not mathreac:
                    correct_rate_list_1 = ['-', [None, 0.0]]
                    
                    if (reversible): correct_rate_list_2 = ['-', [None, None]]
                    
                    # just a reminder reacts_f is a list of the reactant species by string id for the forward reaction
                    # and reacts_b is a list of reactant species for the backward reaction (if there is one)
                    # reac_comp is the compartment id
                    
                    # The 'forward' reaction first: the compartment, parameter0, and the species
                    for_vars_1 = [reac_comp1_id, params_id[0]]+reacts_f
                    correct_rate_list_1[1][0] = make_rate_list(for_vars_1, self.__species_subst_units, reac_comp1_id)            
                    
                    if (reversible):
                        back_vars_1 = [reac_comp1_id, params_id[1]]+reacts_b
                        correct_rate_list_1[1][1] = make_rate_list(back_vars_1, self.__species_subst_units, reac_comp1_id)
                    
                    if (reversible):
                        for_vars_2 = [reac_comp1_id, params_id[1]]+reacts_f
                        correct_rate_list_2[1][0] = make_rate_list(for_vars_2, self.__species_subst_units, reac_comp1_id)
                        back_vars_2 = [reac_comp1_id, params_id[0]]+reacts_b
                        correct_rate_list_2[1][1] = make_rate_list(back_vars_2, self.__species_subst_units, reac_comp1_id)
                    correct_rate_1 = MLfunc(correct_rate_list_1, all_variables)
                    if (reversible): 
                        correct_rate_2 = MLfunc(correct_rate_list_2, all_variables)
                    
                    setparam_f = False
                    if (reversible): setparam_b = False
                    if _float_approx_equal(rate, correct_rate_1):
                        # We have found the right parameters. Set them to the correct reaction
                        # We've got the ids, but need to find the right parameter in the params list
                        for p in params:
                            if not setparam_f and p[0] == params_id[0]:
                                if (p[1] == 'not reac'):
                                    if strict: raise NotImplementedError("Reaction rate parameter for reaction '%s' has unexpected units."%reac.getId())
                                    else: 
                                        mathreac = True
                                        break
                                else:
                                    order = r_f.getOrder()
                                    # CHANGED: OK, now values are stored in SBML units. We should be able to convert just by multipliying by the factor
                                    # and converting the volume to litres, BUT we might not know the factor if no units (so can't use p[3])
                                    
                                    kvalue = p[2]
                                    
                                    factor = (math.pow(self.__volume_units*1000.0, order-1))/(self.__time_units*math.pow(self.__substance_units, order-1))
                                    
                                    # Not sure what to do about the mole - there seems to be no global substance units
                                    p[2] = kvalue*factor
                                    p[3] = factor
                                    
                                r_f.setKName(p[0])
                                r_f.setKValue(p[2])
                                #r_f.setUnitsConv(p[3])  
                                reactions.append(r_f)
                                setparam_f  = True
                            elif (reversible):
                                if (not setparam_b and p[0] == params_id[1]):
                                    if (p[1] == 'not reac'): # or p[1] != len(r_b.getLhs())):
                                        if strict: raise NotImplementedError("Reaction rate parameter for reversible reaction '%s' has unexpected units."%reac.getId())
                                        else:
                                            mathreac = True
                                            break
                                    else:
                                        order = r_b.getOrder()
                                        kvalue = p[2]
                                        
                                        factor = (math.pow(self.__volume_units*1000.0, order-1))/(self.__time_units*math.pow(self.__substance_units, order-1))
                                        
                                        p[2] = kvalue*factor
                                        p[3] = factor    
                                    
                                    r_b.setKName(p[0])
                                    r_b.setKValue(p[2])
                                    #r_b.setUnitsConv(p[3])
                                    reactions.append(r_b)
                                    setparam_b = True
                    elif (reversible):
                        if _float_approx_equal(rate, correct_rate_2):
                            # Yatta, we have found the right params
                            for p in params:
                                if not setparam_f and p[0] == params_id[1]:
                                    if (p[1] == 'not reac'): # or p[1] != len(r_f.getLhs())):
                                        if strict: raise NotImplementedError("Reaction rate parameter for reaction '%s' has unexpected units."%reac.getId())
                                        else:
                                            mathreac = True
                                            break
                                    # 'unknown' basically means that the units were not specified, so we have to use model units
                                    else:
                                        order = r_f.getOrder()
                                        kvalue = p[2]
                                        
                                        factor = (math.pow(self.__volume_units*1000.0, order-1))/(self.__time_units*math.pow(self.__substance_units, order-1))
                                        p[2] = kvalue*factor
                                        p[3] = factor    
                                        
                                    r_f.setKName(p[0])
                                    r_f.setKValue(p[2])
                                    #r_f.setUnitsConv(p[3])  
                                    reactions.append(r_f)
                                    setparam_f  = True
                                elif (not setparam_b and p[0] == params_id[0]):
                                    if (p[1] == 'not reac'):# or p[1] != len(r_b.getLhs())):
                                        if strict: raise NotImplementedError("Reaction rate parameter order for reversible reaction '%s' has unexpected units."%reac.getId())
                                        else:
                                            mathreac = True
                                            break
                                    # 'unknown' basically means that the units were not specified, so we have to use model units
                                    else:
                                        order = r_b.getOrder()
                                        kvalue = p[2]
                                        
                                        factor = (math.pow(self.__volume_units*1000.0, order-1))/(self.__time_units*math.pow(self.__substance_units, order-1))
                                        #kvalue*= p[3]*math.pow(1000.0, order-1)
                                        
                                        p[2] = kvalue*factor 
                                        p[3] = factor

                                    r_b.setKName(p[0])
                                    r_b.setKValue(p[2])
                                    #r_b.setUnitsConv(p[3])
                                    reactions.append(r_b)
                                    setparam_b = True   
                        else:
                            if strict: raise NotImplementedError("Reaction rate maths in reaction %s is not of expected rate."%reac.getId())
                            else: mathreac = True
                    else:
                        if strict: raise NotImplementedError("Reaction rate maths in reaction %s is not of expected rate."%reac.getId())
                        else: mathreac = True
                
                # Separate treatment for "math reacs"
                if mathreac:
                    # Assuming that it makes no difference whether the reaction is reversilbe or not if the maths is strange
                    r_f.setKName('unknown')
                    r_f.setKValue(0.0)
                    
                    # This "correct rate" will be just the spcies and the comp- no parameter.
                    # The changing parameter will then be the actual maths divided by this rate list
                    for_vars_1 = [reac_comp1_id]+reacts_f
                                        
                    # Add the "bad parameters"                
                    for p in bad_params:
                        paramsloc[p[0]] = [p[2], p[3]]
                    
                    # Save the factor for conversion to STEPS units
                    ord_tmp = r_f.getOrder()
                    # Factor is going to convert the returned rate in base units to SBML units
                    factor = 1.0*math.pow(conv_factor_space, ord_tmp-1)
                    factor /= math.pow(conv_factor_subs, ord_tmp-1)                
                    
                    math_reactions.append([r_f, rate_list, make_rate_list(for_vars_1, self.__species_subst_units, reac_comp1_id), paramsloc, factor])
                    reactions.append(r_f)
            else:
                ############### A SURFACE REACTION ##################
                r_f = SurfaceReaction()
                if (reversible): 
                    r_f.setName(str(reac.getId())+"_f") 
                    r_b = SurfaceReaction()
                else : r_f.setName(str(reac.getId())) 
                if (reversible): r_b.setName(str(reac.getId())+"_b")
                
                
                # Tricky situation where the species is a boundary spec, but not constant. This means the 
                # quantity can change in rules and assignments but NOT reactions. No direct support for 
                # that in STEPS really (a spec is either clamped or not) but can do a little trick here 
                # where I equate the rhs stoichiometry to whatever's on the lhs.
                # Extra difficulty is that the reaction may be reversible - treat separately.
                
                # First create a temporary structure to map lhs spec id -> stoichiometry
                lhs_stoich = {}
                for l in reacts_f:
                    if l in lhs_stoich:
                        lhs_stoich[l] += 1
                    else:
                        lhs_stoich[l] = 1   
                
                # Now do the same to map rhs spec id -> stoichiometry
                rhs_stoich = {}
                for r in prods_f:
                    if r in rhs_stoich:
                        rhs_stoich[r] += 1
                    else:
                        rhs_stoich[r] = 1  
                                
                for lsid in lhs_stoich: 
                    if (self.__species_const_bc[lsid][1] == True):
                        # First find the number on the lhs
                        n_lhs = lhs_stoich[lsid]
                        # And rhs
                        if lsid in rhs_stoich: n_rhs = rhs_stoich[lsid]
                        else: n_rhs = 0
                        # Now calculate the difference
                        lhs_diff_rhs = int(n_lhs - n_rhs)
                        
                        # If we have more on the lhs than the right we need to add to the rhs
                        if (lhs_diff_rhs > 0):
                            while (lhs_diff_rhs > 0):
                                #rhsList.append(self.__mdl.getSpec(id))
                                prods_f.append(lsid)
                                lhs_diff_rhs -= 1
                        elif (lhs_diff_rhs < 0):
                            while (lhs_diff_rhs < 0):
                                #rhsList.remove(self.__mdl.getSpec(id))
                                prods_f.remove(lsid)
                                lhs_diff_rhs += 1
                
                # Now recalculate the stoichiometries to check the rhs
                lhs_stoich = {}
                for l in reacts_f:
                    if l in lhs_stoich:
                        lhs_stoich[l] += 1
                    else:
                        lhs_stoich[l] = 1 
                
                rhs_stoich = {}
                for r in prods_f:
                    if r in rhs_stoich:
                        rhs_stoich[r] += 1
                    else:
                        rhs_stoich[r] = 1                            
                                
                for rsid in rhs_stoich:
                    if (self.__species_const_bc[rsid][1] == True):
                        n_rhs = rhs_stoich[rsid]
                        if (rsid in lhs_stoich): n_lhs = lhs_stoich[rsid]
                        else: n_lhs = 0
                        rhs_diff_lhs = int(n_rhs-n_lhs)
                        
                        # If we have more on the rhs than the left we need to remove from the rhs
                        if (rhs_diff_lhs > 0):
                            while (rhs_diff_lhs > 0):
                                #rhsList.remove(self.__mdl.getSpec(id))
                                prods_f.remove(rsid)
                                rhs_diff_lhs -= 1
                        # If less on the rhs we have to add some
                        elif (rhs_diff_lhs < 0):
                            while (rhs_diff_lhs < 0):
                                #rhsList.append(self.__mdl.getSpec(id))
                                prods_f.append(rsid)
                                rhs_diff_lhs += 1
                                
                if (reversible):
                    # First create a temporary structure to map lhs spec id -> stoichiometry
                    lhs_stoich = {}
                    for l in reacts_b:
                        if l in lhs_stoich:
                            lhs_stoich[l] += 1
                        else:
                            lhs_stoich[l] = 1   
                    
                    # Now do the same to map rhs spec id -> stoichiometry
                    rhs_stoich = {}
                    for r in prods_b:
                        if r in rhs_stoich:
                            rhs_stoich[r] += 1
                        else:
                            rhs_stoich[r] = 1  
                    
                    for lsid in lhs_stoich: 
                        if (self.__species_const_bc[lsid][1] == True):
                            # First find the number on the lhs
                            n_lhs = lhs_stoich[lsid]
                            # And rhs
                            if lsid in rhs_stoich: n_rhs = rhs_stoich[lsid]
                            else: n_rhs = 0
                            # Now calculate the difference
                            lhs_diff_rhs = int(n_lhs - n_rhs)
                            
                            # If we have more on the lhs than the right we need to add to the rhs
                            if (lhs_diff_rhs > 0):
                                while (lhs_diff_rhs > 0):
                                    #rhsList.append(self.__mdl.getSpec(id))
                                    prods_b.append(lsid)
                                    lhs_diff_rhs -= 1
                            elif (lhs_diff_rhs < 0):
                                while (lhs_diff_rhs < 0):
                                    #rhsList.remove(self.__mdl.getSpec(id))
                                    prods_b.remove(lsid)
                                    lhs_diff_rhs += 1
                    
                    # Now recalculate the stoichiometries to check the rhs
                    lhs_stoich = {}
                    for l in reacts_b:
                        if l in lhs_stoich:
                            lhs_stoich[l] += 1
                        else:
                            lhs_stoich[l] = 1 
                    
                    rhs_stoich = {}
                    for r in prods_b:
                        if r in rhs_stoich:
                            rhs_stoich[r] += 1
                        else:
                            rhs_stoich[r] = 1                            
                                    
                    for rsid in rhs_stoich:
                        if (self.__species_const_bc[rsid][1] == True):
                            n_rhs = rhs_stoich[rsid]
                            if (rsid in lhs_stoich): n_lhs = lhs_stoich[rsid]
                            else: n_lhs = 0
                            rhs_diff_lhs = int(n_rhs-n_lhs)
                            
                            # If we have more on the rhs than the left we need to remove from the rhs
                            if (rhs_diff_lhs > 0):
                                while (rhs_diff_lhs > 0):
                                    #rhsList.remove(self.__mdl.getSpec(id))
                                    prods_b.remove(rsid)
                                    rhs_diff_lhs -= 1
                            # If less on the rhs we have to add some
                            elif (rhs_diff_lhs < 0):
                                while (rhs_diff_lhs < 0):
                                    #rhsList.append(self.__mdl.getSpec(id))
                                    prods_b.append(rsid)
                                    rhs_diff_lhs += 1                        
                
                # Lets call comp1 the inner comp and comp2 the outer comp. Impossible to check against other
                # surface reactions and all possible ways to set which comp is 'inner' and which 'outer' to
                # a comp, so that will have to be done later
                
                olhsList_f = []
                ilhsList_f = []
                slhsList_f = []
                
                orhsList_f = []
                irhsList_f = []
                srhsList_f = []
                
                for rct_f in reacts_f:
                    if (reac_comp1.values() and rct_f in reac_comp1.values()[0]): ilhsList_f.append(self.__mdl.getSpec(rct_f))
                    elif (reac_comp2.values() and rct_f in reac_comp2.values()[0]): olhsList_f.append(self.__mdl.getSpec(rct_f))
                    elif (reac_surf.values() and rct_f in reac_surf.values()[0]): slhsList_f.append(self.__mdl.getSpec(rct_f))
                    else: assert(False)
                for pro_f in prods_f:
                    if (reac_comp1.values() and pro_f in reac_comp1.values()[0]): irhsList_f.append(self.__mdl.getSpec(pro_f))
                    elif (reac_comp2.values() and pro_f in reac_comp2.values()[0]): orhsList_f.append(self.__mdl.getSpec(pro_f))
                    elif (reac_surf.values() and pro_f in reac_surf.values()[0]): srhsList_f.append(self.__mdl.getSpec(pro_f))
                    else: assert(False)
                
                if  reacts_f == prods_f: 
                    print "WARNING: Reaction '%s' will do nothing and will be ignored (does it involve all boundary species?)."%(reac.getId()) 
                    continue
                
                if (reversible):
                    olhsList_b = []
                    ilhsList_b = []
                    slhsList_b = []
                    
                    orhsList_b = []
                    irhsList_b = []
                    srhsList_b = []

                    for rct_b in reacts_b:
                        if (reac_comp1.values() and rct_b in reac_comp1.values()[0]): ilhsList_b.append(self.__mdl.getSpec(rct_b))
                        elif (reac_comp2.values() and rct_b in reac_comp2.values()[0]): olhsList_b.append(self.__mdl.getSpec(rct_b))
                        elif (reac_surf.values() and rct_b in reac_surf.values()[0]): slhsList_b.append(self.__mdl.getSpec(rct_b))
                        else: assert(False)
                    for pro_b in prods_b:
                        if (reac_comp1.values() and pro_b in reac_comp1.values()[0]): irhsList_b.append(self.__mdl.getSpec(pro_b))
                        elif (reac_comp2.values() and pro_b in reac_comp2.values()[0]): orhsList_b.append(self.__mdl.getSpec(pro_b))
                        elif (reac_surf.values() and pro_b in reac_surf.values()[0]): srhsList_b.append(self.__mdl.getSpec(pro_b))
                        else: assert(False)
                
                r_f.setReacts(reacts_f)
                r_f.setILhs(ilhsList_f)
                r_f.setOLhs(olhsList_f)
                r_f.setSLhs(slhsList_f)
                
                order_f = len(ilhsList_f)+len(olhsList_f)+len(slhsList_f)
                
                if (reversible): 
                    r_b.setProds(prods_b)
                    r_b.setORhs(orhsList_b)
                    r_b.setSRhs(srhsList_b)
                    r_b.setIRhs(irhsList_b)
                                
                r_f.setProds(prods_f)
                r_f.setORhs(orhsList_f)
                r_f.setIRhs(irhsList_f)
                r_f.setSRhs(srhsList_f)
                if (reversible):
                    r_b.setReacts(reacts_b)
                    r_b.setILhs(ilhsList_b)
                    r_b.setOLhs(olhsList_b)
                    r_b.setSLhs(slhsList_b)
                    order_b = len(ilhsList_b)+len(olhsList_b)+len(slhsList_b)
                
                #######################################
                
                # Kinetic constants
                params = {}
                kLaw = reac.getKineticLaw()
                if not kLaw: 
                    print "WARNING: Reaction '%s' has undefined kinetic law " \
                    "and will be ignored."%reac.getId()
                    continue
                parameters = kLaw.getListOfParameters()
                
                # Attempt to match the parameters to the right reaction (if reversible, 
                # i.e. forward or back - very difficult from SBML). This 
                # involves looking at the parameters' units and trying to match
                # to the reaction by order, if orders are different.
                # Even if reaction is not reversible, a number of parameters may exist
                # possibly representing things we don't support right now, current,
                # temperature etc 
                
                # Start with a list of params by [id, order, value(STEPS units)]
                params = []
                # Store a list of bad parameters for error checking
                bad_params = []
                
                for p in parameters:
                    p_id = str(p.getId())
                    p_value = p.getValue()
                    p_units = p.getUnits()
                    if (not p_value): 
                        if p_id in self.__globalParameters:
                            p_value, p_factor  = self.__globalParameters[p_id]     
                            p_order = self.__glob_params_order[p_id]
                            if (p_order == 'not reac'): 
                                if (strict): raise NotImplementedError("Reaction '%s' parameter '%s' not the right units for a reaction parameter."%(reac.getId()), p_id)

                            params.append([p_id, p_order, p_value, p_factor])
                            # Reaction params already converted to STEPS, Molar units.   --- NO. Units are NOT STEPS units yet!
                            # Will ckeck later if this param is not a reaction param
                            continue  
                    else:
                        pass
                    # For local parameters we are looking specifically for things involved in a reaction. If 
                    # they are other things (checked in the units if specified) we simply don't add them.
                    if (p_units): 
                        # Using a for loop that can be exited with a break, 
                        # for want of a better alternative
                        for j in [0]:
                            badunit = False
                        
                            unitdef = self.__model.getUnitDefinition(p_units) 
                            if not (unitdef):
                                if strict: 
                                    print "WARNING: Local parameter '%s' has unknown unit definition " \
                                "and will cause error if used in rate expression."%p_id
                                    break
                                else:
                                    # This seems to happen with dimensionless units now
                                    badunit = True
                                    bad_params.append([p_id, -1, p_value, p_factor])
                                    break
                            units = unitdef.getListOfUnits()  
                            
                            # Units for a reaction should be (conc)^-(n-1).(time)^-1
                            # Where n gives the order of the reaction
                            # NOTE: If a unit is bad it is imply not added to the list of parameters for the search.
                            if not reac_type == '2D surface reaction':
                                gotTime = False
                                gotLitre = False
                                gotMetre = False
                                gotMole = False
                                gotSecond = False
                                order = -1
                                p_factor = 1
                                for u in units:
                                    mult = u.getMultiplier()
                                    sca = u.getScale()
                                    exp = u.getExponent()
                                    p_factor *= unitsDef_to_STEPS(mult, sca, exp)
                                    p_value *= unitsDef_to_STEPS(mult, sca, exp)
                                    if (u.isDimensionless() or u.isItem()):
                                        # OK, do nothing
                                        # Don't think we'll get here with new libsbml, bug?
                                        continue
                                    elif (u.isLitre()):
                                        if (gotLitre or gotMetre): badunit = True
                                        elif (exp%1 != 0): badunit = True
                                        elif (order != exp+1) and (order != -1): badunit = True # If we already have an order from substance, break if different
                                        order = exp+1
                                        if (order < 0): badunit = True
                                        # Convert to meters^3
                                        p_factor *= math.pow(0.001, exp)
                                        p_value *= math.pow(0.001, exp)
                                        
                                        p_factor/=math.pow(self.__volume_units, exp)
                                        p_value/=math.pow(self.__volume_units, exp)
                                        
                                        gotLitre = True
                                    elif (u.isMetre()): 
                                        if (gotLitre or gotMetre): badunit = True
                                        elif (exp%3 != 0): badunit = True
                                        elif (order != (exp/3)+1) and (order != -1): badunit = True
                                        order  = (exp/3)+1
                                        if (order < 0): badunit = True
                                        gotMetre = True
                                        
                                        p_factor/=math.pow(self.__length_units, exp)
                                        p_value/=math.pow(self.__length_units, exp)
                                        
                                    elif (u.isMole()):
                                        if (gotMole): badunit = True
                                        elif (exp%1 != 0): badunit  = True
                                        elif (order != -(exp-1)) and (order != -1): badunit = True
                                        order = -(exp-1)
                                        if (order < 0): badunit = True
                                        gotMole = True
                                        
                                        p_factor/=math.pow(self.__substance_units, exp)
                                        p_value/=math.pow(self.__substance_units, exp) 
                                        
                                    elif (u.isSecond()):
                                        if (gotSecond): badunit = True
                                        elif (exp != -1):  badunit = True
                                        gotSecond = True
                                        
                                        p_factor/=math.pow(self.__time_units, exp)
                                        p_value/=math.pow(self.__time_units, exp)
                                        
                                    else: 
                                        raise NotImplementedError("Reaction '%s' parameter '%s' contains unsupported unit."%(reac.getId(), p_id))
                                        badunit = True
                                    # Breaking here means the local parameter is not a reaction parameter and will not be added to params
                                    if(badunit): 
                                        if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                        "parameter units and will cause error if in rate expression."%p_id
                                        bad_params.append([p_id, order, p_value, p_factor])    
                                        break
                                # Set first order if so
                                if (gotSecond and order == -1): order = 1    
                                # Some sanity checks: break if true and don't add param to params (print warning??)
                                if (badunit or order < 0 or gotSecond == False or p_value < 0): 
                                    if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                    "parameter units and will cause error if in rate expression."%p_id
                                    bad_params.append([p_id, order, p_value, p_factor])
                                    break
                                if order != 1: 
                                    if not (gotLitre or gotMetre): 
                                        if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                        "parameter units and will cause error if in rate expression."%p_id
                                        bad_params.append([p_id, order, p_value, p_factor])                            
                                        break
                                
                                # We may still have say (for second order reaction): m^3/s but substance units are molar. Need to convert in this case:
                                if (gotMole): params.append([p_id, order, p_value, p_factor])
                                else:
                                    # First convert to moles:
                                    p_value*=math.pow(AVOGADRO, order-1)
                                    p_factor*=math.pow(AVOGADRO, order-1)
                                    # And into SBML units, remmeber here we don't have the (usually negative) exp:
                                    p_factor*=math.pow(self.__substance_units, order-1)
                                    p_value*=math.pow(self.__substance_units, order-1)                             
                                    params.append([p_id, order, p_value, p_factor])
                                break     
                            
                            else:
                                # A 2D reaction.
                                gotTime = False
                                gotMetre = False
                                gotMole = False
                                gotSecond = False
                                order = -1
                                p_factor = 1
                                for u in units:
                                    mult = u.getMultiplier()
                                    sca = u.getScale()
                                    exp = u.getExponent()
                                    p_factor *= unitsDef_to_STEPS(mult, sca, exp)
                                    p_value *= unitsDef_to_STEPS(mult, sca, exp)
                                    if (u.isDimensionless() or u.isItem()):
                                        # OK, do nothing
                                        continue
                                    elif (u.isMetre()): 
                                        if (gotMetre): badunit = True
                                        elif (exp%2 != 0): badunit = True
                                        elif (order != (exp/2)+1) and (order != -1): badunit = True
                                        order  = (exp/2)+1
                                        if (order < 0): badunit = True
                                        gotMetre = True
                                        
                                        p_factor/=math.pow(self.__length_units, exp)
                                        p_value/=math.pow(self.__length_units, exp)
                                        
                                    elif (u.isMole()):
                                        if (gotMole): badunit = True
                                        elif (exp%1 != 0): badunit  = True
                                        elif (order != -(exp-1)) and (order != -1): badunit = True
                                        order = -(exp-1)
                                        if (order < 0): badunit = True
                                        gotMole = True
                                        
                                        p_factor/=math.pow(self.__substance_units, exp)
                                        p_value/=math.pow(self.__substance_units, exp) 
                                        
                                    elif (u.isSecond()):
                                        if (gotSecond): badunit = True
                                        elif (exp != -1):  badunit = True
                                        gotSecond = True
                                        
                                        p_factor/=math.pow(self.__time_units, exp)
                                        p_value/=math.pow(self.__time_units, exp)                                        
                                    
                                    else: 
                                        raise NotImplementedError("Reaction '%s' parameter '%s' contains unsupported unit."%(reac.getId(), p_id))
                                        badunit = True
                                    # Breaking here means the local parameter is not a reaction parameter and will not be added to params
                                    if(badunit): 
                                        if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                        "parameter units and will cause error if in rate expression."%p_id
                                        bad_params.append([p_id, order, p_value, p_factor])    
                                        break
                                # Set first order if so
                                if (gotSecond and order == -1): order = 1    
                                # Some sanity checks: break if true and don't add param to params (print warning??)
                                if (badunit or order < 0 or gotSecond == False or p_value < 0): 
                                    if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                    "parameter units and will cause error if in rate expression."%p_id
                                    bad_params.append([p_id, order, p_value, p_factor])
                                    break
                                if order != 1: 
                                    if not gotMetre: 
                                        if strict: print "WARNING: Local parameter '%s' does not have correct reaction " \
                                        "parameter units and will cause error if in rate expression."%p_id
                                        bad_params.append([p_id, order, p_value, p_factor])                            
                                        break
                                
                                # We may still have say (for second order reaction): m^3/s but substance units are molar. Need to convert in this case:
                                if (gotMole): params.append([p_id, order, p_value, p_factor])
                                else:
                                    # First convert to moles:
                                    p_value*=math.pow(AVOGADRO, order-1)
                                    p_factor*=math.pow(AVOGADRO, order-1)
                                    # And into SBML units, remmeber here we don't have the (usually negative) exp:
                                    p_factor*=math.pow(self.__substance_units, order-1)
                                    p_value*=math.pow(self.__substance_units, order-1)                             
                                    params.append([p_id, order, p_value, p_factor])
                                break
                    
                    else: 
                        # No units, assume default, which are MODEL UNITS OF VOLUME, SUBSTANCE, TIME 
                        order = 'unknown'
                        params.append([p_id, order, p_value, 1]) 
                
                kMath = kLaw.getMath()
                # Attempt to check the form of the rate in kMath to what is expected
                # Parameters are kind of a free thing - we don't know which parameter is 
                # which (if it is reversible)
                # First provide a new function to read the list and return parameters (from the 
                # global or local parameters) involved with a return value. Then construct two lists-
                # one for each combination of parameters (for a reversible reaction) of the 
                # form that is expected. If the return from one of the lists is the same as the original return 
                # we have a good reaction. 
                # One complication, as usual, is that species quantities in the expression may be 
                # amount or concentration. The concentration form is quite straight forward, but
                # if amounts the compartment volume will take a power depending on the order.
                rate_list = MLtoList(kMath, self.__function_defs)
                
                # Supply a dictionary of math parameters, but give random values to work out rates. 
                # This is because any parameter could be initialized with zero giving zero rates which would compare as equal
                # Temporary structures must be used otherwise global values will get overwritten
                
                all_variables = {}
                params_variables = {}
                
                # First the local parameters 
                paramsloc_temp = {}
                paramsloc = {}
                for p in params:
                    paramsloc[p[0]] = [p[2], p[3]]
                    paramsloc_temp[p[0]] = [random.random(), p[3]]                
                # Note: paramsloc_temp should come after global parameters because it may be a shadower
                
                params_temp = {}
                for p in self.__globalParameters:
                    params_temp[p] = [random.random(), self.__globalParameters[p][1] ]
                    
                species_temp = {}
                for s in self.__species:
                    species_temp[s] = [random.random(), self.__species[s][1] ]
                
                comps_temp = {}
                for c in self.__comps:
                    comps_temp[c] = [random.random(), self.__comps[c][1] ]
                
                patches_temp = {}
                for p in self.__patches:
                    patches_temp[p] = [random.random(), self.__patches[p][1] ]
                
                all_variables.update(params_temp)
                all_variables.update(paramsloc_temp)
                all_variables.update(species_temp)
                all_variables.update(comps_temp)
                all_variables.update(patches_temp)
                
                params_variables.update(self.__globalParameters)
                params_variables.update(paramsloc_temp)                
                
                mathreac = False
                params_id = []
                
                if strict:
                    rate = rate_func(rate_list, all_variables, params_variables, params_id, bad_params)
                    # Separate out the parameters
                    if not reversible:
                        if len(params_id) != 1:
                            if (strict): raise NotImplementedError("Reaction rate maths in reaction '%s' is not of supported form."%reac.getId())
                    else:
                        if len(params_id) != 2:
                            raise NotImplementedError("Reaction rate maths in reaction '%s' is not of supported form."%reac.getId())
                else:
                    try: rate = rate_func(rate_list, all_variables, params_variables, params_id, bad_params)
                    except:
                        mathreac = True
                    if not mathreac:
                        if not reversible:
                            if len(params_id) != 1: mathreac = True
                        else:
                            if len(params_id) != 2: mathreac = True
                
                # Update the params list
                for p_id in params_id:
                    if (p_id in bad_params):
                        if strict: raise NotImplementedError("Reaction rate maths in reaction '%s' includes" \
                            "parameter '%s' with unexpected units."%(reac.getId(), p_id))
                    if (p_id in self.__globalParameters):
                        gotP = False
                        for p in params:
                            if (p[0] == p_id): gotP = True 
                        if not gotP:
                            p_value, p_factor  = self.__globalParameters[p_id]    
                            p_order = self.__glob_params_order[p_id] 
                            if (p_order == 'not reac'): 
                                if strict: raise NotImplementedError("Reaction '%s' parameter '%s' not the right units for a reaction parameter."%(reac.getId()), p_id)
                            elif (type(p_order) in (int, float)):
                                pass
                            params.append([p_id, p_order, p_value, p_factor])
                
                # OK, now we are here we may have some parameters that have not specified units, the flag is p_order = 'unknown'
                # otherwise they have been converted already
                
                # Right, hopefully that worked (error will have been thrown if not). 
                # Now the tricky part of constructing the expected list
                
                # Lets treat forward and backward separately (for reversible reactions), but create a holder
                # for a reverse reaction regardless (initialised to 0, ignored for non-reversible reactions)
                
                # We may already know that the raection has bad parameters etc and is a 'math reac'
                if not mathreac:
                    correct_rate_list_1 = ['-', [None, 0.0]]
                    
                    if (reversible): correct_rate_list_2 = ['-', [None, None]]
                    
                    # just a reminder reacts_f is a list of the reactant species by string id for the forward reaction
                    # and reacts_b is a list of reactant species for the backward reaction (if there is one)
                    # reac_comp is the compartment id
                    
                    for_vars_1 = [reactant_space_id, params_id[0]]+reacts_f
                    correct_rate_list_1[1][0] = make_rate_list(for_vars_1, self.__species_subst_units, reactant_space_id) 
                    if (reversible):
                        back_vars_1 = [product_space_id, params_id[1]]+reacts_b
                        correct_rate_list_1[1][1] = make_rate_list(back_vars_1, self.__species_subst_units, product_space_id)                    
                    if (reversible):
                        for_vars_2 = [reactant_space_id, params_id[1]]+reacts_f
                        correct_rate_list_2[1][0] = make_rate_list(for_vars_2, self.__species_subst_units, reactant_space_id)
                        back_vars_2 = [product_space_id, params_id[0]]+reacts_b
                        correct_rate_list_2[1][1] = make_rate_list(back_vars_2, self.__species_subst_units, product_space_id)                    
                    
                    correct_rate_1 = MLfunc(correct_rate_list_1, all_variables)

                    if (reversible): 
                        correct_rate_2 = MLfunc(correct_rate_list_2, all_variables)
                    
                    setparam_f = False
                    if (reversible): setparam_b = False
                    
                    if _float_approx_equal(rate, correct_rate_1):
                        # We have found the right parameters. Set them to the correct reaction
                        # We've got the ids, but need to find the right parameter in the params list
                        for p in params:
                            if not setparam_f and p[0] == params_id[0]:
                                if (p[1] == 'not reac'):
                                    if strict: raise NotImplementedError("Reaction rate parameter for reaction '%s' has unexpected units."%reac.getId())
                                    else: 
                                        mathreac = True
                                        break
                                # 'unknown' basically means that the units were not specified, so we have to use model units
                                else:
                                    order = r_f.getOrder()
                                    # I am not clear about what to do here - what are the default units for say a second order reaction?
                                    # Guess I am going to assume volume units are the same as the model volume units, though 
                                    # documentation is unclear on this
                                    kvalue = p[2]
                                    factor = 1.0 /(self.__time_units*math.pow(self.__substance_units, order-1))
                                    if not reac_type == '2D surface reaction':
                                        factor *= (math.pow(self.__volume_units*1000.0, order-1))
                                    else: 
                                        factor *= (math.pow(self.__area_units, order-1))
                                    # Not sure what to do about the mole - there seems to be no global substance units
                                    
                                    p[2] = kvalue*factor
                                    p[3] = factor
                                    
                                r_f.setKName(p[0])
                                r_f.setKValue(p[2])
                                #r_f.setUnitsConv(p[3]) 
                                if reac_type == '2D surface reaction': r_f.setType('2D')
                                else: r_f.setType('3D') 
                                surface_reactions.append(r_f)
                                setparam_f  = True
                            elif (reversible):
                                if (not setparam_b and p[0] == params_id[1]):
                                    if (p[1] == 'not reac'): # or p[1] != len(r_b.getLhs())):
                                        if strict: raise NotImplementedError("Reaction rate parameter for reversible reaction '%s' has unexpected units."%reac.getId())
                                        else:
                                            mathreac = True
                                            break
                                    # 'unknown' basically means that the units were not specified, so we have to use model units
                                    else:
                                        order = r_b.getOrder()
                                        kvalue = p[2]
                                        
                                        factor = 1.0 / (self.__time_units*math.pow(self.__substance_units, order-1))
                                        if not reac_type == '2D surface reaction':
                                            factor *= (math.pow(self.__volume_units*1000.0, order-1))
                                        else: 
                                            factor *= (math.pow(self.__area_units, order-1))
                                        
                                        p[2] = kvalue*factor    
                                        p[3] = factor
                                    
                                    r_b.setKName(p[0])
                                    r_b.setKValue(p[2])
                                    #r_b.setUnitsConv(p[3])
                                    if reac_type == '2D surface reaction': r_b.setType('2D')
                                    else: r_b.setType('3D')
                                    surface_reactions.append(r_b)
                                    setparam_b = True
                    elif (reversible):
                        if _float_approx_equal(rate, correct_rate_2):
                            # We have found the right params
                            for p in params:
                                if not setparam_f and p[0] == params_id[1]:
                                    if (p[1] == 'not reac'): # or p[1] != len(r_f.getLhs())):
                                        if strict: raise NotImplementedError("Reaction rate parameter for reaction '%s' has unexpected units."%reac.getId())
                                        else:
                                            mathreac = True
                                            break
                                    # 'unknown' basically means that the units were not specified, so we have to use model units
                                    else:
                                        order = r_f.getOrder()
                                        kvalue = p[2]
                                        
                                        factor = 1.0 / (self.__time_units*math.pow(self.__substance_units, order-1))
                                        if not reac_type == '2D surface reaction':
                                            factor *= (math.pow(self.__volume_units*1000.0, order-1))
                                        else: 
                                            factor *= (math.pow(self.__area_units, order-1))
                                        
                                        p[2] = kvalue*factor
                                        p[3] = factor
                                        
                                    r_f.setKName(p[0])
                                    r_f.setKValue(p[2])
                                    #r_f.setUnitsConv(p[3])
                                    if reac_type == '2D surface reaction': r_f.setType('2D')
                                    else: r_f.setType('3D')
                                    surface_reactions.append(r_f)
                                    setparam_f  = True
                                elif (not setparam_b and p[0] == params_id[0]):
                                    if (p[1] == 'not reac'):# or p[1] != len(r_b.getLhs())):
                                        if strict: raise NotImplementedError("Reaction rate parameter order for reversible reaction '%s' has unexpected units."%reac.getId())
                                        else:
                                            mathreac = True
                                            break
                                    # 'unknown' basically means that the units were not specified, so we have to use model units
                                    else:
                                        order = r_b.getOrder()
                                        kvalue = p[2]
                                        
                                        factor = 1.0 /(self.__time_units*math.pow(self.__substance_units, order-1))
                                        if not reac_type == '2D surface reaction':
                                            factor *= (math.pow(self.__volume_units*1000.0, order-1))
                                        else:
                                            factor *= (math.pow(self.__area_units, order-1))
                                        
                                        p[2] = kvalue*factor
                                        p[3] = factor
                                    
                                    r_b.setKName(p[0])
                                    r_b.setKValue(p[2])
                                    #r_b.setUnitsConv(p[3])
                                    if reac_type == '2D surface reaction': r_b.setType('2D')
                                    else: r_b.setType('3D')
                                    surface_reactions.append(r_b)
                                    setparam_b = True    
                        else:
                            if strict: raise NotImplementedError("Reaction rate maths in reaction %s is not of expected rate."%reac.getId())
                            else: mathreac = True
                    else:
                        if strict: raise NotImplementedError("Reaction rate maths in reaction %s is not of expected rate."%reac.getId())
                        else: mathreac = True
                
                # Separate treatment for "math reacs"
                if mathreac:
                    # Assuming that it makes no difference whether the reaction is reversilbe or not if the maths is strange
                    r_f.setKName('unknown')
                    r_f.setKValue(0.0)
                    if reac_type == '2D surface reaction': r_f.setType('2D')
                    else: r_f.setType('3D')
                    
                    # This "correct rate" will be just the spcies and the comp- no parameter.
                    # The changing parameter will then be the actual maths divided by this rate list
                    if not reac_type == '2D surface reaction':
                        for_vars_1 = [reac_comp1_id]+reacts_f
                    else: 
                        for_vars_1 = [reac_surf_id]+reacts_f

                    # Add the "bad parameters"                
                    for p in bad_params:
                        paramsloc[p[0]] = [p[2], p[3]]
                        
                    if not reac_type == '2D surface reaction': 
                        ord_tmp = r_f.getOrder()
                        # Factor is going to convert the returned rate in base units to SBML units
                        factor = 1.0*math.pow(conv_factor_space, ord_tmp-1)
                        factor /= math.pow(conv_factor_subs, ord_tmp-1)
                        surface_math_reactions.append([r_f, rate_list, make_rate_list(for_vars_1, self.__species_subst_units, reac_comp1_id), paramsloc, factor])
                    else:
                        ord_tmp = r_f.getOrder()
                        factor = 1.0*math.pow(conv_factor_space, ord_tmp-1)
                        factor /= math.pow(conv_factor_subs, ord_tmp-1)
                        surface_math_reactions.append([r_f, rate_list, make_rate_list(for_vars_1, self.__species_subst_units, reac_surf_id), paramsloc, factor])
                    
                    surface_reactions.append(r_f)
                            
        return reactions, math_reactions, surface_reactions, surface_math_reactions
    
    ################################################################################################
    
    def _parseGlobalParameters(self):
        """
        Import the global parameters.
        
        Note::  Will read all global parameters which may or may not be reaction parameters
                Frustratingly, there is no specification as to whether a parameter is a reaction
                parameter in SBML so we have to find out. The 'order' attribute for a 
                parameter will be set to 'not reac' if it is thought not to be a reaction parameter
                and 'unknown' if units are unspecified basically. This info is then avaialble for reactions.
        """
        
        listOfParameters = self.__model.getListOfParameters()
        globPar = {}
        globPar_order = {}
        for par in listOfParameters:
            p_id = str(par.getId())
            p_value = par.getValue()
            p_units = par.getUnits() 
            # A factor to convert the value to STEPS units
            p_factor = 1 
            if (p_units): 
                # Units for a reaction should be (conc)^-(n-1).(time)^-1
                # Where n gives the order of the reaction
                gotTime = False
                gotLitre = False
                gotMetre = False
                gotMole = False
                gotSecond = False
                order = -1
                notreacunit = False
                
                # The problem here is that parameters don't have to have units specified and can 
                # be anything! Ok, we can assume and store the value in SBML units, but cannot
                # convert to anything useful for us until we know the dimensions of the parameter, 
                # e.g. if it is a second-order reaction parameter we need to convert by substance, volume 
                # and time model units. If the conversion cannot be done yet we leave it blank and check later. 
                pfactor = 1
                unitdef = self.__model.getUnitDefinition(p_units) 
                if (unitdef):
                    units = unitdef.getListOfUnits()  
                    # The following effort is all because of a simple sbml simplicity for modellers:
                    # the reversible reactions, without having to specify which parameter is which. Try to find 
                    # the order if it is a reaction constant so as to make it possible to assign.
                    for u in units:
                        mult = u.getMultiplier()
                        sca = u.getScale()
                        exp = u.getExponent()
                        p_value *= unitsDef_to_STEPS(mult, sca, exp)
                        p_factor*= unitsDef_to_STEPS(mult, sca, exp)
                        if (u.isDimensionless() or u.isItem()):
                            # OK, do nothing
                            continue
                        elif (u.isLitre()):
                            # Convert to m
                            p_factor*=math.pow(0.001, exp)
                            p_value*=math.pow(0.001, exp)
                            
                            # now convert to SBML units
                            p_factor/=math.pow(self.__volume_units, exp)
                            p_value/=math.pow(self.__volume_units, exp)    
                                                    
                            if (gotLitre or gotMetre): notreacunit = True
                            if (exp%1 != 0): notreacunit = True
                            if (order != exp+1) and (order != -1): notreactunit = True
                            order = exp+1
                            if (order < 0): notreacunit = True
                            gotLitre = True
                        elif (u.isMetre()):
                            if (gotLitre or gotMetre): notreacunit = True
                            if (exp%3 != 0): notreacunit = True
                            if (order != (exp/3)+1) and (order != -1): notreacunit = True
                            order  = (exp/3)+1
                            if (order < 0): notreacunit = True
                            gotMetre = True
                            
                            p_factor/=math.pow(self.__length_units, exp)
                            p_value/=math.pow(self.__length_units, exp)
                            
                        elif (u.isMole()):
                            if (gotMole): notreacunit = True
                            if (exp%1 != 0): notreacunit  = True
                            if (order != -(exp-1)) and (order != -1): notreacunit = True
                            order = -(exp-1)
                            if (order < 0): notreacunit = True
                            gotMole = True
                            
                            p_factor/=math.pow(self.__substance_units, exp)
                            p_value/=math.pow(self.__substance_units, exp) 
                            
                        elif (u.isSecond()):
                            if (gotSecond): notreacunit = True
                            if (exp != -1): notreacunit = True
                            gotSecond = True
                            
                            p_factor/=math.pow(self.__time_units, exp)
                            p_value/=math.pow(self.__time_units, exp) 
                            
                        else: 
                            raise NotImplementedError("Parameter '%s' contains unsupported unit"%p_id)
                            notreacunit = True
                # Set first order if so
                if (gotSecond and order == -1): order =1
                if (notreacunit): 
                    globPar[p_id] = [p_value, p_factor]
                    globPar_order[p_id] = 'not reac'                   
                elif (order < 0 or gotSecond == False or p_value < 0): 
                    # Also not a reac unit
                    globPar[p_id] = [p_value, p_factor]
                    globPar_order[p_id] = 'not reac' 

                elif (order != 1 and not (gotLitre or gotMetre)):
                        # Also not a reac unit
                        globPar[p_id] = [p_value, p_factor]
                        globPar_order[p_id] = 'not reac'                  
                else:
                    # We may still have say (for second order reaction): m^3/s but substance units are molar. Need to convert in this case:
                    if (gotMole): 
                        globPar[p_id] = [p_value, p_factor]
                        globPar_order[p_id] = order
                    else:
                        # First convert to moles if this is specifically a reaction parameter:
                        p_value*=math.pow(AVOGADRO, order-1)
                        p_factor*=math.pow(AVOGADRO, order-1)
                        # And into SBML units, remmeber here we don't have the (usually negative) exp:
                        p_factor*=math.pow(self.__substance_units, order-1)
                        p_value*=math.pow(self.__substance_units, order-1)   
                        globPar[p_id] = [p_value, p_factor]  
                        globPar_order[p_id] = order             

            else: 
                # No units, assume default, whatever they are 
                globPar[p_id] = [p_value, 1.0]
                globPar_order[p_id] = 'unknown'
                print "WARNING: Parameter '%s' units are undefined and will be assumed to be model units."%p_id
        
        return globPar, globPar_order
    
    ################################################################################################
    
    def _setupReactions(self):
        """
        Create the reactions in STEPS
        """
        for r in self.__reactions:
            kreac = smodel.Reac(r.getName(), self.__mdl.getVolsys(r.getVsys()), lhs = r.getLhs(), rhs = r.getRhs())
            # This is a special case where the reaction parameter is already in STEPS units
            kreac.kcst = r.getKValue()
        
        for sr in self.__surface_reactions:
            kreac = smodel.SReac(sr.getName(), self.__mdl.getSurfsys(sr.getSsys()), \
                ilhs = sr.getILhs(), slhs = sr.getSLhs(), olhs = sr.getOLhs(), \
                irhs = sr.getIRhs(), srhs = sr.getSRhs(), orhs = sr.getORhs())
            kreac.kcst = sr.getKValue()
        
    ################################################################################################
    
    def _parseevnts(self):
        """
        Import the sbml Events.
        """
        
        listOfevnts = self.__model.getListOfEvents()
        
        #  Events: evnt_id:[var_id, math_list]
        evnts_trig = {}
        # Assignment rules for evnts: evnt_id:[[var_id, math_list], ... ]
        evnts_ass = {}
        # Delays for event - [the unchanging maths, delay time calculated when event fires]
        evnts_dl = {}
        # To store whether a trigger has gone from "False" to "True"
        evnts_flip = {}   
        # Stores whether to evaluate assignments as trigger time, or execution time
        evnts_trigvals = {}     

        for evnt in listOfevnts:
            # First a table storing whether the event has 'flipped' from false to true - 
            # required to fire an event
            evnts_flip[evnt.getId()] = False
            trigger = evnt.getTrigger()
            trig_math = trigger.getMath()
            # Convert math into list, don't forget to give the function defs
            trig_list = MLtoList(trig_math, self.__function_defs)
            evnts_trig[evnt.getId()] = trig_list
            evnts_trigvals[evnt.getId()] = evnt.getUseValuesFromTriggerTime()
            evnts_ass[evnt.getId()] = []
            for eass in evnt.getListOfEventAssignments():
                var = eass.getVariable()
                eass_math = eass.getMath()
                eass_list = MLtoList(eass_math, self.__function_defs)
                evnts_ass[evnt.getId()].append([var, eass_list])
            
            delay = evnt.getDelay()
            if (delay):
                dl_math = delay.getMath()
                dl_list = MLtoList(dl_math, self.__function_defs)
                evnts_dl[evnt.getId()] = [dl_list, 0.0]
            else: evnts_dl[evnt.getId()] = [0.0, 0.0]
            
        return evnts_trig, evnts_ass, evnts_dl, evnts_flip, evnts_trigvals
            
    ################################################################################################

    def _parseRules(self):
        """
        Import the rules.
        """
        
        listOfRules = self.__model.getListOfRules()
        
        # Assignment rules, id:[var_id, value]
        ass_rules = {}
        # Rate rules: id:[var_id, rate (in /time)]
        rate_rules = {}
        
        for rule in listOfRules:
            if (rule.isAlgebraic()):
                raise NotImplementedError("Rule '%s' is Algebraic (unsupported)."%rule.getId())
            var = rule.getVariable()   
            if (var in self.__comps or var in self.__globalParameters or var in self.__species):
                rmath = rule.getMath()
                rmath_list = MLtoList(rmath, self.__function_defs)
            else:
                # Special case for species reference (really what is the point in this silly thing?!)
                # I am not going to allow changing stoichiometry, but allow assignment of the stoichiometry
                # at the beginning
                variables = {}
                variables.update(self.__time)
                variables.update(self.__species)
                variables.update(self.__globalParameters)
                variables.update(self.__comps)    
                variables.update(self.__patches)
                rmath = rule.getMath()
                rmath_list = MLtoList(rmath, self.__function_defs)
                value = MLfunc(rmath_list, variables)
                self.__spec_refs[var] = value
                # Don't want to add to the ass_rules etc
                continue
            if (rule.isRate() == True): rate_rules[rule.getId()] = [var, rmath_list]
            elif (rule.isAssignment() == True): ass_rules[rule.getId()] = [var, rmath_list] 
            else:
                assert(False)
        
        return ass_rules, rate_rules
                
    ################################################################################################

    def _parseInitialAssignments(self):
        """
        Import initial assignments.
        """
        
        listOfInitialAssignments = self.__model.getListOfInitialAssignments()
        
        # Keeping any changed varaibles as they are during loop, then update after
        update_specs={}
        update_comps={}
        update_patches={}
        update_params={}
        update_spec_ref={}

        variables = {}
        variables.update(self.__time)
        variables.update(self.__species)
        variables.update(self.__globalParameters)
        variables.update(self.__comps)
        variables.update(self.__patches)

        for initass in listOfInitialAssignments:
            var = initass.getSymbol()
            math = initass.getMath() 
            # Convert the math to a list
            ia_math_frm = MLtoList(math, self.__function_defs)
            
            # Remember the return will be in SBML units
            value = MLfunc(ia_math_frm, variables)
            if (value!=value): raise NotImplementedError("Initial assignment '%s' value is not-a-number."%initass.getId())
            if (var in self.__comps):
                # Convert to SBML units
                value*=self.__comps[var][1]                
                update_comps[var] = value
                
                # Now update any species that, complicatedly, may be intended to be amounts injected
                # but are hasonlysubstanceunits flaf is false, meaning it is a conc in expressions so stored as concs
                # The problem that a species may also be changed in these initial assignments and the ordering
                # is just too complicated a problem to even attempt to solve 
                
                for specie in self.__species_comps:
                    if self.__species_comps[specie] == var:
                        if (self.__species_amount_flag[specie] == True and self.__species_subst_units[specie] == False):
                            # We have an intended injected amount, but stored as conc, and volume has changed
                            # Don't do anything if already changed by initial assignment
                            if (update_specs.has_key(specie)): continue
                            update_specs[specie]=self.__species[specie][0]*(self.__comps[var][0]/value)
                        if (self.__species_amount_flag[specie] == False and self.__species_subst_units[specie] == True):
                            # We have an intended injected conc, but stored as amount, and volume has changed
                            # Don't do anything if already changed by initial assignment
                            if (update_specs.has_key(specie)): continue
                            update_specs[specie]=self.__species[specie][0]*(value/self.__comps[var][0])
            elif var in self.__patches:
                value*=self.__patches[var][1]
                update_patches[var] = value                
                for specie in self.__species_comps:
                    if self.__species_comps[specie] == var:
                        if (self.__species_amount_flag[specie] == True and self.__species_subst_units[specie] == False):
                            # We have an intended injected amount, but stored as conc, and area has changed
                            # Don't do anything if already changed by initial assignment
                            if (update_specs.has_key(specie)): continue
                            update_specs[specie]=self.__species[specie][0]*(self.__patches[var][0]/value)
                        if (self.__species_amount_flag[specie] == False and self.__species_subst_units[specie] == True):
                            # We have an intended injected conc, but stored as amount, and volume has changed
                            # Don't do anything if already changed by initial assignment
                            if (update_specs.has_key(specie)): continue
                            update_specs[specie]=self.__species[specie][0]*(value/self.__patches[var][0]) 
            elif (var in self.__species):
                # Convert to SBML units
                value*=self.__species[var][1]
                update_specs[var] = value

            elif (var in self.__globalParameters):
                if self.__globalParameters[var][1]: value*=self.__globalParameters[var][1]
                update_params[var] = value
            else:
                # Something else, must be species reference. 
                update_spec_ref[var] = value
        
        for us in update_specs:
            self.__species[us][0] = update_specs[us]
        for uc in update_comps:
            self.__comps[uc][0] = update_comps[uc]
        for up in update_patches:
            self.__patches[up][0] = update_patches[up]
        for up in update_params:
            self.__globalParameters[up][0] = update_params[up]
        for sr in update_spec_ref:
            self.__spec_refs[sr] = update_spec_ref[sr]
            
    ################################################################################################
            
    def getModel(self):
        """
        Returns a reference to the steps.model.Model biocehmical model container object 
        created during SBML import.
        
        Syntax::
        
            getModel()
        
        Arguments:
            None
        
        Return:
            steps.model.Model
        
        """
        return self.__mdl

    ################################################################################################
    
    def getGeom(self):
        """
        Returns a reference to the steps.geom.Geom geoemtry container object 
        created during SBML import.
        
        Syntax::
        
            getGeom()
        
        Arguments:
            None
        
        Return:
            steps.geom.Geom
        
        """
        return self.__geom

    ################################################################################################
        
    def _getEvents(self):
        return self.__evnts_lt, self.__evnts_gt, self.__evnts_ass
    
    ################################################################################################

    def setupSim(self, sim):
        """
        Setup the simulation, initialising molecule counts, compartment volumes and 
        SBML 'boundary conditions" to those specified in the SBML file, 
        which may differ from default values through mechanisms such as Initial Assignments.
        This function is NOT called internally and must be called by the user.
        
        Syntax::
        
            setupSim(sim)
        
        Arguments:
            steps.solver.API sim (steps.solver.Wmdirect or steps.solver.Wmrk4 solver object)
        
        Return 
            None
        
        """
        
        # REMEMBER ALL VALUES ARE STORED IN SBML MODEL UNITS AND MUST BE CONVERTED TO STEPS UNITS
        
        # Need to check the species_amount flag because even though the substanceunits might be false (meaning
        # the substance should be treated as a conc in expressions) it still might be an amount injected initially, 
        # and vice versa
        
        for comp in self.__comps:
            try: 
                sim.setCompVol(comp, self.__comps[comp][0]*self.__volume_units)
            except: 
                raise NotImplementedError("Failed to set compartment volume during"\
                " Initial Assignment (negative volume?).")
        for patch in self.__patches:
            try: 
                sim.setPatchArea(patch, self.__patches[patch][0]*self.__area_units)
            except: 
                raise NotImplementedError("Failed to set patch area during"\
                " Initial Assignment (negative area?).")      
        for specie in self.__species:
            if self.__species_comps[specie] in self.__comps:
                if (self.__species_subst_units[specie] == True):
                    try: sim.setCompAmount(self.__species_comps[specie], specie, self.__species[specie][0]*self.__substance_units)
                    except: 
                        raise NotImplementedError("Cannot set amount of species '%s" \
                        "' in compartment '%s'. Species may not be reactant or "\
                        "product in any comparment reactions, or molecules count may"\
                        " be larger than maximum unsigned integer "\
                        "in stochastic simulation."%(specie, self.__species_comps[specie]) ) 
                else:
                    # Remember the conversion is to meters, steps needs molar
                    try :
                        sim.setCompConc(self.__species_comps[specie], specie, self.__species[specie][0]*(self.__substance_units/(self.__volume_units*1.0e3)))
                    except: 
                        raise NotImplementedError("Cannot set concentration of species '%s" \
                        "' in compartment '%s'. Species may not be reactant or "\
                        "product in any comparment reactions, or molecules count may"\
                        " be larger than maximum unsigned integer "\
                        "in stochastic simulation."%(specie, self.__species_comps[specie]) )
                const_bcs = self.__species_const_bc[specie]
                if (const_bcs[0] == True):
                    try : sim.setCompClamped(self.__species_comps[specie], specie, True)
                    except: pass
            elif self.__species_comps[specie] in self.__patches:
                if (self.__species_subst_units[specie] == True):
                    try: sim.setPatchAmount(self.__species_comps[specie], specie, self.__species[specie][0]*self.__substance_units)
                    except: 
                        raise NotImplementedError("Cannot set amount of species '%s" \
                        "' in compartment '%s'. Species may not be reactant or "\
                        "product in any comparment reactions, or molecules count may"\
                        " be larger than maximum unsigned integer "\
                        "in stochastic simulation."%(specie, self.__species_comps[specie]) ) 
                else:
                    # We want to inject conc, but can't for patch
                    try :sim.setPatchAmount(self.__species_comps[specie], specie, self.__species[specie][0]*self.__substance_units*(sim.getPatchArea(self.__species_comps[specie])/self.__area_units))
                    except: 
                        raise NotImplementedError("Cannot set concentration of species '%s" \
                        "' in patch '%s'. Species may not be reactant or "\
                        "product in any patch reactions, or molecules count may"\
                        " be larger than maximum unsigned integer "\
                        "in stochastic simulation."%(specie, self.__species_comps[specie]) )
                const_bcs = self.__species_const_bc[specie]
                if (const_bcs[0] == True):
                    try : sim.setPatchClamped(self.__species_comps[specie], specie, True)
                    except: pass
            else: assert (False)
                
        
        # NOTE: No need to do reactions because the parameters have been changed to any initial 
        # Assignments by setupReactions()
                
    ################################################################################################

    def updateSim(self, sim, simdt):
        """
        Update the simulation solver state, which may impact any variables in the simulation that 
        can be altered within Rules and Events. Time since last update given as argument dt in seconds. 
        
        Syntax:: 
        
            updateSim(sim, dt)
        
        Arguments:
            * steps.solver.API  sim (steps.solver.Wmdirect or steps.solver.Wmrk4 solver object)
            * float             dt
            
        Return 
            None
        """
        # Update the spcies dictionary with the species concentrations. No need to update comp volumes
        # because these are not changed during simulation
        for s in self.__species:
            if self.__species_comps[s] in self.__comps:
                if (self.__species_subst_units[s]): 
                    # Now directly manipulate the self.__species object
                    # We need to convert to SBML units
                    self.__species[s][0] = sim.getCompAmount(self.__species_comps[s], s)/self.__substance_units
                else: 
                    self.__species[s][0] = (sim.getCompConc(self.__species_comps[s], s)*1.0e3*self.__volume_units)/self.__substance_units
            elif self.__species_comps[s] in self.__patches:
                if (self.__species_subst_units[s]):
                    self.__species[s][0] = sim.getPatchAmount(self.__species_comps[s], s)/self.__substance_units
                else:
                    # We want to get conc, but can't for patch
                    self.__species[s][0] = sim.getPatchAmount(self.__species_comps[s], s)/(self.__substance_units*(sim.getPatchArea(self.__species_comps[s])/self.__area_units))
        # Adding time. 
        newtime = sim.getTime()/self.__time_units
        self.__time['time'][0] = self.__time['t'][0] = self.__time['s'][0] = self.__time['Time'][0] = newtime
        
        # Collect all the variables in a dictionary on the fly-easier to group them for one dict to MLfunc
        variables = {}
        variables.update(self.__time)
        variables.update(self.__species)
        variables.update(self.__globalParameters)
        variables.update(self.__comps)
        variables.update(self.__patches)
        # I guess updating the rates should come first        
        poplist = []
        
        # Not easy to get the parameters from STEPS, so store a dictionary of params and new value
        params_update = {}
        # Lets do it for comps and patches too so as to save some time
        comps_update={}
        patches_update={}
        specs_update={}
        
        # NOTE: Shouldn't update the global value sduring the loop, because of the sequencing thing they should be updated afterwards
        for r_rate in self.__rules_rate:
            var = self.__rules_rate[r_rate][0]
            # MLfunc returned value should be / (timeunits) so convert to /s
            # This is a rate so we multiply by sim dt
            value = MLfunc(self.__rules_rate[r_rate][1], variables)*(simdt/self.__time_units)
            if (var in self.__species):
                # Convert to SBML units
                value*=self.__species[var][1]
                
                # First update the specs:
                specs_update[var] = self.__species[var][0]+value
                # We have to wait see if we have an amount or conc
                if self.__species_comps[var] in self.__comps:
                    if (self.__species_subst_units[var] == True):
                        try: 
                            before_amount = sim.getCompAmount(self.__species_comps[var], var)
                            after_amount = value*self.__substance_units
                            sim.setCompAmount(self.__species_comps[var], var, before_amount+after_amount)
                        except: 
                            raise NotImplementedError("Failed to set amount during rate-rule '%s' (negative amount?)."%r_rate)
                    else:
                        try: 
                            before_conc = sim.getCompConc(self.__species_comps[var], var)
                            after_conc = value*(self.__substance_units/self.__volume_units)*1.0e-3
                            sim.setCompConc(self.__species_comps[var], var, before_conc+after_conc)
                        except: 
                            raise NotImplementedError("Failed to set concentration during rate-rule '%s' (negative conc?)."%r_rate)
                elif self.__species_comps[var] in self.__patches:
                    if (self.__species_subst_units[var] == True):
                        try:
                            before_amount = sim.getPatchAmount(self.__species_comps[var], var)
                            after_amount = value*self.__substance_units
                            sim.setPatchAmount(self.__species_comps[var], var, before_amount+after_amount)
                        except:
                            raise NotImplementedError("Failed to set amount during rate-rule '%s' (negative amount?)."%r_rate)
                    else:
                        try:
                            before_conc = sim.getPatchAmount(self.__species_comps[var], var)/sim.getPatchArea(self.__species_comps[var])
                            after_conc = value*(self.__substance_units/self.__area_units)
                            sim.setPatchAmount(self.__species_comps[var], var, (before_conc+after_conc)*sim.getPatchArea(self.__species_comps[var]))
                        except:
                            raise NotImplementedError("Failed to set concentration during rate-rule '%s' (negative conc?)."%r_rate)                            
            elif (var in self.__comps):
                # Convert to SBML units
                value*=self.__comps[var][1] 
                # Convert returned value in SBML units to STEPS units
                new_val = self.__comps[var][0]+value
                comps_update[var] = new_val
                new_val*=self.__volume_units
                try: sim.setCompVol(var, new_val)
                except:
                    raise NotImplementedError("Failed to set compartment volume during rate-rule '%s' (negative volume?)."%r_rate)
            elif (var in self.__patches):
                # we want to keep in SBML units now:
                value*=self.__patches[var][1]
                # Convert returned value in SBML units to STEPS units
                new_val = self.__patches[var][0]+value
                patches_update[var] = new_val
                new_val*=self.__area_units
                try: sim.setPatchArea(var, new_val)
                except:
                    raise NotImplementedError("Failed to set patch area during rate-rule '%s' (negative area?)."%r_rate)                
            elif (var in self.__globalParameters):                
                if self.__globalParameters[var][1]: value*=self.__globalParameters[var][1]
                
                # Now need to convert to molar, STEPS units from metre units
                order = self.__glob_params_order[var]

                setorder = False
                if (type(order) in (int, float)): setorder = True
                # Need to try and get all the reactions with this parameter and change- store in a list first
                # Store the order too for scaling from metres to litres
                k_reacs_comp = {}
                for r in self.__reactions:
                    if (r.getKName() == var): 
                        k_reacs_comp[r.getName()] = r.getCompID() 
                        ord_tmp = r.getOrder()
                        if setorder: assert ord_tmp == order 
                        else: 
                            order = ord_tmp
                            setorder = True

                srtype=''
                k_sreacs_patch = {}
                for sr in self.__surface_reactions:
                    if (sr.getKName() == var): 
                        k_sreacs_patch[sr.getName()] = sr.getPatchID()
                        ord_tmp = sr.getOrder()
                        if setorder: assert ord_tmp == order 
                        else: 
                            order = ord_tmp
                            setorder = True
                        type_temp = sr.getType()
                        if  srtype: assert type_temp == srtype
                        else:
                            srtype = type_temp
                
                # Value should already be in SBML units
                params_update[var] = [value+self.__globalParameters[var][0]]
                
                if (k_reacs_comp): 
                    assert setorder

                    value*= self.__param_converter_vol[order]
                    for k_r in k_reacs_comp:
                        compid = k_reacs_comp[k_r]
                        before_val = sim.getCompReacK(compid, k_r)
                        try: sim.setCompReacK(compid, k_r, before_val+value)
                        except: raise NotImplementedError("Negative constant of Reaction '%s'"%k_r)
                
                if (k_sreacs_patch): 
                    assert setorder
                    assert srtype
                    # Finally can convert to STEPS units

                    if srtype == '3D': value*= self.__param_converter_vol[order]
                    else: value*= self.__param_converter_area[order]
                    for k_sr in k_sreacs_patch:
                        patchid = k_sreacs_patch[k_sr]
                        before_val = sim.getPatchSReacK(patchid, k_sr)
                        try: sim.setPatchSReacK(patchid, k_sr, before_val+value)
                        except: raise NotImplementedError("Negative constant of Surface Reaction '%s'"%k_sr)

        
        for p in poplist:
            self.__rules_rate.pop(p)
        # Now assignment rules: just like rate rules, but no need to multiply by a dt
        poplist = []
        
        for r_ass in self.__rules_ass:
            var = self.__rules_ass[r_ass][0]
            # Returned value is in given units
            value = MLfunc(self.__rules_ass[r_ass][1], variables)
            if (var in self.__species):
                # First convert to SBML units
                value*=self.__species[var][1]
                specs_update[var] = value
                if self.__species_comps[var] in self.__comps:
                    if (self.__species_subst_units[var] == True):
                        try: 
                            sim.setCompAmount(self.__species_comps[var], var, value*self.__substance_units)
                        except: 
                            raise NotImplementedError("Failed to set amount during assignment-rule '%s' (negative amount?)."%r_ass)
                    else:
                        try: 
                            sim.setCompConc(self.__species_comps[var], var, value*(self.__substance_units/self.__volume_units)*1.0e-3)
                        except: 
                            raise NotImplementedError("Failed to set concentration during assignment-rule '%s' (negative conc?)."%r_ass)
                elif self.__species_comps[var] in self.__patches:
                    if (self.__species_subst_units[var] == True):
                        try:
                            sim.setPatchAmount(self.__species_comps[var], var, value*self.__substance_units)
                        except:
                            raise NotImplementedError("Failed to set amount during assignment-rule '%s' (negative amount?)."%r_ass)
                    else:
                        try:
                            sim.setPatchAmount(self.__species_comps[var], var, value*(self.__substance_units/self.__area_units)*sim.getPatchArea(self.__species_comps[var]))
                        except:
                            raise NotImplementedError("Failed to set concentration during assignment-rule '%s' (negative conc?)."%r_ass)         
            elif (var in self.__comps):
                # Convert returned value to SBML units
                value*=self.__comps[var][1]
                comps_update[var] = value
                # Convert  value in SBML units to STEPS units
                value*=self.__volume_units
                try: sim.setCompVol(var, value)
                except:
                    raise NotImplementedError("Failed to set compartment volume during assignment-rule '%s' (negative volume?)."%r_ass)
            elif (var in self.__patches):
                # Convert returned value to SBML units
                value*=self.__patches[var][1]
                patches_update[var] = value
                # Convert returned value in SBML units to STEPS units
                value*=self.__area_units
                try: sim.setPatchArea(var, value)
                except:
                    raise NotImplementedError("Failed to set patch area during assignment-rule '%s' (negative area?)."%r_ass)
            elif (var in self.__globalParameters):
                # Convert value from given units to SBML units
                if self.__globalParameters[var][1]: value*=self.__globalParameters[var][1]
                        
                # Note: Not updating the global parameters here because they may appear in other rate rules
                # params[var] = value
                # Now need to convert to molar, STEPS units from metre units
                order = self.__glob_params_order[var]

                setorder = False
                if (type(order) in (int, float)): setorder = True
                # Need to try and get all the reactions with this parameter and change- store in a list first
                # Store the order too for scaling from metres to litres
                k_reacs_comp = {}
                for r in self.__reactions:
                    if (r.getKName() == var): 
                        k_reacs_comp[r.getName()] = r.getCompID() 
                        ord_tmp = r.getOrder()
                        if setorder: assert ord_tmp == order 
                        else: 
                            order = ord_tmp
                            setorder = True
                
                srtype=''
                k_sreacs_patch = {}
                for sr in self.__surface_reactions:
                    if (sr.getKName() == var): 
                        k_sreacs_patch[sr.getName()] = sr.getPatchID()
                        ord_tmp = sr.getOrder()
                        if setorder: assert ord_tmp == order 
                        else: 
                            order = ord_tmp
                            setorder = True
                        type_temp = sr.getType()
                        if  srtype: assert type_temp == srtype
                        else:
                            srtype = type_temp                
                
                params_update[var] = [value]

                if (k_reacs_comp): 
                    assert setorder
                    value*= self.__param_converter_vol[order]
                    for k_r in k_reacs_comp:
                        try: sim.setCompReacK(k_reacs_comp[k_r], k_r, value)
                        except: raise NotImplementedError("Negative constant of Reaction '%s'"%k_r)
                if (k_sreacs_patch): 
                    assert setorder
                    assert srtype
                    # Finally can convert to STEPS units
                    if srtype == '3D': value*= self.__param_converter_vol[order]
                    else: value*= self.__param_converter_area[order]
                    for k_sr in k_sreacs_patch:
                        try: sim.setPatchSReacK(k_sreacs_patch[k_sr], k_sr, value)
                        except: raise NotImplementedError("Negative constant of Surface Reaction '%s'"%k_sr)
        
        for p in poplist:
            self.__rules_rate.pop(p)
        
        for s in specs_update:
            self.__species[s][0] = specs_update[s]
        for c in comps_update:
            self.__comps[c][0] = comps_update[c]
        for p in patches_update:
            self.__patches[p][0] = patches_update[p]
        for p in params_update:
            self.__globalParameters[p][0] = params_update[p][0]
        
        variables.update(self.__time)
        variables.update(self.__species)
        #variables.update(self.__globalParameters)
        variables.update(self.__comps)
        variables.update(self.__patches)
        
        #####################################
        ######## "MATH REACTIONS" ###########
        #####################################
        # Don't work if the actual rate of the rection is zero- i.e. if there are no chemical species to turn into anything
        ## VOLUME REACTIONS
        
        # These are reactions with unexpected mathmeatics. 
        for mr in self.__math_reactions:
        
            # Might have a local paramter, so have to get the global ones then
            # possibly override with any local ones
            variables.update(self.__globalParameters)
            variables.update(mr[3])

            # Parameter can depend on anything and is actual rate/expected rate
            # of course expected rate is without parameter
            exp_rate = MLfunc(mr[2], variables)
            actual_rate =MLfunc(mr[1], variables)
                        
            # Occasional case where both are zero
            if exp_rate == 0.0: kconst = 0.0
            else: kconst = actual_rate/exp_rate
            if (kconst < 0.0): kconst = 0.0
            
            ord_tmp = mr[0].getOrder()
            
            # The kconst will actually be in units in terms of the species involved in the reaction. 
            # actual_rate has units (substance)/(time)
            # exp_rate units depends on the reaction:
            #   order 1 : (substance)
            #   order 2 : (substance)^2/(volume)
            
            # So actual rate/exp_rate will be units:
            #   order 1: /(time)
            #   order 2: (volume)/(substance).(time)
            # As we expect for the parameter. 
            #
            # Anyway, the substance unit should be in terms of the species involved in the model
            # so may well be in units other from model units. Convert depending on the order:
            kconst *= mr[4]
            
            # THe kconst is in SBML units, but we have nothing to convert to STEPS units
            # because it might not be in the global parameters. Need to do fancy tricks with the order
            # Value is in SBML units and should be converted to STEPS units
            
            kconst *= self.__param_converter_vol[ord_tmp]
            
            #kconst *= (math.pow(1000, ord_tmp-1))
            
            compid = mr[0].getCompID() 
            try: sim.setCompReacK(compid, mr[0].getName(),kconst)
            except: raise NotImplementedError("Negative rate of Reaction '%s'"%mr[0].getName())
        
        #####################################
        ##### "SURFACE MATH REACTIONS" ######
        #####################################
        ## SURFACE REACTIONS
        
        # These are reactions with unexpected mathmeatics. 
        for smr in self.__surface_math_reactions:
        
            # Might have a local paramter, so have to get the global ones then
            # possibly override with any local ones
            variables.update(self.__globalParameters)
            variables.update(smr[3])

            # Parameter can depend on anything and is actual rate/expected rate
            # of course expected rate is without parameter
            
            exp_rate = MLfunc(smr[2], variables)
            actual_rate =MLfunc(smr[1], variables)
            
            # Occasional case where both are zero
            if exp_rate == 0.0: kconst = 0.0
            else: kconst = actual_rate/exp_rate
            
            if (kconst < 0.0): kconst = 0.0
            
            
            ord_tmp = smr[0].getOrder()
            
            # The kconst will actually be in units in terms of the species involved in the reaction. 
            # actual_rate has units (substance)/(time)
            # exp_rate units depends on the reaction:
            #   order 1 : (substance)
            #   order 2 : (substance)^2/(volume)
            
            # So actual rate/exp_rate will be units:
            #   order 1: /(time)
            #   order 2: (volume)/(substance).(time)
            # As we expect for the parameter. 
            #
            # Anyway, the substance unit should be in terms of the species involved in the model
            # so may well be in units other from model units. Convert depending on the order:
            kconst*=smr[4]
            
            if (smr[0].getType()=='3D'): kconst *= self.__param_converter_vol[ord_tmp]
            else: kconst *= self.__param_converter_area[ord_tmp]
            
            patchid = smr[0].getPatchID() 
            try: sim.setPatchSReacK(patchid, smr[0].getName(), kconst)  
            except: raise NotImplementedError("Negative rate during of '%s'"%smr[0].getName())

        # Remove any local math parameters from variables, you never know what names they might have
        # and could for example shadow a species, so update all for safety:
        variables = {}
        variables.update(self.__time)
        variables.update(self.__species)
        variables.update(self.__globalParameters)
        variables.update(self.__comps)
        variables.update(self.__patches)
                
        
        #######################
        ####### EVENTS ########
        #######################   
        
        # NOTE: delay in event is assumed to be in model time units, as specified in SBML documentation
        for ev_trig in self.__evnts_trig:
            if (MLfunc(self.__evnts_trig[ev_trig], variables)):
                if (self.__evnts_flip[ev_trig] == False):
                    self.__evnts_flip[ev_trig] = True
                    # Store the time the event 'kicked-in', only if it's not already in the queue
                    if not (self.__evnts_fire.has_key(ev_trig)): 
                        self.__evnts_fire[ev_trig] = self.__time['time'][0]
                        self.__evnts_dl[ev_trig][1] = MLfunc(self.__evnts_dl[ev_trig][0], variables)*self.__time_units
                        # It might be necessary to evaluate the values at trigger time
                        if (self.__evnts_trigvals[ev_trig] == True):
                            self.__evnts_vals[ev_trig] = {}
                            for ass in self.__evnts_ass[ev_trig]:
                                var = ass[0]
                                self.__evnts_vals[ev_trig][var] = MLfunc(ass[1], variables)
            else: self.__evnts_flip[ev_trig] = False
        
        # Now execute event if delay has been passed: delay may be zero
        poplist = []
        # Storing a dictionary of any altered parameters, to be updated
        params_update = {}
        # Lets do it for comps and patches too so as to save some time
        comps_update={}
        patches_update={}
        specs_update={}
        
        for ev in self.__evnts_fire:
            stime = sim.getTime()
            fire_time = self.__evnts_fire[ev]
            delay_time = self.__evnts_dl[ev][1]
            if (stime >= fire_time + delay_time):
                trigvals = self.__evnts_trigvals[ev]
                for ass in self.__evnts_ass[ev]:
                    # Assignment variable could be a reaction constant or a species conc
                    var = ass[0]
                    if (trigvals): 
                        value = self.__evnts_vals[ev][var]
                    else: value = MLfunc(ass[1], variables)
                    if (var in self.__species):
                        # Convert to STEPS units
                        value*=self.__species[var][1]
                        # First update the new values
                        specs_update[var] = value
                        if self.__species_comps[var] in self.__comps:
                            if (self.__species_subst_units[var]):
                                try: sim.setCompAmount(self.__species_comps[var], var, value*self.__substance_units)
                                except: 
                                    raise NotImplementedError("Cannot set amount of species '%s'" \
                                    "in compartment '%s' during event assignment."%(var, self.__species_comps[var]))
                            else:
                                try: sim.setCompConc(self.__species_comps[var], var, value*(self.__substance_units/self.__volume_units)*1.0e-3)
                                except: 
                                    raise NotImplementedError("Cannot set concentration of species '%s'" \
                                    "in compartment '%s' during event assignment."%(var, self.__species_comps[var]))
                        elif self.__species_comps[var] in self.__patches:
                            if (self.__species_subst_units[var] == True):
                                try:
                                    sim.setPatchAmount(self.__species_comps[var], var, value*self.__substance_units)
                                except: 
                                    raise NotImplementedError("Cannot set amount of species '%s'" \
                                    "in patch '%s' during event assignment."%(var, self.__species_comps[var]))
                            else:
                                try:
                                    sim.setPatchAmount(self.__species_comps[var], var, value*(self.__substance_units/self.__area_units)*sim.getPatchArea(self.__species_comps[var]))
                                except: 
                                    raise NotImplementedError("Cannot set concentration of species '%s'" \
                                    "in patch '%s' during event assignment."%(var, self.__species_comps[var]))
                    elif (var in self.__comps):
                        # Convert to SBML units
                        value*=self.__comps[var][1]
                        comps_update[var] = value
                        # Convert to STEPS units
                        value*=self.__volume_units
                        try: sim.setCompVol(var, value)
                        except:
                            raise NotImplementedError("Failed to set compartment volume during event assignment (negative volume?).")
                    elif (var in self.__patches):
                        # Convert to SBML units
                        value*=self.__patches[var][1]
                        patches_update[var] = value
                        # Convert returned value in SBML units to STEPS units
                        value*=self.__area_units
                        try: sim.setPatchArea(var, value)
                        except:
                            raise NotImplementedError("Failed to set patch area during event assignment (negative area?).")
                    elif (var in self.__globalParameters):
                        # Convert to SBML units
                        if (self.__globalParameters[var][1]): value*=self.__globalParameters[var][1]
                        # Note: Not updating the global parameters here because they may appear in other rate rules
                        #       params[var] = value  
                        order = self.__glob_params_order[var]
                        
                        setorder = False
                        if (type(order) in (int, float)): setorder = True
                        
                        # Need to try and get all the reactions with this parameter and change- store in a list first
                        k_reacs_comp = {}
                        for r in self.__reactions:
                            if (r.getKName() == var): 
                                k_reacs_comp[r.getName()] = r.getCompID() 
                                ord_tmp = r.getOrder()
                                if setorder: assert ord_tmp == order 
                                else: 
                                    order = ord_tmp
                                    setorder = True
                        
                        srtype=''
                        k_sreacs_patch = {}
                        for sr in self.__surface_reactions:
                            if (sr.getKName() == var): 
                                k_sreacs_patch[sr.getName()] = sr.getPatchID()
                                ord_tmp = sr.getOrder()
                                if setorder: assert ord_tmp == order 
                                else: 
                                    order = ord_tmp
                                    setorder = True
                                type_temp = sr.getType()
                                if  srtype: assert type_temp == srtype
                                else:
                                    srtype = type_temp
                        
                        params_update[var] = [value]
                        
                        if (k_reacs_comp): 
                            assert setorder
                            value*= self.__param_converter_vol[order]
                            for k_r in k_reacs_comp:
                                try: sim.setCompReacK(k_reacs_comp[k_r], k_r, value)
                                except: raise NotImplementedError("Negative constant of Reaction '%s'"%k_r)

                        if (k_sreacs_patch): 
                            assert setorder
                            assert srtype
                            # Finally can convert to STEPS units
                            if srtype == '3D': value*= self.__param_converter_vol[order]
                            else: value*= self.__param_converter_area[order]
                            for k_sr in k_sreacs_patch:
                                try: sim.setPatchSReacK(k_sreacs_patch[k_sr], k_sr, value)   
                                except: raise NotImplementedError("Negative constant of Surface Reaction '%s'"%k_sr)
                    else: 
                        raise NotImplementedError("Currently only assignment rules referencing species" \
                            ", compartments and parameters are supported.")
                
                poplist.append(ev) 
                self.__evnts_vals[ev] = {}
        
        # Get rid of fired events
        for pevent in poplist:
            self.__evnts_fire.pop(pevent)
        
        for s in specs_update:
            self.__species[s][0] = specs_update[s]
        for c in comps_update:
            self.__comps[c][0] = comps_update[c]
        for p in patches_update:
            self.__patches[p][0] = patches_update[p]
        for p in params_update:
            self.__globalParameters[p][0] = params_update[p][0]
        
    ################################################################################################

####################################################################################################
####################################################################################################

class Reaction(object):
    def __init__(self):
        self.__name = ""
        self.__reacts = []
        self.__lhs = []
        self.__prods = []
        self.__rhs = []
        self.__kName = ""
        self.__vsys= ""
        self.__compid = ""
        self.__kValue = 0
        self.__unitsConv = 1

    def getKName(self):
        return self.__kName

    def getKValue(self):
        return self.__kValue

    def setKName(self, value):
        self.__kName = value

    def setKValue(self, value):
        self.__kValue = value

    def setVsys(self, value):
        self.__vsys = value
    
    def getVsys(self):
        return self.__vsys

    def setCompID(self, value):
        self.__compid = value
        
    def getCompID(self):
        return self.__compid

    def getName(self):
        return self.__name
    
    def setName(self, value):
        self.__name = value
    
    def getLhs(self):
        return self.__lhs
    
    def getOrder(self):
        return len(self.__lhs)
    
    def getRhs(self):
        return self.__rhs
    
    def setLhs(self, value):
        self.__lhs = value
    
    def setRhs(self, value):
        self.__rhs = value
            
    def delLhs(self):
        del self.__lhs
    
    def delRhs(self):
        del self.__rhs
    
    def getReacts(self):
        return self.__reacts
    
    def getProds(self):
        return self.__prods
    
    def setReacts(self, value):
        self.__reacts = value
    
    def setProds(self, value):
        self.__prods = value
    
    def delReacts(self):
        del self.__reacts
    
    def delProds(self):
        del self.__prods
    
    reacts = property(getReacts, setReacts, delReacts, "The reactants of the reaction")
    
    prods = property(getProds, setProds, delProds, "The products of the reaction")
    
    lhs = property(getLhs, setLhs, delLhs, "Left side STEPS")
    
    rhs = property(getRhs, setRhs, delRhs, "Right side STEPS")
    
    name = property(getName, setName, None, "Unique name of the reaction")
    
    kName = property(getKName, setKName, None, "The name of the costant")
    
    kValue = property(getKValue, setKValue, None, "The Value of the costant")
        
    vsys = property(getVsys, setVsys, None, "The volume system string")

####################################################################################################

class SurfaceReaction(object):
    def __init__(self):
        self.__name = ""
        self.__reacts = []
        self.__olhs = []
        self.__slhs = []
        self.__ilhs = []
        self.__prods = []
        self.__orhs = []
        self.__srhs = []
        self.__irhs = []        
        self.__kName = ""
        self.__ssys = ""
        self.__patchid = ""
        self.__type = ""
        self.__kValue = 0
        self.__unitsConv = 1

    def getKName(self):
        return self.__kName

    def getKValue(self):
        return self.__kValue

    def setKName(self, value):
        self.__kName = value

    def setKValue(self, value):
        self.__kValue = value

    def setSsys(self, value):
        self.__ssys = value
    
    def getSsys(self):
        return self.__ssys

    def setPatchID(self, value):
        self.__patchid = value
    
    def getPatchID(self):
        return self.__patchid
    
    def setType(self, value):
        self.__type = value
        
    def getType(self):
        return self.__type
    
    def getName(self):
        return self.__name
    
    def setName(self, value):
        self.__name = value
    
    def getOLhs(self):
        return self.__olhs

    def getSLhs(self):
        return self.__slhs

    def getILhs(self):
        return self.__ilhs
    
    def getORhs(self):
        return self.__orhs

    def getSRhs(self):
        return self.__srhs

    def getIRhs(self):
        return self.__irhs
    
    def getOrder(self):
        return len(self.__slhs)+len(self.__ilhs)+len(self.__olhs)
    
    def setOLhs(self, value):
        self.__olhs = value

    def setSLhs(self, value):
        self.__slhs = value

    def setILhs(self, value):
        self.__ilhs = value
    
    def setORhs(self, value):
        self.__orhs = value

    def setSRhs(self, value):
        self.__srhs = value

    def setIRhs(self, value):
        self.__irhs = value
            
    def delOLhs(self):
        del self.__olhs

    def delSLhs(self):
        del self.__slhs

    def delILhs(self):
        del self.__ilhs
    
    def delORhs(self):
        del self.__orhs

    def delSRhs(self):
        del self.__srhs

    def delIRhs(self):
        del self.__irhs
    
    def getReacts(self):
        return self.__reacts
    
    def getProds(self):
        return self.__prods
    
    def setReacts(self, value):
        self.__reacts = value
    
    def setProds(self, value):
        self.__prods = value
    
    def delReacts(self):
        del self.__reacts
    
    def delProds(self):
        del self.__prods
    
    def flip(self):
        ilhs_temp = self.__ilhs
        self.__ilhs = self.__olhs
        self.__olhs = ilhs_temp
        
        irhs_temp = self.__irhs
        self.__irhs = self.__orhs
        self.__orhs = irhs_temp
    
    reacts = property(getReacts, setReacts, delReacts, "The reactants of the reaction")
    
    prods = property(getProds, setProds, delProds, "The products of the reaction")
    
    olhs = property(getOLhs, setOLhs, delOLhs, "Left side STEPS")
    slhs = property(getSLhs, setSLhs, delSLhs, "Left side STEPS")
    ilhs = property(getILhs, setILhs, delILhs, "Left side STEPS")
    
    orhs = property(getORhs, setORhs, delORhs, "Right side STEPS")
    srhs = property(getSRhs, setSRhs, delSRhs, "Right side STEPS")
    irhs = property(getIRhs, setIRhs, delIRhs, "Right side STEPS")
    
    name = property(getName, setName, None, "Unique name of the reaction")
    
    kName = property(getKName, setKName, None, "The name of the costant")
    
    kValue = property(getKValue, setKValue, None, "The Value of the costant")
        
    ssys = property(getSsys, setSsys, None, "The surface system string")

####################################################################################################    
####################################################################################################

