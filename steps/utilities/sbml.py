# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import libsbml
import sys
import os
import math

import steps.model as smodel
import steps.geom as sgeom

####################################################################################################
#                                           Constants                                              #
####################################################################################################

# Default volume of one femtolitre to avoid the situation of running a stochastic simulation
# on the scale of litres
VOL_DEFAULT = 1.0e-18

# Source: http://physics.nist.gov/cgi-bin/cuu/Value?na 10/11/2010
AVOGADRO = 6.02214179e23

####################################################################################################
#                                       Auxilliary functions                                       #
####################################################################################################

def is_num(n):
    """
    Returns whether argument is a number or not
    """
    return isinstance(n, (float, int, long, complex))

def _float_approx_equal(x, y, tol=1e-16, rel=1e-7):
    """
    Return whether two floats are approximately equal
    """
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
            Note::  factor is the factor for converting between STEPS and SBML units: 
                    values will be in STEPS units so the values should be DIVIDED by factor to get 
                    SBML units in expression.
    """
    
    # This is why we need to do the units conversion, the number is in SBML units
    if (is_num(level)): return level
    
    if (isinstance(level, (str))):
        # Dividing by factor to keep in SBML units
        if (level in refs): return refs[level][0]/refs[level][1]
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
        return MLfunc(level[1][0], refs) / MLfunc(level[1][1], refs)    
    elif (chrct == 'power'):
        return math.pow(MLfunc(level[1][0], refs), MLfunc(level[1][1], refs))
    elif (chrct == 'root'):
        return math.pow(MLfunc(level[1][0], refs), (1.0/MLfunc(level[1][1], refs)))
    elif (chrct == 'lt'):
        return MLfunc(level[1][0], refs) < MLfunc(level[1][1], refs)
    elif (chrct == 'leq'):
        return MLfunc(level[1][0], refs) <= MLfunc(level[1][1], refs)
    elif (chrct == 'gt'):
        return MLfunc(level[1][0], refs) > MLfunc(level[1][1], refs)
    elif (chrct == 'geq'):
        return MLfunc(level[1][0], refs) >= MLfunc(level[1][1], refs)
    elif (chrct == 'eq'):
        return MLfunc(level[1][0], refs) == MLfunc(level[1][1], refs)
    elif (chrct == 'neq'):
        return MLfunc(level[1][0], refs) != MLfunc(level[1][1], refs)
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
            if (MLfunc(level[1][(i*2)+1], refs)): return MLfunc(level[1][i*2], refs)
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
            return refs[level][0]/refs[level][1]
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

# Converts a mathMl object to a more readable Python layered list. Easiest to explain with examples:
#
# 1 + 2  ->  ['+', [1,2]]
#
# 1 + (2 * 3)  ->  ['+', [1, ['*', [2,3] ] ] ]
#
# etc etc

# Replace_dict is an optional dictionary of function arguments to real arguments, i.e model variables
# This is useful for lambda functions

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
    concentration are molar units, i.e. bsaed on litres not cubic metres, but this appears to be 
    the assumption in SBML too.
    """
    return (math.pow((mult *math.pow(10, sca)), exp))

####################################################################################################
####################################################################################################

class Interface(object):
    """
    The interface to the SBML importer.
    """
    
    def __init__(self, filename, defvolunits_litre = True):
        """
        Construction::
        
            iSbml = steps.utilities.sbml.Interface(sbmlFile, defvolunits_litre = True)
            
        Construct an SBML interface object, and optionally declare whether default volume units are
        litres rather than cubic meters. 
        
        Arguments:
            * string sbmlFile
            * bool defvolunits_litre (default = True)
        """
        self.__reader = libsbml.SBMLReader()
        self.__document = self.__reader.readSBML(filename)
        self.__model = self.__document.getModel()
        if self.__document.getNumErrors() > 0 :
            pm = self.__document.getError(0)
            if pm.getErrorId() == libsbml.XMLFileUnreadable:
                raise IOError("Sbml File '%s' not found or unreadable in dir %s" %(filename, os.getcwd())) 
        
        # Sets 
        # __globalParamters: a dictionary with elements list len<2>; {'parameter_id': [value, factor]
        #       Value is, as will be standard, in STEPS units
        #       Factor * SBML value => value in STEPS (s.i.) units ; STEPS value/factor converts to SBML units
        # __glob_params_order: a dictionary mapping: {'parameter_id': order}
        #       Order may be 'unknown', if units are not specified or parameter is not a reaction constant
        #
        self.__globalParameters, self.__glob_params_order = self._parseGlobalParameters()
        
        # Basic function parsing. May be used in lots of places in the model.
        self.__function_defs = self._parseFunctionDefs()
        
        # Get the time structure, with units, and the volume units from the model in one fell swoop
        # NOTE: I am assuming the SBML documentation is wrong here and that volume units has a bearing
        # on reactions other than order 1
        # The simulation time, with conversion units to STEPS 
        self.__time, self.__volume_units, self.__substance_units = self._parseModelUnits(defvolunits_litre) 
        
        # Returns a dictionaray with elements list len<2>; {'compr_id': [value, factor]
        # Value is, as will be standard, in STEPS units (here m^3)
        # Factor * SBML value => value in STEPS (s.i.) units ; STEPS value/factor converts to SBML units
        self.__comps = self._parseComps()
        
        # Create the geometry and model container objects in STEPS
        self.__geom = sgeom.Geom()
        self.__mdl = smodel.Model()
        
        # Create a volume system for each compartment, and store dictionary- {comp_id:volsys_id}
        # volsys_id will just be comp_id + 'volsys'
        self.__comp_volsys = {}
        # Simply creates the compartments in STEPS and adds 'volsys' volume system.
        self._setupGeom()
        
        # A lot of information about the species in the model, so stored in separate objects:
        #
        # self.__species is a dictionary mapping: {'species_id': [value, factor]}
        #   however the value can be concentration or amount: 
        #       If concentration; STEPS units are MOLAR and factor*SBML_value will give Molar concentration
        #       If amount; STEPS units are MOLE and factor * SBML_value will give amount in Moles
        # 
        # self.__species_amount_flag: Used to be important but now kind of defunct:
        #   is True if quanitity is Amount, False if Concentration
        #
        # self.__species_comps: maps species id to the id of its compartment in SBML. 
        #   NOTE: species can only belong to one compartment in SBML, not the case in STEPS
        #
        # self.__species_const_bc: maps id to tuple of two booleans: (constant, boundary_condition)
        #
        self.__species, self.__species_amount_flag, self.__species_subst_units,\
            self.__species_comps, self.__species_const_bc, = self._parseSpecies()
        
        # Create the chemical species in STEPS
        # NOTE: This must come after _setupGeom because it needs acces to the volume systems
        self._setupModel()
        
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
        # Reaction object contains method getUnitsConv which returns a factor by which to convert
        #   reaction constant from SBML units to STEPS units
        self.__reactions = self._parseReactions()
        
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
        
    def _parseModelUnits(self, vol_units_litre):
        """
        Attempt to read model units for time, volume and substance, though often absent. 
        """
        
        ret_time = {}
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
                    # TODO: more checks the units aren't silly?
                    ret_time['time'] = [0.0, unitsDef_to_STEPS(multiplier, scale, expo)]
                    ret_time['t'] = [0.0, unitsDef_to_STEPS(multiplier, scale, expo)] 
                else:
                    raise NotImplementedError("Time base units are not seconds!")
            else:
                "WARNING: Failed to read model units for time. Default units (second) will be used."
                ret_time['time'] = [0.0, 1]
                ret_time['t'] = [0.0, 1]
        else:
            print "WARNING: Failed to read model units for time. Default units (second) will be used."
            ret_time['time'] = [0.0, 1]
            ret_time['t'] = [0.0, 1]
        
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
                # TODO: more checks the units aren't silly?
                if (unit.isLitre()):
                    if (expo != 1): raise NotImplementedError("Volume units in model are not volume units!")
                    ret_vol_units = unitsDef_to_STEPS(multiplier, scale, expo)*0.001
                elif (unit.isMetre()):
                    if (expo != 3): raise NotImplementedError("Volume units in model are not volume units!")
                    ret_vol_units = unitsDef_to_STEPS(multiplier, scale, expo)
                else: raise NotImplementedError("Volume units in model are not volume units!")
            else:
                if (vol_units_litre):
                    print "WARNING: Failed to read model units for volume. Litre units will be used as specified."
                    ret_vol_units =  0.001
                else:
                    print "WARNING: Failed to read model units for volume. Default units (m^3) will be used."
                    ret_vol_units =  1
        else:
            if (vol_units_litre):
                print "WARNING: Failed to read model units for volume. Litre units will be used as specified."
                ret_vol_units =  0.001
            else:
                print "WARNING: Failed to read model units for volume. Default units (m^3) will be used."
                ret_vol_units =  1
        
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
                # TODO: more checks the units aren't silly?
                if (unit.isMole()):
                    ret_subs_units = unitsDef_to_STEPS(multiplier, scale, expo)
                elif (unit.isDimensionless() or unit.isItem()):
                    ret_subs_units = unitsDef_to_STEPS(multiplier, scale, expo)/AVOGADRO
                else: raise NotImplementedError("Substance units in model are not supported units.")
            else:
                print "WARNING: Failed to read model units for substance. Default units (mole) will be used."
                ret_subs_units =  1
        else:
            print "WARNING: Failed to read model units for substance. Default units (mole) will be used."
            ret_subs_units =  1
        
        return ret_time, ret_vol_units, ret_subs_units
    
    ################################################################################################
    
    def _parseComps(self):
        """
        Import the compartments.
        """
        
        ListOfComp = self.__model.getListOfCompartments()
        
        comps = {}
        for comp in ListOfComp:
        
            idComp = str(comp.getId())
            
            # It appears now in level 3 that a user doesn;t even need to set the dimensions of a compartment. 
            # How I wish they would consider developers and not just modellers!!
            if (comp.isSetSpatialDimensions()):
                if (comp.getSpatialDimensions() != 3):
                    raise NotImplementedError("Compartments must be 3D.", \
                    "Compartment'%s' in this model %d dimensions." %(idComp, comp.getSpatialDimensions()))
            else:
                # What choice do we have? Assume 3d?? Oh well, will get checked if units are set.
                # If not just have to give some default volume
                pass
            
            # Need a bool to keep track if we have already found sufficient volume for comp
            added = False
            
            volComp = comp.getVolume()
            if (volComp):
                # Find the units:
                unitsComp = comp.getUnits()
                # Try to get from the model
                if not unitsComp: unitsComp = self.__model.getVolumeUnits()
                if not unitsComp: unitsComp = 'volume'
                if (unitsComp):
                    # Using a for loop that can be exited with a break, 
                    # for want of a better alternative
                    for j in [0]:
                        unitdef = self.__model.getUnitDefinition(unitsComp)
                        if not (unitdef):
                            print "WARNING: Failed to read units for compartment '%s'. " \
                            "Model volume units will be assumed."%idComp
                            break
                        if (unitdef.isVariantOfVolume() != True):
                            raise NotImplementedError("Compartment '%s' dimensions not 3D."%idComp)
                        
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
                            volComp *= unitsDef_to_STEPS(multiplier, scale, 3)
                            # self.__volComp[idComp] = volComp
                            # self.__unitsComp[idComp] = unitsDef_to_STEPS(multiplier, scale, 3)
                            comps[idComp] = [volComp, unitsDef_to_STEPS(multiplier, scale, 3)]
                            added = True
                            break
                        elif (unit.isLitre()):
                            if (unit.getExponent() != 1):
                                raise NotImplementedError("Compartment '%s' dimensions not 3D."%idComp)
                            scale = unit.getScale()
                            multiplier = unit.getMultiplier()
                            # Native units of steps are m^3
                            volComp *= (unitsDef_to_STEPS(multiplier, scale, 1)*1.0e-3)
                            comps[idComp] = [volComp, (unitsDef_to_STEPS(multiplier, scale, 1)*1.0e-3)]
                            added = True
                            break
                        else: 
                            print "WARNING: Failed to read units for compartment '%s'. " \
                            "Model volume units will be assumed."%idComp
                            break
                        break
                
                    # If we got here we breaked. 
            else:
                print "WARNING: No volume specified for compartment '%s'. " \
                        "Default value (%s) will be used"%(idComp, str(VOL_DEFAULT))
                volComp = VOL_DEFAULT  
                comps[idComp] = [volComp, self.__volume_units]
                added = True
                continue
            
            # Special case if size == 1 without units, a common occurunce in SBML
            if (comp.getSize() == 1.0 and not comp.getUnits()):
                print "WARNING: Compartment '%s' size = 1.0 wih no units. " \
                "Default volume value (%sm^3) will be used."%(idComp, str(VOL_DEFAULT))
                volComp = VOL_DEFAULT  
                comps[idComp] = [volComp , self.__volume_units]
                added = True
                continue  
            
            # Special case if size == 1 without default units (litre), also a common occurance
            if (comp.getSize() == 1.0 and idComp not in comps):
                print "WARNING: Compartment '%s' size = 1.0 wih default units. " \
                "Default volume value (%sm^3) will be used."%(idComp, str(VOL_DEFAULT))
                volComp = VOL_DEFAULT  
                comps[idComp] = [volComp , self.__volume_units]
                added = True
                continue   
            
            if not added:
                # If we got here we have a volume, which is not equal to 1.0, 
                # but couldn't read the units. Assume model units
                comps[idComp] = [volComp*self.__volume_units, self.__volume_units]
                added = True
                continue
        
        return comps
  
    ################################################################################################

    def _setupGeom(self):
        """ 
        Create the compartments in STEPS.
        """
        
        for c in self.__comps:
            comp = sgeom.Comp(c, self.__geom, vol = self.__comps[c][0])
            vsysid = c+'_volsys'
            comp.addVolsys(vsysid)
            # Doing the volsys here just to keep it in one place
            vsys = smodel.Volsys(vsysid, self.__mdl)
            self.__comp_volsys[c] = vsysid

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
        species_const_bcs = {}
        
        factor = 1
        
        for specie in ListOfSpecies:
            # Any volume units should match the compartment units
            compId = specie.getCompartment()
            species_comps[str(specie.getId())] = compId
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
            if (specie.isSetHasOnlySubstanceUnits()): hasAmountUnits = specie.getHasOnlySubstanceUnits()
            else: hasAmountUnits = isAmount
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
                        print "WARNING: Failed to read units of substance " \
                        "for Species '%s'. Default units (mole) will be assumed."%specie.getId()
                        break
                    # Note: All units here should be of substance, i.e. dimensionless
                    # Units for the volume (if concentration) are taken from 
                    # the units for volume of the compartment
                    units = unitdef.getListOfUnits()
                    if (len(units) != 1):
                        print "WARNING: Failed to read units of substance " \
                        "for Species '%s'. Default units (mole) will be assumed."%specie.getId()
                        break
                    unit = units[0]
                    if (unit.getExponent() != 1):
                        print "WARNING: Failed to read units of substance " \
                        "for Species '%s'. Default units (mole) will be assumed."%specie.getId()
                        break                                
                    multiplier = unit.getMultiplier()
                    scale = unit.getScale()
                    factor = unitsDef_to_STEPS(multiplier, scale, 1)
                    value  = value*factor
                    # Convert to mole
                    if (unit.isMole() != True): 
                        value /= AVOGADRO
                        factor /= AVOGADRO
                    break
            else: 
                print "WARNING: Failed to read units of substance " \
                "for Species '%s'. Default units (mole) will be assumed."%specie.getId()
            
            # Confusing, but value could still be an amount or conc, but the amount part is converted to 
            # Mole. If a conc we still need to convert the volume part. 
            # Another complication is that we are going to set to conc or amount by the hasSubstanceUnits flag. 
            if (isAmount):
                if (hasAmountUnits): 
                    # Easy, we have an amount already converted to mole and we want to store an amount.
                    species[str(specie.getId())] = [value, factor] 
                else: 
                    # We have an amount in mole, but we want to store the concentration. Simply 
                    # divide by the comp volume (converted from^3 to litre) to get the conc in mole/litre.
                    # The units will be factor for the amount divided by the comp units, because it SBML conc -> STEPS conc
                    conc = value/(self.__comps[compId][0]*1.0e3)
                    species[str(specie.getId())] = [conc, factor/(self.__comps[compId][1]*1.0e3)]
            else:
                if (hasAmountUnits):
                    # We've got a conc (the volume part still in SBML units), but have to convert to amount for expressions
                    # We have to multiply by the volume of the compartment converted to SBML units
                    amount = value*((self.__comps[compId][0]*1.0e3)/self.__comps[compId][1])
                    # Units are just going to be factor
                    species[str(specie.getId())] = [amount, factor]
                else:
                    # We've got a conc and want to store a conc. Just need to convert the units for comp.
                    # Concentration units for compartment may be other than s.i.
                    # e.g. if comp units are microns, concentration should be {substance units}/micron^3
                    # Convert to MOLAR concentration, so we need litres
                    conc = value/(self.__comps[compId][1]*1.0e3)
                    species[str(specie.getId())] = [conc, factor/(self.__comps[compId][1]*1.0e3)]
            
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
    
    def _setupModel(self):
        """ 
        Setup the model in STEPS.
        """
        
        for specie in self.__species:
            mol = smodel.Spec(specie, self.__mdl)
            # This is a bit of a hack, but a few occasions when a species does not appear in
            # reactions, or stoichiometry changed becase of BCs so doesn't appear in reactions, 
            # and therefore not added to comp by STEPS.
            comp_id= self.__species_comps[specie]
            reac = smodel.Reac(specie+'reac', self.__mdl.getVolsys(self.__comp_volsys[comp_id]), lhs = [self.__mdl.getSpec(specie)], \
                rhs = [self.__mdl.getSpec(specie)], kcst = 0.0)
        
    ################################################################################################
            
    def _parseReactions(self):
        """
        Import the reactions.
        """
        listOfReactions = self.__model.getListOfReactions()
        reactions = []
        for reac in listOfReactions:
            r_f = Reaction()
            reversible = reac.getReversible()
            if (reversible): 
                        r_f.setName(str(reac.getId())+"_f") 
                        r_b = Reaction()
            else : r_f.setName(str(reac.getId())) 
            if (reversible): r_b.setName(str(reac.getId())+"_b")
            
            # Reactants: list of ids of the reactants, separated for reversible reactions
            reacts_f = []
            if (reversible): reacts_b = []
            #lhsList = [] # Left side STEPS
            
            # Products: list of ids of the products, separated for reversible reactions  
            prods_f = []
            if (reversible): prods_b = []
            #rhsList = [] # Right side STEPS
            
            products = reac.getListOfProducts()            
            reactants = reac.getListOfReactants()
            
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
                    sto = self.__spec_refs[spec_ref]
                # stoichiometry must be a whole number (allowing a small error)
                if (sto % 1 > 0.001) : 
                    raise NotImplementedError("Partial stoichiometry (%f) in reaction '%s'."%(sto, reac.getId()))
                sto = int(round(sto))
                # More than one molecule of each reactant species may be present in the reaction
                for j in range(sto):
                    if(reactant.getSpecies() != "Empty"):
                        reacts_f.append(reactant.getSpecies())
                        if (reversible): prods_b.append(reactant.getSpecies())
                        # Adding the mol Obj from STEPS
                        #lhsList.append(self.__mdl.getSpec(reactant.getSpecies()))  
            
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
                    # Hmm. Probably a species reference. 
                    # Bizarrely, getId returns the species reference 
                    spec_ref = product.getId()
                    sto = self.__spec_refs[spec_ref]
                if (sto % 1 ) : 
                    raise NotImplementedError("Partial stoichiometry (%f) in reaction '%s'."%(sto, reac.getId()))
                sto = int(round(sto))
                for j in xrange(sto):
                    if(product.getSpecies() != "Empty"):
                        prods_f.append(product.getSpecies())
                        if (reversible): reacts_b.append(product.getSpecies())
                        #rhsList.append(self.__mdl.getSpec(product.getSpecies()))
            
            # Modifiers
            mods = reac.getListOfModifiers()
            for mod in mods:
                if(mod.getSpecies() != "Empty"):
                    # Adding the modifiers to lhs and rhs since we don't want any to be destroyed
                    reacts_f.append(mod.getSpecies())
                    if (reversible): prods_b.append(mod.getSpecies())
                    #lhsList.append(self.__mdl.getSpec(mod.getSpecies())) 
                    prods_f.append(mod.getSpecies())
                    if (reversible): reacts_b.append(mod.getSpecies())
                    #rhsList.append(self.__mdl.getSpec(mod.getSpecies()))
            
            # A CHECK AS TO WHETHER SPECIES ARE IN DIFFERENT COMPS OR NOT, currently not supported
            # This could possibly be built upon in the future to have proper surface reactions - a lot of work??
            reac_comp = ''
            for r in reacts_f:
                reac_comp_next = self.__species_comps[r]
                if (reac_comp ==''): reac_comp = self.__species_comps[r]
                if (reac_comp != reac_comp_next):
                    raise NotImplementedError("Reaction '%s': reactants and products do not belong to one compartment."%reac.getId())
                reac_comp = reac_comp_next
            for r in prods_f:
                reac_comp_next = self.__species_comps[r]
                if (reac_comp ==''): reac_comp = self.__species_comps[r]
                if (reac_comp != reac_comp_next):
                    raise NotImplementedError("Reaction '%s': reactants and products do not belong to one compartment."%reac.getId())
                reac_comp = reac_comp_next    
            assert (reac_comp != '')
            
            # Tell the reaction which volsys it belongs to: compid+'_volsys'
            r_f.setVsys(reac_comp+'_volsys')
            if (reversible): 
                r_b.setVsys(reac_comp+'_volsys')
            
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
            # In the following I am assuming units of (molar)^n/s if not defined so help me god
            for p in parameters:
                p_id = str(p.getId())
                p_value = p.getValue()
                p_units = p.getUnits()
                if (not p_value): 
                    if p_id in self.__globalParameters:
                        p_value, p_factor  = self.__globalParameters[p_id]     
                        p_order = self.__glob_params_order[p_id]
                        if (p_order == 'not reac'): 
                            raise NotImplementedError("Reaction '%s' parameter '%s' not the right units for a reaction parameter."%(reac.getId()), p_id)
                        elif (type(p_order) in (int, float)):
                            # We must convert to STEPS (Molar) units here if we have a known parameter
                            fac = math.pow(0.001, p_order-1)
                            p_value = p_value/fac
                            p_factor = p_factor/fac                        
                        
                        params.append([p_id, p_order, p_value, p_factor])
                        # Reaction params already converted to STEPS, Molar units. 
                        # Will ckeck later is this param is not a reaction param
                        continue  
                else:
                    """
                    # Give it a local only name incase it is a shadower
                    p_id+=reac.getId()
                    """
                    pass
                # For local parameters we are looking specifically for things involved in a reaction. If 
                # they are other things (checked in the units if specified) we simply don't add them.
                if (p_units): 
                    # Using a for loop that can be exited with a break, 
                    # for want of a better alternative
                    for j in [0]:
                        unitdef = self.__model.getUnitDefinition(p_units) 
                        if not (unitdef):
                            print "WARNING: Local parameter '%s' has unknown unit definition " \
                            "and will cause error if used in rate expression."%p_id
                            break
                        units = unitdef.getListOfUnits()  
                        
                        # Units for a reaction should be (conc)^-(n-1).(time)^-1
                        # Where n gives the order of the reaction
                        # NOTE: If a unit is bad it is imply not added to the list of parameters for the search.
                        gotTime = False
                        gotLitre = False
                        gotMetre = False
                        gotMole = False
                        gotSecond = False
                        order = -1
                        badunit = False
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
                            elif (u.isLitre()):
                                if (gotLitre or gotMetre): badunit = True
                                elif (exp%1 != 0): badunit = True
                                elif (order != exp+1) and (order != -1): badunit = True # If we already have an order from substance, break if different
                                order = exp+1
                                if (order < 0): badunit = True
                                gotLitre = True
                            elif (u.isMetre()): 
                                if (gotLitre or gotMetre): badunit = True
                                elif (exp%3 != 0): badunit = True
                                elif (order != (exp/3)+1) and (order != -1): badunit = True
                                order  = (exp/3)+1
                                if (order < 0): badunit = True
                                # Convert to litres
                                p_factor *= math.pow(1000, order-1)
                                p_value *= math.pow(1000, order-1)
                                gotMetre = True
                            elif (u.isMole()):
                                if (gotMole): badunit = True
                                elif (exp%1 != 0): badunit  = True
                                elif (order != -(exp-1)) and (order != -1): badunit = True
                                order = -(exp-1)
                                if (order < 0): badunit = True
                                gotMole = True
                            elif (u.isSecond()):
                                if (gotSecond): badunit = True
                                elif (exp != -1):  badunit = True
                                gotSecond = True
                            else: 
                                badunit = True
                            # Breaking here means the local parameter is not a reaction parameter and will not be added to params
                            if(badunit): 
                                print "WARNING: Local parameter '%s' does not have correct reaction " \
                                "parameter units and will cause error if in rate expression."%p_id
                                bad_params.append(p_id)    
                                break
                        # Set first order if so
                        if (gotSecond and order == -1): order = 1    
                        # Some sanity checks: break if true and don't add param to params (print warning??)
                        if (badunit or order < 0 or gotSecond == False or p_value < 0): 
                            print "WARNING: Local parameter '%s' does not have correct reaction " \
                            "parameter units and will cause error if in rate expression."%p_id
                            bad_params.append(p_id)
                            break
                        if order != 1: 
                            if not (gotMole and (gotLitre or gotMetre)): 
                                print "WARNING: Local parameter '%s' does not have correct reaction " \
                                "parameter units and will cause error if in rate expression."%p_id
                                bad_params.append(p_id)                            
                                break
                        # NOTE: AS this is specifically a reaction parameter we have already converted to STEPS (molar) units
                        params.append([p_id, order, p_value, p_factor])
                        break
                    
                else: # No units, assume default, whatever they are
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
            
            # First convert our local parameters to a dictionary for the maths function
            params_temp = {}
            for p in params:
                params_temp[p[0]] = [p[2], p[3]]
            
            all_variables = {}
            params_variables = {}
            # Note: params_temp should come after global parameters because it may be a shadower
            all_variables.update(self.__globalParameters)
            all_variables.update(params_temp)
            # all_variables.update(self.__species)
            # Not using species from model now because there may be zero concentrations messing up
            # the rate comparison
            species_temp = {}
            for s in self.__species:
                species_temp[s] = [None, None]
                species_temp[s][0] = self.__species[s][0]
                species_temp[s][1] = self.__species[s][1]
                
                # Set the value to some non-zero arbitrary number
                if (species_temp[s][0] == 0.0 or _float_approx_equal(species_temp[s][0], 1.0)):
                    species_temp[s][0] = 2.5
            
            all_variables.update(species_temp)
            all_variables.update(self.__comps)
            # Note: params_temp should come after global parameters because it may be a shadower
            params_variables.update(self.__globalParameters)
            params_variables.update(params_temp)
            
            params_id = []
            rate = rate_func(rate_list, all_variables, params_variables, params_id, bad_params)
            # Separate out the parameters
            if not reversible:
                if len(params_id) != 1:
                    raise NotImplementedError("Reaction rate maths in reaction '%s' is not of supported form."%reac.getId())
            else:
                if len(params_id) != 2:
                    raise NotImplementedError("Reaction rate maths in reaction '%s' is not of supported form."%reac.getId())
            
            # Update the params list
            for p_id in params_id:
                if (p_id in bad_params):
                    raise NotImplementedError("Reaction rate maths in reaction '%s' includes" \
                        "parameter '%s' with unexpected units."%(reac.getId(), p_id))
                if (p_id in self.__globalParameters):
                    gotP = False
                    for p in params:
                        if (p[0] == p_id): gotP = True 
                    if not gotP:
                        p_value, p_factor  = self.__globalParameters[p_id]    
                        p_order = self.__glob_params_order[p_id] 
                        if (p_order == 'not reac'): 
                            raise NotImplementedError("Reaction '%s' parameter '%s' not the right units for a reaction parameter."%(reac.getId()), p_id)
                        elif (type(p_order) in (int, float)):
                            # We must convert to STEPS units here if we have a known parameter
                            fac = math.pow(0.001, p_order-1)
                            p_value = p_value/fac
                            p_factor = p_factor/fac    
                        params.append([p_id, p_order, p_value, p_factor])
            
            # OK, now we are here we may have some parameters that have not specified units, the flag is p_order = 'unknown'
            # otherwise they have been converted already
            
            # Right, hopefully that worked (error will have been thrown if not). 
            # Now the tricky part of constructing the expected list
            
            # Lets treat forward and backward separately (for reversible reactions), but create a holder
            # for a reverse reaction regardless (initialised to 0, ignored for non-reversible reactions)
            correct_rate_list_1 = ['-', [None, 0.0]]
            
            if (reversible): correct_rate_list_2 = ['-', [None, None]]
            
            # just a reminder reacts_f is a list of the reactant species by string id for the forward reaction
            # and reacts_b is a list of reactant species for the backward reaction (if there is one)
            # reac_comp is the compartment id
            
            # The 'forward' reaction first: the compartment, parameter0, and the species
            for_vars_1 = [reac_comp, params_id[0]]+reacts_f
            correct_rate_list_1[1][0] = make_rate_list(for_vars_1, self.__species_subst_units, reac_comp)            
            
            if (reversible):
                back_vars_1 = [reac_comp, params_id[1]]+reacts_b
                correct_rate_list_1[1][1] = make_rate_list(back_vars_1, self.__species_subst_units, reac_comp)
            
            if (reversible):
                for_vars_2 = [reac_comp, params_id[1]]+reacts_f
                correct_rate_list_2[1][0] = make_rate_list(for_vars_2, self.__species_subst_units, reac_comp)
                back_vars_2 = [reac_comp, params_id[0]]+reacts_b
                correct_rate_list_2[1][1] = make_rate_list(back_vars_2, self.__species_subst_units, reac_comp)
            
            variables={}
            #variables.update(self.__species)
            # Replacing with the arbitrary non-zero initial values
            variables.update(species_temp)
            
            variables.update(self.__globalParameters)
            variables.update(self.__comps)
            variables.update(params_temp)
            correct_rate_1 = MLfunc(correct_rate_list_1, variables)
            if (reversible): 
                correct_rate_2 = MLfunc(correct_rate_list_2, variables)
            
            setparam_f = False
            if (reversible): setparam_b = False
            
            if _float_approx_equal(rate, correct_rate_1):
                # We have found the right parameters. Set them to the correct reaction
                # We've got the ids, but need to find the right parameter in the params list
                for p in params:
                    if not setparam_f and p[0] == params_id[0]:
                        if (p[1] == 'not reac'):
                            raise NotImplementedError("Reaction rate parameter for reaction '%s' has unexpected units."%reac.getId())
                        # 'unknown' basically means that the units were not specified, so we have to use model units
                        if (p[1] == 'unknown'):
                            order = len(r_f.getLhs())
                            # I am not clear about what to do here - what are the default units for say a second order reaction?
                            # Guess I am going to assume volume units are the same as the model volume units, though 
                            # documentation is unclear on this
                            # We want molar units for STEPS, but remember model volume units are in metres:
                            kvalue = p[2]
                            kvalue /= (self.__time['time'][1]*math.pow(self.__substance_units, order-1))
                            kvalue *= (math.pow(self.__volume_units*1000, order-1))
                            # Not sure what to do about the mole - there seems to be no global substance units
                            p[2] = kvalue
                        r_f.setKName(p[0])
                        r_f.setKValue(p[2])
                        r_f.setUnitsConv(p[3])  
                        reactions.append(r_f)
                        setparam_f  = True
                    elif (reversible):
                        if (not setparam_b and p[0] == params_id[1]):
                            if (p[1] == 'not reac'): # or p[1] != len(r_b.getLhs())):
                                raise NotImplementedError("Reaction rate parameter for reversible reaction '%s' has unexpected units."%reac.getId())
                            # 'unknown' basically means that the units were not specified, so we have to use model units
                            if (p[1] == 'unknown'):
                                order = len(r_b.getLhs())
                                kvalue = p[2]
                                kvalue /= (self.__time['time'][1]*math.pow(self.__substance_units, order-1))
                                kvalue *= (math.pow(self.__volume_units*1000, order-1))
                                p[2] = kvalue                                
                            r_b.setKName(p[0])
                            r_b.setKValue(p[2])
                            r_b.setUnitsConv(p[3])
                            reactions.append(r_b)
                            setparam_b = True
            elif (reversible):
                if _float_approx_equal(rate, correct_rate_2):
                    # Yatta, we have found the right params
                    for p in params:
                        if not setparam_f and p[0] == params_id[1]:
                            if (p[1] == 'not reac'): # or p[1] != len(r_f.getLhs())):
                                raise NotImplementedError("Reaction rate parameter for reaction '%s' has unexpected units."%reac.getId())
                            # 'unknown' basically means that the units were not specified, so we have to use model units
                            if (p[1] == 'unknown'):
                                order = len(r_f.getLhs())
                                kvalue = p[2]
                                kvalue /= (self.__time['time'][1]*math.pow(self.__substance_units, order-1))
                                kvalue *= (math.pow(self.__volume_units*1000, order-1))
                                p[2] = kvalue                            
                            r_f.setKName(p[0])
                            r_f.setKValue(p[2])
                            r_f.setUnitsConv(p[3])  
                            reactions.append(r_f)
                            setparam_f  = True
                        elif (not setparam_b and p[0] == params_id[0]):
                            if (p[1] == 'not reac'):# or p[1] != len(r_b.getLhs())):
                                raise NotImplementedError("Reaction rate parameter order for reversible reaction '%s' has unexpected units."%reac.getId())
                            # 'unknown' basically means that the units were not specified, so we have to use model units
                            if (p[1] == 'unknown'):
                                order = len(r_b.getLhs())
                                kvalue = p[2]
                                kvalue /= (self.__time['time'][1]*math.pow(self.__substance_units, order-1))
                                kvalue *= (math.pow(self.__volume_units*1000, order-1))
                                p[2] = kvalue 
                            r_b.setKName(p[0])
                            r_b.setKValue(p[2])
                            r_b.setUnitsConv(p[3])
                            reactions.append(r_b)
                            setparam_b = True                            
            else:
                raise NotImplementedError("Reaction rate maths in reaction %s is not of expected rate."%reac.getId())
            
            # The following would be a coding error- replace with assert once tested??
            if (setparam_f == False):
                raise NotImplementedError("Failed to match reaction rate maths to parameter for reaction %s"%reac.getId())
            elif (reversible):
                if setparam_b == False:
                    raise NotImplementedError("Failed to match reaction rate maths to parameter for reaction %s"%reac.getId())
        
        return reactions
    
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
                        p_factor *= unitsDef_to_STEPS(mult, sca, exp)
                        p_value *= unitsDef_to_STEPS(mult, sca, exp)
                        if (u.isDimensionless() or u.isItem()):
                            # OK, do nothing
                            continue
                        elif (u.isLitre()):
                            # Convert to m
                            p_factor*=math.pow(0.001, exp)
                            p_value*=math.pow(0.001, exp)
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
                        elif (u.isMole()):
                            if (gotMole): notreacunit = True
                            if (exp%1 != 0): notreacunit  = True
                            if (order != -(exp-1)) and (order != -1): notreacunit = True
                            order = -(exp-1)
                            if (order < 0): notreacunit = True
                            gotMole = True
                        elif (u.isSecond()):
                            if (gotSecond): notreacunit = True
                            if (exp != -1): notreacunit = True
                            gotSecond = True
                        else: 
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
                elif (order != 1 and not (gotMole and (gotLitre or gotMetre))):
                        # Also not a reac unit
                        globPar[p_id] = [p_value, p_factor]
                        globPar_order[p_id] = 'not reac'     
                else:
                    # lets keep in metres
                    globPar[p_id] = [p_value, p_factor]
                    globPar_order[p_id] = order
            else: 
                # No units, assume default, whatever they are 
                print "WARNING: No units specified for parameter '%s'. Default units will be assumed (moles, seconds, meters, etc)."%p_id
                globPar[p_id] = [p_value, p_factor]
                globPar_order[p_id] = 'unknown'
        
        return globPar, globPar_order
    
    ################################################################################################
    
    def _setupReactions(self):
        """
        Create the reactions in STEPS
        """
        for r in self.__reactions:
            kreac = smodel.Reac(r.getName(), self.__mdl.getVolsys(r.getVsys()), lhs = r.getLhs(), rhs = r.getRhs())
            kreac.kcst = r.getKValue()
        
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
                raise NotImplementedError("Rule '%s' is Algebraic."%rule.getId())
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
        update_params={}
        update_spec_ref={}

        variables = {}
        variables.update(self.__time)
        variables.update(self.__species)
        variables.update(self.__globalParameters)
        variables.update(self.__comps)

        for initass in listOfInitialAssignments:
            var = initass.getSymbol()
            math = initass.getMath() 
            # Convert the math to a list
            ia_math_frm = MLtoList(math, self.__function_defs)
            
            # Remember the return will be in SBML units
            value = MLfunc(ia_math_frm, variables)
            
            if (var in self.__species):
                value*=self.__species[var][1]
                update_specs[var] = value
            elif (var in self.__comps):
                value*=self.__comps[var][1]
                update_comps[var] = value
            elif (var in self.__globalParameters):
                value*=self.__globalParameters[var][1]
                update_params[var] = value
            else:
                # Something else, must be species reference. 
                update_spec_ref[var] = value
        
        for us in update_specs:
            self.__species[us][0] = update_specs[us]
        for uc in update_comps:
            self.__comps[uc][0] = update_comps[uc]
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
        
        for comp in self.__comps:
            try: sim.setCompVol(comp, self.__comps[comp][0])
            except: 
                raise NotImplementedError("Failed to set compartment volume during"\
                " Initial Assignment (negative volume?).")
        for specie in self.__species:
            if (self.__species_subst_units[specie] == True):
                try: sim.setCompAmount(self.__species_comps[specie], specie, self.__species[specie][0])
                except: 
                    raise NotImplementedError("Cannot set amount of species '%s" \
                    "' in compartment '%s'. Species may not be reactant or "\
                    "product in any comparment reactions (fractional stoichiometry?), or molecules count may"\
                    " be larger than maximum unsigned integer "\
                    "in stochastic simulation."%(specie, self.__species_comps[specie]) ) 
            else:
                try :sim.setCompConc(self.__species_comps[specie], specie, self.__species[specie][0])
                except: 
                    raise NotImplementedError("Cannot set concentration of species '%s" \
                    "' in compartment '%s'. Species may not be reactant or "\
                    "product in any comparment reactions (fractional stoichiometry?), or molecules count may"\
                    " be larger than maximum unsigned integer "\
                    "in stochastic simulation."%(specie, self.__species_comps[specie]) )
            const_bcs = self.__species_const_bc[specie]
            if (const_bcs[0] == True):
                try : sim.setCompClamped(self.__species_comps[specie], specie, True)
                except: pass
        
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
            if (self.__species_subst_units[s]): 
                # Now directly manipulate the self.__species object
                try: self.__species[s][0] = sim.getCompAmount(self.__species_comps[s], s)
                except: pass
            else: 
                try: self.__species[s][0] = sim.getCompConc(self.__species_comps[s], s)
                except: pass
                
        # Adding time. 
        self.__time['time'][0] = sim.getTime()
        self.__time['t'][0] = sim.getTime()
        
        # Collect all the variables in a dictionary on the fly-easier to group them for one dict to MLfunc
        variables = {}
        variables.update(self.__time)
        variables.update(self.__species)
        variables.update(self.__globalParameters)
        variables.update(self.__comps)
                        
        # I guess updating the rates should come first        
        poplist = []
        # Not easy to get the parameters from STEPS, so store a dictionary of params and new value
        params = {}
        for r_rate in self.__rules_rate:
            var = self.__rules_rate[r_rate][0]
            # MLfunc returned value should be / (timeunits) so convert to /s
            # This is a rate so we multiply by sim dt
            value = (MLfunc(self.__rules_rate[r_rate][1], variables)/self.__time['time'][1])*simdt
            # Bit of a hack, but otherwise can't get the assignment and rate rules to match up
            if (sim.getTime() == 0): value = 0.0
            if (var in self.__species):
                # Convert value returned in SBML units to STEPS untis
                value*=self.__species[var][1]
                if (self.__species_subst_units[var] == True):
                    try: 
                        before_amount = sim.getCompAmount(self.__species_comps[var], var)
                        sim.setCompAmount(self.__species_comps[var], var, value+before_amount)
                    except: 
                        raise NotImplementedError("Failed to set amount during rate-rule '%s' (negative amount?)."%r_rate)
                else:
                    try: 
                        before_conc = sim.getCompConc(self.__species_comps[var], var)
                        sim.setCompConc(self.__species_comps[var], var, value+before_conc)
                    except: 
                        raise NotImplementedError("Failed to set concentration during rate-rule '%s' (negative conc?)."%r_rate)
            elif (var in self.__comps):
                # Convert returned value in SBML units to STEPS units
                value*=self.__comps[var][1]
                before_vol = sim.getCompVol(var)
                try: sim.setCompVol(var, before_vol+value)
                except:
                    raise NotImplementedError("Failed to set compartment volume during rate-rule '%s' (negative volume?)."%r_rate)
            elif (var in self.__globalParameters):
                # Convert value from SBML units to STEPS units
                value*=self.__globalParameters[var][1]
                if (value < 0):
                    raise NotImplementedError("Parameter value evaluates to negative value in rate-rule '%s'."%r_rate)
                    
                # Note: Not updating the global parameters here because they may appear in other rate rules
                params[var] = value
                # Now need to convert to molar, STEPS units from metre units
                order = self.__glob_params_order[var]

                setorder = False
                if (type(order) in (int, float)): setorder = True
                # Need to try and get all the reactions with this parameter and change- store in a list first
                # Store the order too for scaling from metres to litres
                k_reacs = []
                for r in self.__reactions:
                    if (r.getKName() == var): 
                        k_reacs.append(r.getName())
                        ord_tmp = len(r.getLhs())
                        if setorder: assert ord_tmp == order 
                        else: 
                            order = ord_tmp
                            setorder = True
                if (k_reacs): 
                    assert setorder
                    # Finally can convert to STEPS units
                    value *= (math.pow(self.__volume_units*1000, order-1))
                    #fac = math.pow(0.001, order-1)
                    #value = value/fac
                    for k_r in k_reacs:
                        for c in self.__comps:
                            before_val = sim.getCompReacK(c, k_r)
                            sim.setCompReacK(c, k_r, value+before_val)
        
        for p in poplist:
            self.__rules_rate.pop(p)
        # Now assignment rules: just like rate rules, but no need to multiply by a dt
        poplist = []
    
        for a_rate in self.__rules_ass:
            var = self.__rules_ass[a_rate][0]
            # MLfunc given stuff in STEPS units, so returned value should be in correct units
            # This is a rate so we multiply by sim dt
            value = MLfunc(self.__rules_ass[a_rate][1], variables)
            if (var in self.__species):
                # Convert value returned in SBML units to STEPS untis
                value*=self.__species[var][1]
                if (self.__species_subst_units[var] == True):
                    try: 
                        sim.setCompAmount(self.__species_comps[var], var, value)
                    except: 
                        raise NotImplementedError("Failed to set amount during assignment-rule '%s' (negative amount?)."%a_rate)
                else:
                    try: 
                        sim.setCompConc(self.__species_comps[var], var, value)
                    except: 
                        raise NotImplementedError("Failed to set concentration during assignment-rule '%s' (negative conc?)."%a_rate)
            elif (var in self.__comps):
                # Convert returned value in SBML units to STEPS units
                value*=self.__comps[var][1]
                try: sim.setCompVol(var, value)
                except:
                    raise NotImplementedError("Failed to set compartment volume during assignment-rule '%s' (negative volume?)."%a_rate)
            elif (var in self.__globalParameters):
                # Convert value from SBML units to STEPS units
                value*=self.__globalParameters[var][1]
                if (value < 0):
                    raise NotImplementedError("Parameter value evaluates to negative value in  assignment-rule '%s'."%a_rate)
                # Note: Not updating the global parameters here because they may appear in other rate rules
                params[var] = value
                # Now need to convert to molar, STEPS units from metre units
                order = self.__glob_params_order[var]

                setorder = False
                if (type(order) in (int, float)): setorder = True
                # Need to try and get all the reactions with this parameter and change- store in a list first
                # Store the order too for scaling from metres to litres
                k_reacs = []
                for r in self.__reactions:
                    if (r.getKName() == var): 
                        k_reacs.append(r.getName())
                        ord_tmp = len(r.getLhs())
                        if setorder: assert ord_tmp == order 
                        else: 
                            order = ord_tmp
                            setorder = True
                if (k_reacs): 
                    assert setorder
                    # Finally can convert to STEPS units
                    #fac = math.pow(0.001, order-1)
                    #value = value/fac
                    value *= (math.pow(self.__volume_units*1000, order-1))
                    for k_r in k_reacs:
                        for c in self.__comps:
                            sim.setCompReacK(c, k_r, value)
        
        for p in poplist:
            self.__rules_rate.pop(p)
        
        # Update the simVars dictionary again, this time with any new compvols too. This second update is necessary
        # because they can't be updated on the fly
        for s in self.__species:
            if (self.__species_subst_units[s]): 
                # Now directly manipulate the self.__species object
                try: self.__species[s][0] = sim.getCompAmount(self.__species_comps[s], s)
                except: pass
            else: 
                try: self.__species[s][0] = sim.getCompConc(self.__species_comps[s], s)  
                except: pass
        for c in self.__comps:
            self.__comps[c][0] = sim.getCompVol(c)
        for p in params:
            self.__globalParameters[p][0] = params[p]
        
        variables.update(self.__time)
        variables.update(self.__species)
        variables.update(self.__globalParameters)
        variables.update(self.__comps)
        
        #######################
        ####### EVENTS ########
        #######################        
        for ev_trig in self.__evnts_trig:
            if (MLfunc(self.__evnts_trig[ev_trig], variables)):
                if (self.__evnts_flip[ev_trig] == False):
                    self.__evnts_flip[ev_trig] = True
                    # Store the time the event 'kicked-in', only if it's not already in the queue
                    if not (self.__evnts_fire.has_key(ev_trig)): 
                        self.__evnts_fire[ev_trig] = self.__time['time'][0]
                        self.__evnts_dl[ev_trig][1] = MLfunc(self.__evnts_dl[ev_trig][0], variables)*self.__time['time'][1]
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
        params = {}
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
                        if (self.__species_subst_units[var]):
                            try: sim.setCompAmount(self.__species_comps[var], var, value)
                            except: 
                                raise NotImplementedError("Cannot set amount of species '%s'" \
                                "in compartment '%s' during event assignment."%(var, self.__species_comps[var]))
                        else:
                            try: sim.setCompConc(self.__species_comps[var], var, value)
                            except: 
                                raise NotImplementedError("Cannot set concentration of species '%s'" \
                                "in compartment '%s' during event assignment."%(var, self.__species_comps[var]))
                    elif (var in self.__comps):
                        # Convert to STEPS units
                        value*=self.__comps[var][1]
                        try: sim.setCompVol(var, value)
                        except:
                            raise NotImplementedError("Failed to set compartment volume during event assignment (negative volume?).")
                    elif (var in self.__globalParameters):
                        # Convert to STEPS units
                        value*=self.__globalParameters[var][1]
                        # Note: Not updating the global parameters here because they may appear in other rate rules
                        params[var] = value  
                        order = self.__glob_params_order[var]
                        
                        setorder = False
                        if (type(order) in (int, float)): setorder = True
                        # Need to try and get all the reactions with this parameter and change- store in a list first
                        k_reacs = []
                        for r in self.__reactions:
                            if (r.getKName() == var): 
                                k_reacs.append(r.getName())
                                ord_tmp = len(r.getLhs())
                                if setorder: assert ord_tmp == order 
                                else: 
                                    order = ord_tmp
                                    setorder = True
                        if (k_reacs): 
                            assert setorder
                            # Finally can convert to STEPS units
                            #fac = math.pow(0.001, order-1)
                            #value = value/fac
                            value *= (math.pow(self.__volume_units*1000, order-1))
                            # NOTE: Don't need to do time, because that will have been taken care of
                            for k_r in k_reacs:
                                for c in self.__comps:
                                    sim.setCompReacK(c, k_r, value)                                             
                    else: 
                        raise NotImplementedError("Currently only assignment rules referencing species" \
                            ", compartments and parameters are supported.")
                
                poplist.append(ev) 
                self.__evnts_vals[ev] = {}
        
        # Get rid of fired events
        for pevent in poplist:
            self.__evnts_fire.pop(pevent)
        
        # Update stuff again ready for the next time 
        for c in self.__comps:
            self.__comps[c][0] = sim.getCompVol(c)
        for p in params:
            self.__globalParameters[p][0] = params[p]
        
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

    def getName(self):
        return self.__name
    
    def setName(self, value):
        self.__name = value
    
    def getLhs(self):
        return self.__lhs
    
    def getRhs(self):
        return self.__rhs
    
    def setLhs(self, value):
        self.__lhs = value
    
    def setRhs(self, value):
        self.__rhs = value
    
    def setUnitsConv(self, value):
        self.__unitsConv = value
    
    def getUnitsConv(self):
        return self.__unitsConv
    
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
####################################################################################################
