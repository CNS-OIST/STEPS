import sys
import os
import re
import string
import datetime
import steps

if(sys.version_info >= (3,0)):
    import configparser as cfg_parser
else:
    import ConfigParser as cfg_parser

class Parser(object):
    """
        Matlab to STEPS parser to convert a Matlab simbiology model into a STEPS model.

        The Matlab simbiology syntax is only partially supported. So far:
        - addspecies
        - addparameter
        - addreaction
        - addkineticlaw
        - addcompartment
        - sbiomodel

        All reactions and species are added to a single STEPS volume system (volsys) and
        a well mixed geometry.

        The generated code is stored internally and not returned as soon as a line of the
        Matlab model is parsed due to differences in how Matlab and STEPS are handling models.
        For details please see the example file.

        generate steps geometry
        add comprtament from the Matlab model
        can a volsys have multiple compartments?
        - compartment belongs to a volsys -> a volys aggregates reaction and diffusion rules within a volume

        set the volume of a compartment?


        Extension mechanism:

            how to make this easily extendable?
            - have a build in set of supported keywords,
            if keyword is unknown, look in a config file for
            user defined keywords

        """
    def __init__(self, config_file='matlab_parser.cfg'):
        """
            Construction arguments:
            * config_file the name of the configuration file

            """

        # load the config file
        cfg = cfg_parser.ConfigParser()
        cfg.read(os.path.dirname(__file__) + '/' + config_file)

        # read the log file configuration
        if(cfg.has_section('logging')):
            self.log = cfg.get('logging', 'destination')
            self.log_name = cfg.get('logging', 'name')
            self.log_timestamp = cfg.getboolean('logging', 'use_timestamp')
        else:
            self.log = None

        # read user defined Matlab keyword
        # section: keywords
        if(cfg.has_section('keywords')):
            keywords = dict(cfg.items('keywords')) 
            print('read additional keywords: ' + str(keywords))
        else:
            keywords = {}

        self.steps_model_module = 'smod'
        self.volsys_name = 'vsys'
        self.model_name = ''

        # dictionary with keyword lists and function references
        self.call = {'addspecies': self.addSpecies, 'addparameter': self.addParameter, 'addreaction': self.addReaction, 'addkineticlaw': self.addKineticLaw, 'set':self.setParameter, 'addcompartment': self.addCompartment, 'sbiomodel': self.addModel}

        # add user defined keywords
        # convert filenames to executable code
        for k in keywords:
            print('keyword: ' + str(k) + ', filename: ' + str(keywords[k]))
            globals()[keyword] = __import__(keywords[k])
            # to call: ret = global()[keyword].function_name(args)   

            # add to the call dictionary
            self.call[keyword] = globals()[keyword].handle            

        # define the STEPS volsys name
        self.volsys = 'volsys'

        # supported units for molecule counts or concentrations

        # supported units for reaction rates

      
        # internal data strutures - should only be accessed
        # via the accessor methods
        self.species = {'null':None}
        self.specie_names = []
        self.parameters = {}
        self.reactions = {}
        self.compartments = []

        # mapping of compartment references to names
        self.compartment_ref_to_name = {}

        self.code = ''
        self.preambleCode = 'import steps.model as ' + self.steps_model_module + '\n'
        self.preambleCode += 'import steps.geom as swm\n'
        self.speciesCode = ''
        self.parameterCode = ''

        # the code to set the initial amounts and to clamp species
        # is a list of strings since we still need to add the STEPS solver reference
        self.setCountCode = []

        # check if log file directory exists
        # if not, create it
        # and open the log file
        if(self.log == 'file'):
            log_directory = os.path.dirname(steps.__file__) + '/../.logs'
            if not os.path.exists(log_directory):
                os.makedirs(log_directory)

            timestamp_str = ''
            if(self.log_timestamp):
                timestamp_str = str(datetime.datetime.now())    
            self.log_file = open(self.log_name + '_' + timestamp_str, 'w')

        # done init

    # ----------------------------
    # clean-up
    # ----------------------------

    def __del__(self):
        """
            Closes the log file.
            """
        self.log_file.close()

        # done __del__

    # ------------------------------------
    # internal helper functions
    # -----------------------------------

    # TODO: move them out off the class
    def argmax(self, lst):
        if(not type(lst) == type([])):
            raise Exception('argmax requires list as input.')

        return lst.index(max(lst))

    def argmin(self, lst):
        if(not type(lst) == type([])):
            raise Exception('argmin requires list as input.')

        return lst.index(min(lst))

    def logger(self, message):
        if(not self.log):
            pass
        elif(self.log == 'file'):
            self.log_file.write(str(message) + '\n')        
        elif(self.log == 'screen'):
            print(str(message))
        else:
            raise Exception('Unknown log level: ' + str(self.log))

        # done logger

    # -------------------------------------
    # interface functions
    # -------------------------------------

    def parse(self, lines):
        """
            One of the two main entry points.

            Parse a list of lines such as an entire Matlab
            simbiology model.

            The generated STEPS code is not immediately returned but stored internally.

            Arguments:
            * lines list of lines
            """
        for line in lines:
            self.parseLine(line)
        # done parse

    # a line seems to be
    #    LINE: Python reference EQUALS EXPR
    #    EXPR: KEYWORD ( args );
    #    KEYWORD: addspecies addparameter addreaction

    # this may return an empty line
    def parseLine(self, line):
        """
            Parse a single line of a Matlab simbiology model.

            The generated STEPS code is not immediately returned but stored internally.

            Arguments:
            * line a single line of the Matlab simbiology model.
            """
        # forward comments - we would like to preserve them 
        if(line[0] == '#'):
            return(line)
        elif(line[0] == '%'):
            return('#' + line[1:])
        elif(line[0] == '\n'):
            return(line)
        elif(line == ''):
            return(line)

        # remove whitespaces in front of the delimiters = . (
        line = re.sub(r'\s([=.(](?:\s|$))', r'\1', line)
        
        # strip trailing comments
        tokens = line.split('%')
        line = tokens[0]
        
        # find the delimiters
        indices = []
        indices.append(line.find('='))
        indices.append(line.find('('))
        indices.append(line.find(' ')) # does not work
        indices.append(line.find('.'))

        # replace all -1 entries with max int
        for i in range(4):
            if(indices[i] == -1):
                indices[i] = sys.maxint
        delimiter = self.argmin(indices)
        
        # filter for a single keyword without args
        if(indices[delimiter] == sys.maxint):
            delimiter = 4

        reference = None
        keyword = None
        rhstokens = None

        # case 1: reference = keyword
        if(delimiter == 0): # =
            tokens = line.split('=')

            # separate the reference to be returned
            reference = tokens[0].strip()
            
            # split at the first ( - separate keyword and args
            tmp = string.split(tokens[1], '(', 1) 
            keyword = tmp[0].strip()
            
            if(len(tmp) > 1):
                # split the args further
                tmp = re.split(',|;|\n', tmp[1])
            
                rhstokens = filter(bool, tmp)

                # strip \' and \"
                for i in range(len(rhstokens)):
                    rhstokens[i] = rhstokens[i].strip(' \)\'\"')
                rhstokens = filter(bool, rhstokens)

        # case 2: keyword(args)
        elif(delimiter == 1): # (
 
            # split at the first (
            tmp = string.split(line, '(', 1)
            keyword = tmp[0].strip()
            
            if(len(tmp) > 1):
                tmp = re.split(',|;|\n', tmp[1])
            
                rhstokens = filter(bool, tmp)

                # strip \' and \"
                for i in range(len(rhstokens)):
                    rhstokens[i] = rhstokens[i].strip(' \)\'\"')
                rhstokens = filter(bool, rhstokens)

        # case 3: keyword args
        elif(delimiter == 2): # white space
            # args are seperated by a whitespace from the keyword
            # there can only be one arg
            tmp = string.split(line, ' ', 1)

            keyword = tmp[0].strip()
            
            if(len(tmp) > 1):
                tmp = re.split(',|;|\n', tmp[1])
            
                rhstokens = filter(bool, tmp)

                # strip \' and \"
                for i in range(len(rhstokens)):
                    rhstokens[i] = rhstokens[i].strip(' \)\'\"')
                rhstokens = filter(bool, rhstokens)

        # case 4: reference.attribute 
        elif(delimiter == 3):
            self.logger('Access to object attributes not supported. Ignoring.')
            return('# ' + line)
        elif(delimiter == 4): # keyword only
            keyword = line.strip()
        else:
            raise Exception('parseLine unknown delimiter: ' + line)
        
        # call the function to handle the keyword
        # call keyword specific function
        # token[0] is the return value, we need this to identify the reaction, this can be None
        # rhstokens are the arguments of the function call
        try:
            self.logger('dispatch call: ' + str(keyword) + ', reference: ' + str(reference) + ', rhstokens: ' + str(rhstokens))
            self.call[keyword](reference, rhstokens)
        except KeyError:
            self.logger('Matlab keyword ' + str(keyword) + ' not yet supported. Ignoring.')

        # done parseLine
    # ---------------------------------------------------
    # accessor methods
    # ---------------------------------------------------

    def numberOfReactions(self):
        """
            Returns the number of reactions in the model.

            Returns:
            * number of reactions
            """
        return(len(self.reactions))
        # done numberOfReactions

    def numberOfSpecies(self):
        """
            Returns the number of species in the model.

            Returns:
            * number of species
            """
        return(len(self.species))
        # done numberOfSpecies

    def getListOfSpecies():
        """
            Returns the list of species names.
            This is required to query STEPS for the 
            current molecule counts
            without examining the model.

            Returns:
            * list of specie names
            """

        return(self.specie_names)
        # done getListOfSpecies

    def getListOfCompartments():
        return(self.compartments)
        # done getListOfCompartments

    def getPreamble(self):
        return(self.preambleCode + self.parameterCode)
        # done getPreamble

    def getModelCode(self):
        return(self.code)
        # done getModelCode

    def getSetCountCode(self, solver_ref_str):
        if(solver_ref_str == ''):
            raise Exception('getSetCountCode invalid solver string: ' + solver_ref_str)

        code = ''
        for i in self.setCountCode:
            code += str(solver_ref_str) + i

        return(code)
        # done getSetCountCode

    def getSpecieDefinition(self):
       return(self.speciesCode)
       # done getSpecieDefinition)

    # ----------------------------------------------------
    # internal functions, should not be called by the user
    # ----------------------------------------------------


    # adopt to MATLAB behaviour:
    # if species is not known add it automatically
    def addSpecies(self, reference, tokens):
        # Matlab allows to add a species to either
        # a model or a compartment
        # we need to check which one it is.
        
        compartment = tokens.pop(0)

        # get the compartment's name
        compartment = self.compartment_ref_to_name[compartment]

        clamp_all = False
        if(compartment == self.model_name):
            # if the species is clamped, we need to clamp it in all compartments
            # how doe I get if it is clamped at this stage?
            # just set a flag?
            clamp_all = True
        elif(compartment in self.compartments):
            # only clamp it in the given compartment
            clamp_all = False
        else:
            raise Exception('Unknown model or compartment: ' + compartment + ', known: ' + str(self.compartments))

        name = tokens.pop(0)
        
        # from here on it is more flexible
        # it can be the initial amount or a key / value pair
        factor = 1
        count = 0
        clamping_code = []
        molarity = 'notset' # I need a ternary laogic, simply true or false does not work
        while(tokens):  
            token = tokens.pop(0)
            if(token == 'InitialAmount'):
                count = float(tokens.pop(0))
            elif(token == 'InitialAmountUnits'):
                # TODO work out scaling factor
                # we are acepting: micromolarity or molecules
                factorstr = tokens.pop(0)
                if(factorstr == 'micromolarity'):
                    factor = 1e-6
                    molarity = 'True'
                elif(factorstr == 'molecules'):
                    factor = 1
                    molarity = 'False'
                elif(factorstr == 'molecule'):
                    factor = 1
                    molarity = 'False'
                else:
                    raise Exception('Invalid unit: ' + str(factorstr))
            elif(token == 'ConstantAmount'):
                # set clamped flag in steps
                # FIXME: I need the compartment
                
                # set the flag
                flag = tokens.pop(0)
                if(flag == 'true' or flag == '1' or flag == 'TRUE' or flag == 'True'):
                    flag = 'True'
                elif(flag == 'false' or flag == '0' or flag == 'FALSE' or flag == 'False'):
                    flag = 'False'
                else:
                    raise Exception('Invalid flag: ' + str(flag))
                if(clamp_all):
                    for c in self.compartments:
                        clamping_code.append('.setCompSpecClamped(\'' + c + '\', \'' + name + '\' , ' + str(flag) + ')\n')
                else:
                    clamping_code.append('.setCompSpecClamped(\'' + compartment + '\', \'' + name + '\' , ' + str(flag) + ')\n')
            elif(token == 'Name'):
                # same as an unpaired token[1]
                name = tokens.pop(0)
            elif(token == 'BoundaryCondition' or token == 'Notes' or token == 'Parent' or token == 'Tag' or token == 'Type'):
                tokens.pop(0) # ignore
            else:
                # this needs to be the initial amount
                count = tokens.pop(0)

        code = ''
        code += name + ' = ' + self.steps_model_module + '.Spec(\'' + name + '\',' + str(self.model_name) + ')\n'
        if(molarity == 'True'):
            self.setCountCode.append('.setCompSpecAmount(\'' + compartment + '\' ,\'' + name + '\', ' + str(count * factor) + ')\n')
        elif(molarity == 'False'):
            self.setCountCode.append('.setCompSpecCount(\'' + compartment + '\' ,\'' + name + '\', ' + str(count) + ')\n')
        # self.setCountCode += 'print(\'' + name + ' added.\')\n'
        self.setCountCode += clamping_code    

        self.speciesCode += code

        # internal book keeping
        # FIXME: is this ever used?
        self.species[reference] = count
            
        self.specie_names.append(name)

        # done addSpecies

    # add parameter into a dictionary,
    # we need them as reaction rates when adding the reactions
    def addParameter(self, reference, tokens):
        name = tokens.pop(0)
        parameter_name = tokens.pop(0)

        value = 0
        factor = 1
        while(tokens):
            token = tokens.pop(0)
            if(token == 'Value'):
                value = float(tokens.pop(0))
            elif(token == 'ConstantValue'):
                # ignore it for now, just remove the boolean from the stack
                tokens.pop(0)
                pass
            elif(token == 'ValueUnits'):
                parameterUnit = tokens.pop(0)
                print('addParameter unit: ' + str(parameterUnit))
                if(parameterUnit != '1' and parameterUnit != '1/second' and parameterUnit != '1/(micromolarity*second' and parameterUnit != '1/(nanomolarity*second' and parameterUnit != 'micromole/second' and parameterUnit != 'nanomol/second'):
                    raise Exception('addParameter: invalid parameter unit: ' + str(parameterUnit))

                factor = 1
                if(parameterUnit == '1/(micromolarity*second' or parameterUnit == 'micromole/second'):
                    factor = 1e-6
                    molarity = 'True'
                elif(parameterUnit == '1/(nanomolarity*second' or parameterUnit == 'nanomol/second'):
                    factor = 1e-9
                    molarity = 'True'
                # no else path since already covered above

            elif(token == 'Name'):
                name = tokens.pop(0)
            elif(token == 'Notes' or token == 'Parent' or token == 'Tag' or token == 'Type' or token == 'UserData'):
                tokens.pop(0) # ignore
            else:
                # this must be the value)
                value = token

        # scale the value
        value *= factor

        # parameter definition
        self.parameterCode += parameter_name + ' = ' + str(value) + '\n'

        # add to the dictionary of parameters
        self.parameters[parameter_name] = value
        # done addParameter


    # throw exception if anything else but the law of mass action
    # I need to check if there is already a matching parameter around, otherwise throw exception
    
    # Properties:
    #  Active
    #  KineticLaw
    #  Name
    #  Products
    #  Reactants
    #  ReactionRate
    #  Reversible (-> I need to add a second reaction in steps)
    #  Stoichiometry


    # adopt to MATLAB behaviour:
    # if any of the species does not yet exist,
    # add it automatically
    def addReaction(self, reference, tokens):    
        # tokens.pop(0) # remove the keyword

        model_name = tokens.pop(0) # ignore

        definition= tokens.pop(0)

        # ReactionValue
        # see if we got a reversible reaction
        if(definition.find('<->') >= 0):
            raise Exception('Reversible reactions not yet supported')

        if(definition.find('->') == -1):
            raise Exception('Invalid reaction definition: ' + str(token))

        # build the lhs and rhs
        # TODO: can there be a <-> (reversible reaction?)
    
        constituents = definition.split('->')
        reactants = constituents[0].split('+')
        products = constituents[1].split('+')

        # check if we got all the species
        # and replace names with references
        code = 'lhs = ['
        iterator = iter(reactants)
        code1 = ''
        while(True):
            try:
                specie = iterator.next().strip()
                
                # handle the stoichiometric factor
                # stoichiometric factors need to be positive integers
                # according to ethe Matlat simbio manual
                # "spaces are required before and after species names and stoichiometric values"
                tmp = specie.split()
                stoch_factor = 1
                specie_name = ''
                if(len(tmp) == 1):
                    specie_name = specie.strip()
                elif(len(tmp) == 2):
                    if(tmp[0].isdigit()):
                        stoch_factor = int(tmp[0])
                        specie_name = tmp[1]
                    elif(tmp[1].isdigit()):
                        stoch_factor = int(tmp[1])
                        specie_name = tmp[0]
                    else:
                        raise Exception('addReaction, stoichiometric factor needs to be a positive integer: ' + str(specie))
                else:
                    raise Exception('addReaction, invalid species name: ' + str(specie))

                if(not self.species.has_key(specie_name)):
                    # STEPS species require the model, not the compartment as Matlab
                    comp_str = ''
                    if(len(self.compartments) == 1):
                        comp_str = self.compartments[0]
                    else:
                        comp_str = self.model_name

                    # FIXME: we use the same string as name and object reference.
                    # This seem to work for now but might need checking in the future.
                    self.addSpecies(str(specie_name), [comp_str,  str(specie_name)])
                for i in range(stoch_factor):
                    if(code1):
                        code1 += ','
                    if(str(specie_name) != 'null'):
                        code1 += str(specie_name) + ' '
            except StopIteration:
                break
        code += code1 + '], rhs = ['    

        iterator = iter(products)
        code1 = ''
        while(True):
            try:            
                specie = iterator.next().strip()
                
                # handle the stoichiometric factor
                # stoichiometric factors need to be positive integers
                # according to ethe Matlat simbio manual
                # "spaces are required before and after species names and stoichiometric values"
                tmp = specie.split()
                stoch_factor = 1
                specie_name = ''
                if(len(tmp) == 1):
                    specie_name = specie.strip()
                elif(len(tmp) == 2):
                    if(tmp[0].isdigit()):
                        stoch_factor = int(tmp[0])
                        specie_name = tmp[1]
                    elif(tmp[1].isdigit()):
                        stoch_factor = int(tmp[1])
                        specie_name = tmp[0]
                    else:
                        raise Exception('addReaction, stoichiometric factor needs to be a positive integer: ' + str(specie))
                else:
                    raise Exception('addReaction, invalid species name: ' + str(specie))

                if(not self.species.has_key(specie_name)):
                    # STEPS species require the model, not the compartment as Matlab
                    self.logger('addReaction adding missing species: ' + str(specie_name))
                    comp_str = ''
                    if(len(self.compartments) == 1):
                        comp_str = self.compartments[0]
                    else:
                        comp_str = self.model_name

                    self.addSpecies(str(specie_name), [comp_str,  str(specie_name)])
                for i in range(stoch_factor):
                    if(code1):
                        code1 += ','
                    if(str(specie_name) != 'null'):
                        code1 += str(specie_name) + ' '
            except StopIteration:
                break
        code += code1 + ']'
        
        if(not self.reactions.has_key(str(reference).strip())):
            self.reactions[str(reference).strip()] = {}

        self.reactions[str(reference).strip()]['definition'] = code

        # iterate over remaining tokens
        while(tokens):
            token = tokens.pop(0)

            if(token == 'Active'):
                # todo flip reaction active / onsactive
                pass
            elif(token == 'Name'):
                token.pop(0)
            elif(token == 'Reaction'):
                # copy constructor?
                self.logger('addReaction not yet supported copy constructor called: ' + str(tokens.pop(0)))

            elif(token == 'ReactionRate'):
                # FIXME: this can be an equation
                # FIXME: add scaling here?
                self.reactions[reaction_name]['rate'] = tokens.pop(0)
                self.reactions[reaction_name]['state'] = 12
            elif(token == 'Reversible'):
                reversible = tokens.pop(0)
                if(reversible == 'True' or reversible == '1' or reversiblte == True or reversible == 1):
                    raise Exception('Steps does not support reversible reactions.')

            elif(token == 'KineticLaw' or token == 'Notes' or token == 'Parent' or token == 'Products' or token == 'Reactants' or token == 'Stoichiometry' or token == 'Tag' or token == 'Type' or token == 'UserData'):
                tokens.pop(0) # ignore

        ## self.code += code
        # check if we got everything to build the reaction
        #self.reactions[str(reference).strip()]['state'] += 3
        #if(self.reactions[str(reference).strip()]['state'] == 15):
        #    return(self.buildReaction(str(reference).strip()))
        # return('')
        # done addReaction

    def addKineticLaw(self, reference, tokens):
        # tokens.pop(0) # remopve the key word       
 
        name = tokens.pop(0).strip()

        # get the reaction reference
        if(not self.reactions.has_key(name)):
            raise Exception('Reaction ' + str(name) + ' not yet defined.')

        kinetic_law = tokens.pop(0).strip()
        if(not  kinetic_law == 'MassAction'):
            raise Exception('Only the law of mass action is supported.')

        # ignore all the other tokens for now since we only
        # support mass action)

        self.reactions[name]['kineticlaw_name'] = str(reference).strip()
        # done addKineticLaw

    # this is actually "set" in Matlab
    # we need to look for the attribute 'ParameterVariableNames'
    # which links a named parameter to a reaction as rate constant
    def setParameter(self, reference, tokens):
        # tokens.pop(0) # remove key word
        if(reference):
            raise Exception('No object reference expected: ' + str(reference))

        object_ref = tokens.pop(0)
        tokens = filter(bool, tokens)
         
        while(tokens):
            token = tokens.pop(0).strip(' \'\"')
            if(token == 'ParameterVariableNames'): 
                parameter_name = tokens.pop(0).strip(' \'\"\(\)\{\}\{\}')
            else:
                raise Exception('set: ParameterVariableNames expected.')

        # linear search of all reactions
        # print('reactions: ' + str(self.reactions))
        reaction_name = ''
        for r in self.reactions:
            if(self.reactions[r]['kineticlaw_name'] == str(object_ref)):
                self.reactions[r]['parameter_name'] = parameter_name
                reaction_name = r
        return(self.buildReaction(reaction_name))
        # done setParameter

    # actually assemble the reaction and
    # return the steps code
    def buildReaction(self, reaction_name):
        self.logger('buildReaction: ' + str(reaction_name) + ', ' + str(self.reactions[reaction_name]) +', parameters: ' + str(self.parameters) )
        # self.code += str(reaction_name) + ' = ' + str(self.steps_model_module) + '.Reac(\'' + str(reaction_name) + '\', ' + str(self.volsys_name) + ', ' + str(self.reactions[reaction_name]['definition']) + ', kcst = ' + str(self.parameters[self.reactions[reaction_name]['parameter_name']]) + ')\n' 
        self.code += str(reaction_name) + ' = ' + str(self.steps_model_module) + '.Reac(\'' + str(reaction_name) + '\', ' + str(self.volsys_name) + ', ' + str(self.reactions[reaction_name]['definition']) + ', kcst = ' + str(self.parameters[self.reactions[reaction_name]['parameter_name']]) + ')\n' 
        # return(code)
        # done buildReaction

    # add model to get the Matlab reference to it
    # sbiomodel(name)
    # add the STEPS geometry with it
    # TODO: how and when do I add the tetmesh? 
    def addModel(self, reference, tokens):
        # tokens.pop(0) # remove key word

        self.model_name = reference
        #if(reference):
        #    raise Exception('No object reference expected: ' + str(reference))

        if(len(tokens) != 1):
            raise Exception('Only model name expected: ' + str(tokens))

        code = str(self.model_name) + ' = smod.Model()\n'
        code += 'vsys = smod.Volsys(\'vsys\', ' + str(self.model_name) + ')\n'
        code += 'geometry = swm.Geom()\n'
        code += 'model = ' + str(self.model_name) + '\n'
        self.preambleCode += code
        # no further properties currently supported


    # TODO: support multiple compartments
    # comp = swm.Comp(name, geometry)
    # comp.addVolsys('volsys')
    # comp.setVol(vol)
    
    # matlab syntax:
    # modelObj, name
    # modelObj, name, capacity
    # + name, value pairs

    # supported properties:
    # Capaity
    # CapacityUnit (need to convert it to m^3)
    
    # unsupported properties:#
    # ConstantCapacity will throw an exception since the concept is unknown to STEPS

    # all other will be ignored
    # compartments within compartments are not supported

    def addCompartment(self, reference, tokens):

        if(not reference):
            raise Exception('Object reference expected.')
    
        # to add a compartment a model is required
        if(not self.model_name):
            raise Exception('Model reference required.')

        # the first token are the model and the compartment name
        model_ref = tokens.pop(0)
        
        if(model_ref != self.model_name):
            raise Exception('Compartment does not belong to this model.')

        compartment_name = tokens.pop(0)

        tokens = filter(bool, tokens)

        # name value pairs,
        # if a none recognized token is found it is assumed to be the capacity -> exception if not a number
        capacity = 1.0 # default 0 does not seem to make sense here
        factor = 1.0
        while(tokens):
            token = tokens.pop(0).strip(' \'\"')

            if(token == 'Capacity'):
                capacity = float(tokens.pop(0))

            elif(token == 'CapacityUnits'):
                # supported units:
                # l
                # ml
                # nl
                # liter
                # milliliter
                # nanoliter
                unit = tokens.pop(0).strip(' \'\"')
                print('Capacity unit: ' + str(unit))
                if(unit == 'l' or unit == 'liter'):
                    factor = 1e-3
                elif(unit == 'ml' or unit == 'milliliter'):
                    factor = 1e-6
                elif(unit == 'mul' or unit == 'microliter'):
                    factor = 1e-9
                elif(unit == 'nl' or unit == 'nanoliter'):
                    factor = 1e-12
                else:
                    raise Exception('Unsupported capacity unit: ' + str(unit))
            elif(token == 'ConstantCapacity'):
                raise Exception('Compartments with a constant capacity are not supported in STEPS.')

            else:
                # this must be the capacity
                try:
                    capacity = float(token)
                except:
                    raise Exception('Compartment capacity needs to be a number: ' + str(token))
            
        code = str(reference) + ' = swm.Comp(\'' + str(compartment_name) + '\', geometry)\n'
        code += str(reference) + '.addVolsys(\'vsys\')\n'

        if(capacity):
            self.logger('addCompartment scaling factor: ' + str(factor) + ' scaled volume: ' + str(capacity * factor))
            code += str(reference) + '.setVol(' + str(capacity * factor) + ')\n'
        
        # add this compartment to the list of compartments
        self.compartments.append(compartment_name)
        self.compartment_ref_to_name[reference] = compartment_name
        self.code += code
        # done addCompartment

if __name__ == "__main__":
    # give brief usage string
    pass
 

