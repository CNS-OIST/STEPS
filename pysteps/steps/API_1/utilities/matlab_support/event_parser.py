import os
import ast
import operator as op
import logging

# dictionary of supported operators
# ops = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul, ast.Div: op.truediv, ast.Pow: op.pow, ast.USub: op.neg}


# list of operators we support:
#                         Matlab | Python

# Matrix operators:
# matrix multiplication *        | element wise multiplication *
# element wise multiplication .* | element wise multiplication *
#                              / | /
#                             ./ | / 

# Relational operators:
# ~= | !=

# Logical operators
# ~ | !

# Matlab operators that will be ignored
# left division \
# element wise left division .\
# Transpose .'
# complex conjugate transpose '
allowed_operators = {'+': '+', '-':'-', '.*':'*', '*':'*', './':'/', '/':'/', '.^':'**', '^':'**', '==':'==', '~=':'!=', '>':'>', '>=':'>=', '<':'<', '<=':'<=', '&':'and', '&&':'and', '|':'or', '||':'or', '~':'not'}

# excluded operators
excluded_operators =['\\', '.\\', '.\'', '\'']

# get the logging object. The logging config is
# set in he top level module of STEPS for Matlab
logger = logging.getLogger()

# store the absolute path of this module
# in a module level variable such that Matlab
# can read it
abs_path = os.path.dirname(__file__)


class ExprParser(ast.NodeVisitor):
    """Parse a Matlab mathematical expression defining 
    an event and convert it to the corresponding Python
    STEPS code.
    """

    def __init__(self, species, compartments, rates_to_reactions, reac_to_comp=None, sreac_to_patch=None, kcst_si_factor=None):
        self.species = species
        self.compartments = compartments
        self.rates_to_reactions = rates_to_reactions

        self.reac_to_comp = reac_to_comp
        self.sreac_to_patch = sreac_to_patch
        self.kcst_si_factor = kcst_si_factor

        # flag to indicate if we found a valid
        # assignment to a species or reaction rate
        self.assign_flag = False

        # flag to indicate if we shall generate code
        self.generate_flag = False

        self.generated_code = ''
        logger.debug('Expression Parser created.')
        # done __init__

    def generic_visit(self, node):        
        """ Visit a node of the AST of the event string.
        At least one of the nodes needs to be an assignment.
        Otherwise the state of the STEPS simulation is
        not modified.
        """
        logger.debug('generic_visit ' + type(node).__name__)
       
        if(isinstance(node, ast.Name)):
            logger.debug('call isSpecies of ' + str(node.id))
            # FIXME: replace either the isspecies or whatami method
            ##self.isSpecies(node.id)
            self.whatami(node.id)

        # there need to be at least one assignment whose 
        # target is a species (or parameter)
        if(isinstance(node, ast.Assign)):
            # loop over all the targets
            for t in node.targets:
                if(isinstance(t, ast.Name)):
                    logger.debug('Hit variable ' + str(t.id) + ', ' + str(t.__dict__))
                # if(isinstance(t.value, ast.Name)):
                ##comp_spec = self.isSpecies(t.id)
                comp_spec = self.whatami(t.id)
                self.assign_flag = True

        ast.NodeVisitor.generic_visit(self, node)
        # done generic_visit

    def generatePython(self, expr_str):
        """ Generate Python code from the given Matlab string, 
        not from the AST tree.
    
        It is assumed the expression has been parsed successfully.
        """
        logger.debug('generate_from_string expr: ' + expr_str)
        code = ''

        # separate lines at line ends and semicolons
        lines = expr_str.split('\n')

        for l in lines:
            l = l.strip('\n')
            l = l.strip(';')

            # separate assignments
            # into target and expression
            tokens = l.split('=')
            
            # NOTE: tokens[0] needs to be a species
            # but since the expr has been parsed succesfully,
            # we can assume this

            # get species and compartment - this is the target for the assignment
            
            logger.debug('generate_from_string -  tokens: ' + str(tokens))

            # FIXME: I am missing reaction rates here
            # check for reaction rate
            if(self.isRate(tokens[0])):
                rate_flag = True
            else:
                rate_flag = False

            tmp_flag = False
            if(not rate_flag):
                try:
                    # FIXME: why is this true for a reaction rate?
                    comp_spec = self.isSpecies(tokens[0])
                except RuntimeError:
                    tmp_flag = True
            lhs = tokens[0]
            rhs = tokens[1]

            # extrat the compartment and species from tokens[1]
            # (the expression of the assignment), we need to 
            # generate STEPS code to either retrieve the 
            # molecule count or to convert it into an array ref.
            splits = [tokens[1]]
            for a in allowed_operators.values():
                for s in splits:
                    tmp = s.split(a)
                    if(len(tmp) >= 2):
                        splits.remove(s)
                        splits += tmp
            for i in range(len(splits)):
                splits[i] = splits[i].strip()

            # either c__s or s
            assign_ref = {}
            for s in splits:
                if(s.find('__') > 0):
                    tokens = s.split('__')
                    if(not tokens[0].strip() in self.compartments or not tokens[1].strip() in self.species):
                        logger.error('Unknown compartment or species ' + s)
                        raise RuntimeError('Unknown compartment or species ' + s)
                    else:
                        # build the refeence
                        c_index = self.compartments.index(tokens[0])
                        s_index = self.species.index(tokens[1])
                        assign_ref[s] = 'data[' + str(c_index) + '][' + str(s_index) + ']'
                else:
                    if(not s in self.species):
                        # BUG: not yet finished
                        pass # a value
                        # raise RuntimeError('Unknown species ' + s)
                    else:
                        # build the reference to the data returned by STEPS       
                        assign_ref[s] = 'data[0][' + str(self.species.index(s)) + ']'

            # build Python variables
            for a in assign_ref.keys():
                code += a + ' = ' + assign_ref[a] + '\n'
            
            code += str(lhs) + ' = ' + str(rhs) + '\n'
            
            # not tmp_flag: not an assignment but a species
            if not tmp_flag and not rate_flag:
                code += 'sim.setCompSpecCount("' + comp_spec[0].strip() +'", "' + comp_spec[1].strip() + '", ' + str(lhs) + ')\n' 
            elif rate_flag:
                # check if it is a reaction rate
                if lhs.strip() in self.rates_to_reactions:
                    # we need to find the reaction
                    # we can search the defined Python variable names
                    # extract them from the STEPS model? This seems to be too much effort,
                    # since we should know them      
                    # NOTE: reaction indices start at 1 in the reaction system
                    reaction_index = [i for i,val in enumerate(self.rates_to_reactions) if val==lhs.strip()]
                    # there should only be one
                    if(len(reaction_index) != 1):
                        raise RuntimeError('Reaction index ambigous or empty: ' + str(reaction_index))

                    # end fi
                    # check if it is a reaction or surface reaction
                    reaction_index = reaction_index[0]
                    if('r' + str(reaction_index + 1) in self.reac_to_comp):
                        code += 'r' + str(reaction_index + 1) + '_var = ' + str(self.kcst_si_factor['r' + str(reaction_index + 1)]) + ' * (' + str(rhs) + ')\n'
                        code += 'sim.setCompReacK("' + self.reac_to_comp['r' + str(reaction_index + 1)] + '", "r' + str(reaction_index + 1) + '", r' + str(reaction_index + 1) + '_var)\n'
                    elif('r' + str(reaction_index + 1) in self.sreac_to_patch):
                        code += 'r' + str(reaction_index + 1) + '_var = ' + str(self.kcst_si_factor['r' + str(reaction_index + 1)]) + ' * (' + str(rhs) + '\n'
                        code += 'sim.setPatchSReacK("' + self.sreac_to_patch['r' + str(reaction_index + 1)] + '", "r' + str(reaction_index + 1) + '", r' + str(reaction_index + 1) + '_var)\n'
                    else:
                        raise RuntimeError('Compartment or patch not found.')                    

                else:
                    # if not raise an error
                    logger.error('Unknown reaction rate: ' + lhs)
                    raise RuntimeError('Unknown reaction rate: ' + lhs)
            else:
                # FIXME: this should never happen -> inefficient code
                logger.error('NO flag set.')
                raise RuntimeError('No flag set')
        return(code)
        # done generate_from_string

    def getCode(self):
        """Accessor method returning the generated Python code.
        """
        return(self.generated_code)
        # done getCode

    def whatami(self, name):
        tokens = str.split(name, '__')

        if(len(tokens) == 1):
            # can be species or reaction rate
            if(name.strip() in self.rates_to_reactions):
                return({'rate':name.strip()})
            elif(name.strip() in self.species):
                return({'specie':name.strip()})
            else:
                logger.error('Unknown species or rate: ' + name)
                raise RuntimeError('Unknown species or rate: ' + name)
        elif(len(tokens) == 2):
            if(tokens[0].strip() in self.compartments and tokens[1].strip() in self.species):
                return({'specie':[tokens[0].strip(), tokens[1].strip()]})
            else:
                logger.error('Unknown compartment or specie: ' + str(tokens))
                raise RuntimeError('Unknown compartment or specie: ' + str(tokens))
        # done whatami

    def isRate(self, name):
        logger.debug('israte list of reaction rates: ' + str(self.rates_to_reactions))
        if(name.strip() in self.rates_to_reactions):
            return name

        return False
        # done isRate)

    # what happens if we got a reaction rate?
    # also considers reaction rates
    def isSpecies(self, name):
        tokens = str.split(name, '__')

        if(len(tokens)== 1 and name.strip() in self.species):
            return(['', name])
        elif(len(tokens)== 2 and tokens[0].strip() in self.compartments and tokens[1].strip() in self.species):
            return(tokens)
        else:
            logger.error('Unknown compartment or species: ' + name)
            raise RuntimeError('Unknown compartment or species: ' + name)
        # done isSpecies

    def replaceOperators(self, matlab_expr):
        """ Method substituting Matlab operators with Python ones.
        
        Requirement:
        the operators need to be sorted in a way that the ones
        being a prefix or postfix of another come last.
        """

        # remove ; at the end of a Matlab line
        matlab_expr = matlab_expr.strip('\n')
        matlab_expr = matlab_expr.strip(';')
        
        # check for operators we do not support
        # this is a simple list
        for o in excluded_operators:
            if(matlab_expr.find(o) >= 0):
                logger.error('Unsupported operator: ' + str(o))
                raise TypeError('Unsupported operator: ' + str(o))

        # substitute operators
        for o in allowed_operators.keys():
            matlab_expr = matlab_expr.replace(o, allowed_operators[o])

        logger.debug('replaced operators: ' + matlab_expr)

        return(matlab_expr)
        # done replace operators

    def replaceDot(self, matlab_expr):
        """Replace . notation for compartment.species
        since in the Python AST this denotes an object attribute.

        The substitute string is __ which should never occure
        in a model.
        """

        matlab_expr = matlab_expr.strip('\n')
        matlab_expr = matlab_expr.strip(';')        

        # get all combinations of compartments and species
        for c in self.compartments:
            for s in self.species:
                if(matlab_expr.find(c + '.' + s) >= 0):
                    matlab_expr = matlab_expr.replace(c + '.' + s, c + '__' + s)
        
        return(matlab_expr)
        # done replaceDot

    def toPython(self, matlab_expr, generate=False):
        """ Main entry point. This method parses a
        given Matlab string and returns the generated
        Python code.
        """

        logger.debug('toPython: ' + str(matlab_expr))

        # expression needs to be a string
        if(not type(matlab_expr) == type('')):
            logger.error('Expected a string, got ' + str(type(matlab_expr)) + ': ' + str(matlab_expr))
            raise RuntimeError('Expected a string, got ' + str(type(matlab_expr)) + ': ' + str(matlab_expr))

        self.generate_flag = generate

        matlab_expr = self.replaceDot(matlab_expr)
        matlab_expr = self.replaceOperators(matlab_expr)
        matlab_expr = matlab_expr.strip(' =')
        
        tree = ast.parse(matlab_expr)
        super(ExprParser, self).visit(tree)

        if(not self.assign_flag):
            logger.error('No valid assignment found!')
            raise TypeError('No valid assignent found!')

        # if asked to call the code generation part
        if(generate):
            matlab_expr = self.generatePython(matlab_expr)
            # end if
 
        return(matlab_expr)        
        # done toPython

    # end class visitor

# a few simple test cases
if __name__ == "__main__":
    
    # this needs to be handled differently
    # just for testing
    species = ['A', 'B', 'C', 'r1', 'r2']
    compartments = ['c1', 'c2']
    reaction_rates = ['r1', 'r2']

    # dictionary relating the reaction rates to the reactions
    # keys are the python string names of the reactions
    reaction_to_rates = {'reac1':'r1', 'reac2':'r2'}
    rates_to_reaction = {'r1':'reac1', 'r2':'reac2'}
    
    p = ExprParser(species, compartments, rates_to_reaction)

    try:
        p.toPython("3 + 4") # no assignment -> throw exception   
        print('Incorrect.\n')
    except TypeError:
        print('Correct.\n')

    try:
        p.toPython("C = a + B") # unknown species -> throw exception
        print('Incorrect.\n')
    except RuntimeError:
        print('Correct.\n')

    try:
        p.toPython('C = A + B') # valid
        print('Correct.\n')
    except TypeError as e:
        print('Incorrect ' + str(e) + '.\n')

    try:
        p.toPython('C = A + ( B .^ 2.4)') # valid
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')
    
    try:
        p.toPython('C = A + ( B .^ 2.4);\nA = 1 +2') # valid
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')

    # compartment support
    try:
        p.toPython('c1.C = c1.A + c2.B')
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')

    # test the code generatiin
    try:
        # print(p.generate('C = A + B'))
        print(p.toPython('c1.C = c1.A + c1.B', generate=True))
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')

    try:
        print(p.toPython('c1.B = c1.A;\nc1.C = r2 + c1.B', generate=True))
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')

    # change reaction rates
    # reaction rates are per reaction, not per compartment
    # so we don't need to check for the c.r syntax
    try:
        print(p.toPython('c1.B = c1.A;\nc1.C = r2 + c1.B', generate=True))
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')


    try:
        print(p.toPython('c1.B = c1.A;\nc1.C = r2 + c1.B', generate=True))
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')
    try:
        print(p.toPython('c1.B = c1.A;\nc1.C = r2 + c1.B', generate=True))
        print('Correct.\n')
    except TypeError:
        print('Incorrect.\n')
