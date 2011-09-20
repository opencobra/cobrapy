from time import time
from numpy import zeros, matrix, array
from scipy.sparse import lil_matrix, dok_matrix
from cPickle import load, dump
from collections import defaultdict
def remove_reactions(cobra_model, db_cursor):
    """This function is only for ME matrix calculations and
    internal use by the sybep group.
    
    """
    pass

def create_element_matrix(cobra_model, db_cursor=None, me_matrix=False,
                          cas='sympy'):
    """Constructs a matrix of elements x metabolites for the metabolites
    in cobra_model.  If elemental compositions are not available for all
    metabolites then a symbolic matrix is returned with symbols representing
    the unknown compositions.

    cobra_model:  A cobra.Model object.

    db_cursor:  Internal use only.

    me_matrix: False. Internal use only.

    cas: 'sympy' or 'ginac'.  Specifies the computer algebra system to use.
    'ginac' is the more powerful solver however it is accessed through
    swiginac which isn't the easiest thing to install.

    """
    if cas.lower() == 'sympy':
        from sympy import solve, Matrix, Symbol, Add
        from sympy.core.numbers import Zero, Real
    elif cas.lower() == 'ginac':
        #Symbolic.Matrix might be slow to index so we might use
        #the swiginac matrix or numpy.matrix instead
        from Symbolic import Symbol, Matrix, Expr
                
        
    elements = ('c', 'h', 'o', 'n', 'p', 's', 'z') #z is for generics
    element_dict = dict(zip(elements, range(len(elements))))
    element_matrix = dok_matrix((len(elements), len(cobra_model.metabolites)))
    if db_cursor and me_matrix:
        #Used for current incarnation of ME matrix.  
        known_compositions = set()
        #1. Start off by getting known molecular compositions
        db_cursor.execute('Select Id, (c, h, o, n, p, s) from metabolite')
        metabolite_compositions = dict(db_cursor.fetchall())
        for the_metabolite, the_composition in metabolite_compositions.items():
            if the_metabolite in cobra_model.metabolites:
                known_compositions.add(the_metabolite)
                the_column = cobra_model.metabolites.index(the_metabolite)
                the_composition = eval(the_composition)
                element_matrix.update(dict([((i, the_column), the_composition[i])
                                            for i in range(len(the_composition))]))

        #2. Identify the reactions that produce generics and set the stoichiometry
        #to zero to deal with mass imbalances.
        #
        #This isn't a problem any more as the dummy reactions are not added
        #until after balancing.
        generic_production_reactions = dict([(x, cobra_model.reactions.index(x))
                                             for x in cobra_model.reactions\
                                             if x.startswith('generic_rename')])
        #collect all of the generic_X metabolite Ids and then set masses to 1 Z
        db_cursor.execute('SELECT FU_ID_Generic from Generic_FU')
        generic_metabolites = [x[0] for x in db_cursor.fetchall()]

        known_compositions.update(generic_metabolites)
        element_matrix.update(dict([((elements.index('z'),
                                      cobra_model.metabolites.index(x)), 1.)
                                    for x in generic_metabolites]))
        #Remove the generic production reactions
        for column_index in generic_production_reactions.values():
            the_column = cobra_model._S[:, column_index]
            for row_index in the_column.nonzero()[0]:
                the_column[row_index, 0] = 0

        #3. Remove demand reactions.
        #This isn't a problem any more as the demand reactions are not added
        #until after balancing.
        demand_indices = [cobra_model.reactions.index(x)
                          for x in cobra_model.reactions if 'demand' in x]
        for column_index in demand_indices:
            the_column = cobra_model._S[:, column_index]
            for row_index in the_column.nonzero()[0]:
                the_column[row_index, 0] = 0
        #4. Calculate molecular formula for transcripts and peptides.  This isn't
        #necessary, but will probably make solving the problem easier.

        #known_compositions.update(transcripts and peptides)

        #5. For all metabolites not in known compositions.  It will be
        #necessary to add a symbolic value to the matrix for each element
        #excluding z, i.e.  the_metabolite_(c, h, o, n, p, s, z). 
        metabolite_dict = dict(zip(cobra_model.metabolites,
                                   range(len(cobra_model.metabolites))))
        [metabolite_dict.pop(k)
         for k in known_compositions]
        #If there are any metabolites without compositions build a symbolic matrix
        if len(metabolite_dict) > 0:
            #Build the symbolic elemental composition matrix
            if cas.lower() == 'sympy':
                element_matrix = Matrix(element_matrix.todense())
            elif cas.lower() == 'ginac':
                element_matrix = Matrix(element_matrix.todense().tolist())
            #Now populate with symbols for the_metabolites not in known_compositions
            for the_metabolite, metabolite_index in metabolite_dict.items():
                for the_element, element_index in element_dict.items():
                    element_matrix[element_index,
                                   metabolite_index] = Symbol('%s_%s'%(the_metabolite,
                                                                       the_element))
    else:
        print 'Not yet implemented for anything other than ME'

    return({'elements': element_dict,
            'matrix': element_matrix})





def create_balancing_problem(cobra_model, element_vector,
                             system_type='equations',
                             cas='sympy', print_unbalanced=False):
    """Create a symbolic linear algebra problem to solve for unknown metabolite
    compositions in cobra_model.

      element_vector x cobra_model._S.T = 0

      Note: The problems are typically too big to balance for all elements
      at once, thus the vector should only deal with one of the chemical elements.

    cobra_model: A cobra.Model object

    element_vector: A sympy.Matrix with the element counts corresponding to each
    metabolite in cobra_model.reactions.

    return_type: 'equations' or 'matrix'.  Constructing equations is faster, but
    for specific operations a matrix may be desired.  If 'matrix' then the
    unbalanced reactions are not currently returned.

    cas: 'sympy' or 'ginac'.  Specifies the computer algebra system to use.
    'ginac' is the more powerful solver however it is accessed through
    swiginac which isn't the easiest thing to install.

    print_unbalanced: Boolean.  Indicates whether to print the unbalanced reactions.
    
    """
    if cas.lower() == 'sympy':
        from sympy import solve, Matrix, Symbol, Add
        from sympy.core.numbers import Zero, Real
    elif cas.lower() == 'ginac':
        from Symbolic import Symbol, Matrix, Expr
        from swiginac import lsolve, matrix, symbol
        from swiginac import add as Add
    #  If multiple solutions are available, then we'll have to make
    # sure that the sum of each column in element_matrix is >= 1
    #
    #Now deal with multiplying the two matrices.
    unbalanced_dict = {}
    variable_set = set()
    if system_type == 'equations':
        s_matrix_transpose = cobra_model._S.T #Row access is faster for sparse arrays
        #Multiplying symbolic_element_matrix by s_matrix is not an option, due to
        #memory and speed issues.
        the_system = []


        for i in range(s_matrix_transpose.shape[0]):
            the_column = s_matrix_transpose[i, :] 
            #This is faster than multiplying symbolic_element_matrix by
            #the_column by about 3-fold
            the_indices = the_column.nonzero()[1]
            the_factors = [float(the_column[0, x])
                           for x in the_indices]
            if cas.lower() == 'sympy':
                the_variables = [element_vector[0,j]
                                 for j in the_indices]
            elif cas.lower() == 'ginac':
                the_variables = [element_vector[j]
                                 for j in the_indices]

            the_equation = reduce(lambda x,y: x + y,
                                   map(lambda x, y: x*y, the_factors,
                                       the_variables))
            if cas.lower() == 'sympy':
                #this can probably be streamlined for the different CASes
                if isinstance(the_equation, Add):
                    the_system.append(the_equation)
                    the_atoms = list(the_equation.atoms())
                    [variable_set.add(x) for x in the_atoms
                     if isinstance(x, Symbol)]
                elif not isinstance(the_equation, Zero) and \
                         the_equation != 0:
                    unbalanced_dict.update({cobra_model.reactions[i]:
                                            the_equation})
                    if print_unbalanced:
                        print 'Unbalanced Reaction ' +\
                              '%s: Element %s is %s'%(cobra_model.reactions[i],
                                                      the_element,
                                                      repr(the_equation))
            elif cas.lower() == 'ginac':
                if isinstance(the_equation.data, Add):
                    the_system.append(the_equation)
                    [variable_set.add(x) for x in the_variables
                     if isinstance(x, Expr)]
                elif the_equation.eval() != 0:
                    unbalanced_dict.update({cobra_model.reactions[i]:
                                            the_equation})
                    if print_unbalanced:
                        print 'Unbalanced Reaction ' +\
                              '%s: Element %s is %s'%(cobra_model.reactions[i],
                                                      the_element,
                                                      repr(the_equation))

    else:
        print 'Warning this may take 10 Gb RAM and an hour'
        the_system = element_vector * cobra_model._S.todense()

    return({'variables': variable_set,
            'equations': the_system,
            'unbalanced': unbalanced_dict})


def solve_balance_problem(the_equations, the_variables):
    """Solves a systems of linear equations for the variables. Using
    sympy.

    the_equations: A list of sympy.Add equations.

    the_variables: A list of sympy variables (Symbols, Zero, Real, One, ...)
    
    """
    if cas.lower() == 'sympy':
        from sympy import solve
    elif cas.lower() == 'ginac':
        from swiginac import lsolve as solve
    the_solution = solve(the_equations, the_variables)    
    return the_solution


if __name__ == '__main__':
    from sys import argv
    from os.path import lexists
    from time import time
    ## if not len(argv) == 3:
    ##     print 'Need to call the script with the model file name and element'
    ##     print 'e.g. python balance.py cobra_model.pickle c'
    ## model_file = argv[1]
    ## the_element = argv[2]
    cas='ginac'
    me_matrix = True
    the_element = 'c'
    model_file = '/Users/danie/e/builds/cobra_model.pickle'
    system_type = 'equations'
    element_file = '%s.elements.%s'%(model_file, cas)
    problem_file = '%s.%s_problem.%s'%(model_file,
                                       the_element,
                                       cas)
    solution_file = '%s.%s_solution.%s'%(model_file,
                                       the_element,
                                       cas)


    with open(model_file) as in_file:
        cobra_model = load(in_file)

    if lexists(element_file) and not lexists(problem_file):
        #Only load the element file if a problem file does not
        #already exist to reduce memory usage.
        with open(element_file) as in_file:
            the_elements = load(in_file)
    elif not lexists(problem_file):
        #Only build the elements if the problem file doesn't exist
        import pgdb as PgSQL
        start_time = time()
        db_con = PgSQL.connect(database='cobra')
        db_cursor = db_con.cursor()
        the_genus = 'thermotoga'
        db_cursor.execute('Set search_path to ' + the_genus)
        print 'Building element matrix'
        start_time = time()
        the_elements = create_element_matrix(cobra_model, db_cursor=db_cursor,
                                             me_matrix=True, cas=cas)
        if cas.lower() != 'ginac':
            #Can't pickle PySwigObjects
            with open(element_file, 'w') as out_file:
                dump(the_elements, out_file) 
        print 'Element matrix %s created in %f minutes'%(element_file,
                                                       (time()-start_time)/60)
    element_index = the_elements['elements'][the_element]
    element_vector = the_elements['matrix'][element_index, :]
    print 'This problem is symbolic and may take some time to solve'
    if not lexists(problem_file):
        print 'Constructing the problem for element %s. Restart to solve the problem'%the_element
        print 'This process can take 1-100 minutes depending on model size'
        start_time = time()
        
        the_problem = create_balancing_problem(cobra_model,
                                               element_vector,
                                               system_type=system_type,
                                               cas=cas,
                                               print_unbalanced=True)
        if cas.lower() != 'ginac':
            print 'Rerun the script to solve the problem'
            with open(problem_file, 'w') as out_file:
                dump(the_problem, out_file)
        print 'Problem file %s created in %1.2f minutes.'%(problem_file,
                                                            (time() - start_time) / 60)
    
    ## else:
    ##     with open(problem_file) as in_file:
    ##         the_problem = load(in_file)

    ##     print 'Solving the problem %s'%problem_file
    ##     print 'This may take some time'
    ##     start_time = time()
    ##     the_solution = solve_balance_problem(the_problem['equations'],
    ##                                          the_problem['variables'])

    ##     with open(solution_file, 'w') as out_file:
    ##         dump(the_solution, out_file)
    ##     print 'Problem solved. %s created in %1.2f minutes.'%(solution_file,
    ##                                                           (time() - start_time) / 60) 
    
#solve nonsymbolic problems
        ## element_matrix = the_elements['matrix'].tocsr()
        ## elements = the_elements['elements']
        ## reaction_matrix = cobra_model._S.T.tocsr()
        ## print 'Not a symbolic problem'
        ## the_balance = element_matrix * reaction_matrix
        ## the_balance = the_balance.tolil()
        ## for e, i in elements:
        ##     print '%s generated by system: %f'%(e, the_balance[i,:].sum())
