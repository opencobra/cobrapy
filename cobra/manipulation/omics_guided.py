#cobra.manipulation.omics_guided.py
#The functions for omics_guided tailoring will be kept here.
def tailor_model(cobra_model, the_method='GIMME', data_type='mRNA', data_kind='log_ratio',
                 solver='glpk', the_problem='return' ):
    """

    the_method: Type of tailoring to employ.  GIMME or shlomi.
    data_type: 'mRNA', 'protein', 'metabolite', ...
    data_kind: 'p-value','log_ratio': assumed vs control, 'intensity'
    solver: 'glpk' or 'gurobi'
    
    
    """
    cobra_model = cobra_model.copy()
    print 'Under development'
    return

#function [reactionActivity,reactionActivityIrrev,model2gimme,gimmeSolution] = solveGimme(model,objectiveCol,expressionCol,cutoff)
