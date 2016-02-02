import re


def parseReactionFormula(formula=None):
    """Parse formula into a list of metabolites and a list of coefficients
    
     [metaboliteList,stoichCoeffList,revFlag] = parseReactionFormula(formula)
    
    INPUT
    formula           Reaction formula, may contain symbols '+', '->', '<=>' in
                       addition to stoichiometric coefficients and metabolite names
                       examples: 
                       '0.01 cdpdag-SC[m] + 0.01 pg-SC[m]  -> 0.01 clpn-SC[m] + cmp[m] + h[m]' (irreversible reaction)
                       'cit[c] + icit[x]  <=> cit[x] + icit[c] ' (reversible reaction)
                       If no stoichiometric coefficient is provided, it is assumed
                       to be = 1
                       Reaction formula should be a string, not a cell array
                        (default: irreversible reaction above)
    
    OUTPUTS
    metaboliteList    List of metabolite names
    compartmentList   List of compartments 
    stoichCoeffList   List of S coefficients
    revFlag           Indicates whether the reaction is reversible (True) or
                       not (False)
    
    Example:
         
    formula = '0.01 cdpdag-SC[m] + 0.01 pg-SC[m]  -> 0.01 clpn-SC[m] + cmp[m] + h[m]'
    
    [metaboliteList,stoichCoeffList,revFlag] = parseReactionFormula(formula)
    
        metaboliteList = 
          ['cdpdag-SC[m]','pg-SC[m]','clpn-SC[m]','cmp[m]','h[m]']
        compartmentList = 
          ['m','m','m','m','m']
        stoichCoeffList = 
           [-0.01,-0.01,0.01,1,1]
        revFlag =
           False
    
     Markus Herrgard 6/1/07
     Richard Que 1/25/10 Modified to handle '-->' and '<==>' as arrows 
     as well as reactionsformatted as '[compartment] : A --> C'. 
     IT May 2012 Modified to handle '=>'
     WBryant Jan 2016 Converted to Python for COBRApy
     """
    
    formula = formula or '0.01 cdpdag-SC[m] + 0.01 pg-SC[m]  -> 0.01 clpn-SC[m] + cmp[m] + h[m]'
    
    tokens = formula.split()
    
    metaboliteList = []
    compartmentList = []
    stoichCoeffList = []
    
    revFlag = True
    
    # Marks the start of a new stoichiometry + metabolite block
    newMetFlag = True
    # Designates products vs reactants
    productFlag = False
    compartment_default = 'c'
    forwardArrows = ['->','-->','=>','==>']
    reversibleArrows = ['<=>','<==>']
    reverseArrows = ['<-','<--','<=','<==']
    switchDirection = False
    for token in tokens:
        t = token
        if t.startswith('['):
            #set compartment
            try:
                compartment_default = re.search('\[(.+)\]',t).group(1)
            except:
                pass
        elif t == ':':
            pass
        elif t == '+':
            newMetFlag = True
        elif t in forwardArrows:
            # Irreversible
            revFlag = False
            productFlag = True
            newMetFlag = True
        elif t in reversibleArrows:
            # Reversible
            revFlag = True
            productFlag = True
            newMetFlag = True
        elif t in reverseArrows:
            revFlag = False
            productFlag = True
            newMetFlag = True
            switchDirection = True            
        else:
            try:
                # Coefficient   
                sCoeff = float(t)
                if not productFlag:
                    sCoeff = -sCoeff
                if switchDirection:
                    sCoeff = -sCoeff
                stoichCoeffList.append(sCoeff)
                newMetFlag = False
            except:
                # Extract compartment if possible
                metaboliteID = t
                try:
                    compartment = re.search('.+\[(.+)\]$',t).group(1)
                    metaboliteID = re.sub('\[(.+)\]$',r'_\1',t)
                except:
                    try:
                        compartment = re.search('.+_([^_]+)$').group(1)
                    except:
                        compartment = compartment_default
                        metaboliteID += '_' + compartment
                compartmentList.append(compartment)
                
                metaboliteList.append(metaboliteID)
                
                if newMetFlag:
                    if not productFlag:
                        sCoeff = -1
                    else:
                        sCoeff = 1
                    if switchDirection:
                        sCoeff = -sCoeff
                    stoichCoeffList.append(sCoeff)
                    newMetFlag = True
    
    return metaboliteList,compartmentList,stoichCoeffList,revFlag
