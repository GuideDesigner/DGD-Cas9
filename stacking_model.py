# Stacking model of a pair of RNA-RNA base pairs
# return a dictionary of dictionary for all possible pairs with the following format:
# free_energy["AC"]["UG"]= (free_energy_of_AC/UG_base_pair)
def return_stacking_model():
    """
    Return a dictionary of dictionary for all possible pairs with the following format:
    free_energy["AC"]["UG"]= (free_energy_of_AC/UG_base_pair)
    """
    stacking_model = {}

    stacking_model["AA"] = {}
    stacking_model["AC"] = {}
    stacking_model["AG"] = {}
    stacking_model["AU"] = {}

    stacking_model["CA"] = {}
    stacking_model["CC"] = {}
    stacking_model["CG"] = {}
    stacking_model["CU"] = {}

    stacking_model["GA"] = {}
    stacking_model["GC"] = {}
    stacking_model["GG"] = {}
    stacking_model["GU"] = {}

    stacking_model["UA"] = {}
    stacking_model["UC"] = {}
    stacking_model["UG"] = {}
    stacking_model["UU"] = {}   

    # watson-crick base pairs
    stacking_model["AA"]["UU"] = -0.94
    stacking_model["AC"]["UG"] = -2.25
    stacking_model["AG"]["UC"] = -2.01
    stacking_model["AU"]["UA"] = -1.09   
    stacking_model["CA"]["GU"] = -2.07
    stacking_model["CC"]["GG"] = -3.28
    stacking_model["CG"]["GC"] = -2.33
    stacking_model["CU"]["GA"] = -2.01
    stacking_model["GA"]["CU"] = -2.42
    stacking_model["GC"]["CG"] = -3.46
    stacking_model["GG"]["CC"] = -3.28
    stacking_model["GU"]["CA"] = -2.25           
    stacking_model["UA"]["AU"] = -1.29
    stacking_model["UC"]["AG"] = -2.42
    stacking_model["UG"]["AC"] = -2.07
    stacking_model["UU"]["AA"] = -0.94

    # non-watson-crick base pairs
    stacking_model["GC"]["UG"] = -2.23   
    stacking_model["CU"]["GG"] = -1.93
    stacking_model["GG"]["CU"] = -1.80
    stacking_model["CG"]["GU"] = -1.05
    stacking_model["AU"]["UG"] = -0.76
    stacking_model["GA"]["UU"] = -0.60
    stacking_model["UG"]["GU"] = -0.38
    stacking_model["UA"]["GU"] = -0.22
    stacking_model["GG"]["UU"] = -0.20
    stacking_model["GU"]["UG"] = -0.19
    stacking_model["AG"]["UU"] = +0.02

    stacking_model["GU"]["CG"] = -2.23   
    stacking_model["GG"]["UC"] = -1.93
    stacking_model["UC"]["GG"] = -1.80
    stacking_model["UG"]["GC"] = -1.05
    stacking_model["GU"]["UA"] = -0.76
    stacking_model["UU"]["AG"] = -0.60
    stacking_model["UG"]["GU"] = -0.38
    stacking_model["UG"]["AU"] = -0.22
    stacking_model["UU"]["GG"] = -0.20
    stacking_model["GU"]["UG"] = -0.19
    stacking_model["UU"]["GA"] = +0.02   

    return stacking_model