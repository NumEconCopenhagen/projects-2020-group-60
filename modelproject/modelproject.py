from scipy import optimize

def solve_for_ss(n,tau,rho,alpha):
    """ solve for the steady state level of capital-per-worker

    Args:
        tau (float): taxation (fraction of income)
        rho (float): propensity to consume
        n (float): population growth rate
        alpha (float): cobb-douglas parameter

    Returns:
        result (RootResults): the solution represented as a RootResults object

    """ 

    # a. define objective function
    f = lambda k: k**alpha
    obj_kss = lambda kss: kss - ((1/((2+rho)*(1+n)))*(1-tau)*(1-alpha)*f(kss)-((1+rho)/(2+rho))*tau*((1-alpha)/alpha)*kss)

    #. b. call root finder
    result = optimize.root_scalar(obj_kss,bracket=[0.1,100],method='bisect')
    
    return result
