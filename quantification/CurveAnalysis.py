"""
This submodule contains functions using data fiting.
"""
import numpy as np
import cmath
from sklearn.linear_model import LinearRegression
from sklearn.mixture import GaussianMixture
from ..errors import SolutionNotRealError

def simple_linear_regression(X: np.array, Y: np.array) :
    """
    from X, Y distribution returns (slope, intercep) from linear regression curve.
    """
    X = np.array(X).reshape(-1,1)
    Y = np.array(Y)
    lin_model = LinearRegression()
    lin_model.fit(X,Y)

    return lin_model.coef_[0], lin_model.intercept_

def _MultiGaussianfit(distribution:'list[float]') :
    """
    Fit a distribution with bi-modal gaussian curve.

    Returns
    -------
        res : dict
            'mu1', 'mu2', 'sigma1', 'sigma2'
            mu is the expected value and sigma the variance of the individual gaussian distribution.
    """
    distribution = np.array([[distribution],[distribution]], dtype= float).reshape(2,-1)
    Gaussian_fit = GaussianMixture(n_components= 2).fit(distribution)
    mu1, mu2 = Gaussian_fit.means_[0], Gaussian_fit.means_[1]
    sigma1, sigma2 = np.sqrt(Gaussian_fit.covariances_[0], Gaussian_fit.covariances[1])
    res = {'mu1' : mu1,
           'mu2' : mu2,
           'sigma1' : sigma1,
           'signa2' : sigma2
           }
    
    return res

def _Guassians_intersect(mu1:float, mu2:float, sigma1:float, sigma2:float) -> float :
    """
    Finds the x-axis coordinate where 2 gaussians intersect. This can be achieved by solving a 2nd degree equation ax² + bx + c = 0 where a,b and c are defined as below.
    """
    a = np.power(sigma2,2) - np.power(sigma1,2)
    b = 2*(np.power(sigma1,2)*mu2 - np.power(sigma2,2)*mu1)
    c = np.power(sigma2,2)*np.power(mu1,2) - np.power(sigma1,2)*np.power(mu2,2) - 2*np.power(sigma1,2)*np.power(sigma2,2)*np.log(sigma2/sigma1)

    ans1, ans2 = solve_quadratic_equation(a,b,c)
    return ans2

def solve_quadratic_equation(a,b,c, real= False) :
    """
    Solve a quadratic equation ax² + bx + c = 0

    Returns
    -------
    res = (ans1,ans2) 
    """

    # calculating  the discriminant
    dis = (b**2) - (4 * a*c)
    if real and dis < 0 : raise SolutionNotRealError('Equation is set to real set but discriminant was found < 0 which means there are no solutions. Try setting "real" parameter to False.')
    
    # find two results
    ans1 = (-b-cmath.sqrt(dis))/(2 * a)
    ans2 = (-b + cmath.sqrt(dis))/(2 * a)

    return (ans1,ans2)


def gaussian(x,mu,sigma) -> float: 
    res = (1/sigma*np.sqrt(2*np.pi))*np.exp(-np.power((x-mu),2)/(2*np.power(sigma,2)))
    return res


def multi_gaussian(x,mu_list:list, sigma_list:list) -> float :
    
    if isinstance(mu_list, (list, tuple, np.ndarray)) and isinstance(sigma_list, (list, tuple, np.ndarray)):
        if len(mu_list) != len(sigma_list) : raise ValueError("mu and sigma list must have the same length.")
        is_iterable = True
    elif isinstance(mu_list, (int,float)) and isinstance(sigma_list, (int, float)) :
        is_iterable = False
    else :
        raise TypeError("mu and sigma parameters should either be both float like, or iterables (list,tuple,np.ndarray).")

    if is_iterable :
        res = 0
        for mu,sigma in zip(mu_list, sigma_list):
            res += gaussian(x,mu,sigma)
    else : 
        mu, sigma = mu_list, sigma_list
        res = gaussian(x,mu,sigma)
    
    return res


def multi_gaussian_fit(x, mu_list, sigma_list, data_point) : 
    res = multi_gaussian(x= x, mu_list=mu_list, sigma_list=sigma_list)
    return res - data_point

