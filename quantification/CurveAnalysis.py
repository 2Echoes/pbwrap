"""
This submodule contains functions using data fiting.
"""
import numpy as np
import cmath
from sklearn.linear_model import LinearRegression
from sklearn.mixture import GaussianMixture
from ..errors import SolutionNotRealError

def simple_linear_regression(X: np.array, Y: np.array) :
    X = np.array(X).reshape(-1,1)
    Y = np.array(Y)
    lin_model = LinearRegression()
    lin_model.fit(X,Y)

    return lin_model.coef_[0], lin_model.intercept_


def _MultiGaussianfit(distribution:'list[float]', gaussian_number=2) :
    """
    Fit a distribution with 'gaussian_number'-modal gaussian curve.

    Returns
    -------
        res : dict
            'mu1', 'mu2', 'sigma1', 'sigma2'
            mu is the expected value and sigma the variance of the individual gaussian distribution.
    """

    Gaussian_fit = GaussianMixture(n_components= gaussian_number).fit(distribution)
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