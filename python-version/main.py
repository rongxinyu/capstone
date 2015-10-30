import mibian
from numpy import *

def bstar(k, alpha):
    return ((k**(alpha+1) - (k-1)**(alpha+1)) / (alpha+1)) ** (1/alpha)

def BSImpliedVolCall(S0, K, T, r, C):
    return mibian.BS([S0, K, r*100., T*365.], callPrice=C).impliedVolatility

def hybridScheme(params):

    S0, xi, eta, alpha, rho = params

    def covMatrix(n, kappa):
        sigma = zeros((kappa+1, kappa+1))
        sigma[0][0] = 1./n
        for j in xrange(1, kappa+1):
            sigma[0][j] = sigma[j][0] = (j**(alpha+1) - (j-1)**(alpha+1)) / ((alpha+1)* n**(alpha+1))
        for j in xrange(1, kappa+1):
            for k in xrange(1, kappa+1):
                if j == k:
                    sigma[j][k] = (j**(2*alpha+1) - (j-1)**(2*alpha+1)) / ((2*alpha+1)*n**(2*alpha+1))
            else:
                pass
        return sigma 

    def Simulation(n, kappa, T):
        W = array([random.multivariate_normal(mean=[0]*(kappa+1), cov=covMatrix(n, kappa)) for i in xrange(int(n*T))])
        Wperp = random.normal(scale=sqrt(1./n), size=int(n*T))
        Z = rho * W[:,0] + sqrt(1-rho**2) * Wperp

        Gamma = array([(bstar(x+1, alpha)/n) ** alpha for x in xrange(int(n*T))])
        Gamma[:kappa+1] = 0
        Y2 = convolve(Gamma, W[:,0])[:100] # size 199
        Y1 = zeros(int(n*T))
        for i in xrange(int(n*T)):
            Y1[i] = 0
            for k in range(min(i, kappa-1)):
                Y1[i] = Y1[i] + W[i+1-k][k+1]
        Y = Y1 + Y2
        v = xi * exp(eta * sqrt(2*alpha+1)*Y - eta**2/2 * array([(x/n)**(2*alpha+1) for i in xrange(int(n*T))]))
        S = S0 * exp(sum(v**0.5 * Z) - 0.5*sum(v)/n)
        return S

    def MC(N, n, kappa, T):
        return [Simulation(n, kappa, T) for i in xrange(N)]

    return MC

def impvol(k, st, T):
    payoff = (st > exp(k)) * (st - exp(k))
    return BSImpliedVolCall(1, exp(k), T, 0, mean(payoff))

def main():
    params = (1., 0.235**2, 1.9, -.43, -.9)
    finalPrices = hybridScheme(params)(10000, 100, 1, 1.)
    print finalPrices

if __name__ == '__main__':
    main()