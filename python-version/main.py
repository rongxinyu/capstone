import numpy

def convolve(x, y):
	return numpy.convolve(x, y)

def bstar(k, alpha):
	return ((k**(alpha+1) - (k-1)**(alpha+1)) / (alpha+1)) ** (1/alpha)