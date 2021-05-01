# implementation of MARKOV CHAIN MONTE-CARLO optimization on discrete grid


# libraries
import argparse

import numpy as np


# helpers by staff
def get_cmdline_parser():
    """
    method to help parse info from cmd line
    :return: the parser itself
    """
    parser = argparse.ArgumentParser(
        description='Run MCMC on a discrete n*n configuration space.')
    parser.add_argument('n', type=int, default=5, nargs='?',
                        help='number of rows/columns of grid')
    parser.add_argument('m', type=int, default=1000, nargs='?',
                        help='number of iteration of MCMC optimization')
    parser.add_argument('k', type=float, default=1.0, nargs='?',
                        help='kT - the denominator for the metropolis '
                             'criterion, (Boltzmann constant times '
                             'temperature)')
    return parser


def is_valid(c, n):
    """ Return True if c is a valid 2-D coordinate on an n x n grid
   with 0-based indices """
    if len(c) != 2:
        return False
    return (c[0] >= 0 and c[1] >= 0 and c[0] < n and c[1] < n)


def get_p_accept_metropolis(dE, kT, p_forward, p_backward):
    """
    return the probability to accept the metropolis criteria
    for the specified conditions
    dE - change in energy from current to proposed configuration
    kT - the factor of Boltzmann constant (kB) and the temperature (T)
    p_forward - probability to propose a move from current to proposed
    configuration p_backward - probability to propose a move from proposed to
    current configuration """

    p = np.exp(-dE / kT) * p_backward / p_forward
    return min(p, 1.0)


def E(c):
    """
    method to get Energy levels of coordinates C
    :param c: tuple of (x,y) on grid
    :return: energy levels
    """
    assert (len(c) == 2)
    return 1.0 * c[0] + 0.5 * c[1]


def get_neighbours(c, n):
    """
    get up/down/left/right neighbours on an n x n grid with 0-based indices
    """
    assert (is_valid(c, n))
    ret_value = []
    if c[0] > 0:
        ret_value.append((c[0] - 1, c[1]))
    if c[0] < n - 1:
        ret_value.append((c[0] + 1, c[1]))
    if c[1] > 0:
        ret_value.append((c[0], c[1] - 1))
    if c[1] < n - 1:
        ret_value.append((c[0], c[1] + 1))
    return ret_value
