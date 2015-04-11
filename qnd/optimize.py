#-*-coding: utf-8-*-

__whatami__ = 'qnd Design optimizers.'
__author__  = 'Danny Goldstein <dgold@berkeley.edu>'

from .core import random_design
import numpy as np

def designopt_mc(point_type, lower, upper, score_func, size=100, niter=1000):
    """Randomly generate designs. Return the highest-scoring design.
    
    Arguments:
    
    point_type -- (design.Point) The class of points to use for the
    designs. 

    lower -- array like (n_features,) The lower boundaries of the
    Design support space.
    
    upper -- array like (n_features,) The upper boundaries of the
    Design support space.
    
    Keyword Arguments:

    size -- (int; default 100) The size of the designs to generate.
    
    score_func -- (default maximin) Function that takes an instance of
    a design and returns a score.
    
    niter -- (int; default 1000) Number of designs to realize. 
    
    Returns:
    
    The highest-scoring design.
    """ 

    high_score = -np.inf

    for i in range(niter):
        design = random_design(point_type, lower, upper, size)
        score = score_func(design)
        if score > high_score:
            best_design = design
            high_score = score
    
    return best_design
