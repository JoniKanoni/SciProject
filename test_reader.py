'''
Tests the reader
'''
import numpy as np
from modules import reader

def test_unendlich():
    ''' Test 1 - Unendlich tiefer'''
    massetest = 2.0
    xmintest = -2.0
    xmaxtest = 2.0
    npointtest = 1999.0
    firsttest = 1
    lasttest = 5
    interpolationtest = 'linear'
    resttest = np.array([[-2.0, 0.0], [2.0, 0.0]])
    masse, xmin, xmax, xnum, first, last, int_type, pot_points, __ = reader.get_data('results/infinite_well/')
    assert masse == massetest
    assert xmin == xmintest
    assert xmax == xmaxtest
    assert xnum == npointtest
    assert first == firsttest
    assert last == lasttest
    assert int_type == interpolationtest
    assert np.array_equal(pot_points, resttest)
