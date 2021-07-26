'''
Tests the reader.py
'''
import numpy as np
import reader

def test_unendlich():
    ''' Test 1 - Unendlich tiefer'''
    massetest = 2.0
    xmintest = -2.0
    xmaxtest = 2.0
    npointtest = 1999.0
    firsttest = 1
    lasttest = 5
    interpolationtest = 'linear'
    numinterpolationtest = 2
    resttest = np.array([[-2.0, 0.0], [2.0, 0.0]])
    masse, xmin, xmax, xnum, first, last, interpolation, numinterpol, rest = reader.getdata('sciproject')
    assert masse == massetest
    assert xmin == xmintest
    assert xmax == xmaxtest
    assert xnum == npointtest
    assert first == firsttest
    assert last == lasttest
    assert interpolation == interpolationtest
    assert numinterpol == numinterpolationtest
    assert np.array_equal(rest,resttest)
