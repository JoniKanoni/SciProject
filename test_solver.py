"""
This module contains tests for the solver module
"""
import numpy as np
import numpy.random as random
import solver
 
def test_Int():
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xvals = np.linspace(border1,border2,2000)
    a, b, c = random.randint(1,10), random.randint(1,10), random.randint(1,10)
    poly = a + b * xvals + c * np.square(xvals)
    ana_int = a * (border2 - border1) + b * (border2 **2 - border1 **2) / 2 + c * (border2 **3 - border1 ** 3) / 3
    num_int = solver.Integrate(poly,xvals)
    assert np.abs((ana_int - num_int) / ana_int) < 0.01
     
def test_Norm():
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xvals = np.linspace(border1,border2,2000)
    a, b, c = random.randint(1,10), random.randint(1,10), random.randint(1,10)
    poly = a + b * xvals + c * np.square(xvals)
    ana_int_a =a * ( a * (border2 - border1) + b * (border2 **2 - border1 **2) / 2 + c * (border2 **3 - border1 **3) / 3)
    ana_int_b =b * ( a * (border2 **2 - border1 **2) / 2 + b * (border2 **3 - border1 **3) / 3 + c * (border2 **4 - border1 **4) / 4)
    ana_int_c =c * ( a * (border2 **3 - border1 **3) / 3 + b * (border2 **4 - border1 **4) / 4 + c * (border2 **5 - border1 **5) / 5)
    ana_int = ana_int_a + ana_int_b + ana_int_c
    Normed_ana = poly / np.sqrt(ana_int)
    Normed_num = solver.QM_Norming(poly, xvals)
    assert np.all(np.abs((Normed_ana - Normed_num) / Normed_ana) < 0.01)
     
def test_Pos_info():
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xvals = np.linspace(border1,border2,2000)
    a, b, c = random.randint(1,10), random.randint(1,10), random.randint(1,10)
    poly = a + b * xvals + c * np.square(xvals)
    ana_exp1_a =a * ( a * (border2 **2 - border1 **2) / 2 + b * (border2 **3 - border1 **3) / 3 + c * (border2 **4 - border1 **4) / 4)
    ana_exp1_b =b * ( a * (border2 **3 - border1 **3) / 3 + b * (border2 **4 - border1 **4) / 4 + c * (border2 **5 - border1 **5) / 5)
    ana_exp1_c =c * ( a * (border2 **4 - border1 **4) / 4 + b * (border2 **5 - border1 **5) / 5 + c * (border2 **6 - border1 **6) / 6)
    ana_exp2_a =a * ( a * (border2 **3 - border1 **3) / 3 + b * (border2 **4 - border1 **4) / 4 + c * (border2 **5 - border1 **5) / 5)
    ana_exp2_b =b * ( a * (border2 **4 - border1 **4) / 4 + b * (border2 **5 - border1 **5) / 5 + c * (border2 **6 - border1 **6) / 6)
    ana_exp2_c =c * ( a * (border2 **5 - border1 **5) / 5 + b * (border2 **6 - border1 **6) / 6 + c * (border2 **7 - border1 **7) / 7)
    ana_int_a =a * ( a * (border2 - border1) + b * (border2 **2 - border1 **2) / 2 + c * (border2 **3 - border1 **3) / 3)
    ana_int_b =b * ( a * (border2 **2 - border1 **2) / 2 + b * (border2 **3 - border1 **3) / 3 + c * (border2 **4 - border1 **4) / 4)
    ana_int_c =c * ( a * (border2 **3 - border1 **3) / 3 + b * (border2 **4 - border1 **4) / 4 + c * (border2 **5 - border1 **5) / 5)
    ana_int = ana_int_a + ana_int_b + ana_int_c
    poly = poly / np.sqrt(ana_int)
    expval1_ana = (ana_exp1_a + ana_exp1_b + ana_exp1_c) / np.square(ana_int)
    expval2_ana = (ana_exp2_a + ana_exp2_b + ana_exp2_c) / np.square(ana_int)
    expval1_num, variance_num = solver.QM_position_info(poly, xvals)
    variance_ana = np.sqrt(expval2_ana - np.square(expval1_ana))
    assert np.abs((expval1_ana - expval1_num) / expval1_ana) < 0.01 and np.abs((variance_ana - variance_num) / variance_ana) < 0.01
