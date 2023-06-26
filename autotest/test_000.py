# all test modules should be in a file that starts with the name test.
import os
import shutil
import pandas
import numpy as np

def test_py():
    cd = os.getcwd()
    os.chdir(os.path.join('..','python'))
    print("Heavy test: Python All American Canal model")
    os.system('python model1_run.py')
    os.system('cp sim_py.grav ../test/aac/sim_py.grav')
    os.chdir(cd)
    assert True

def test_aac_model():
    cd = os.getcwd()
    print("Heavy test: All American Canal model")
    os.chdir(os.path.join('..','test','aac'))
    os.system(r'../../autotest/heavy model1_hvy.nam')
    assert os.path.exists("./model1.out")
    os.chdir(cd)

def test_py_vs_fortran():
    cd = os.getcwd()
    print("Heavy test: Comparing Python output to Fortran output")
    os.chdir(os.path.join('..','test','aac'))
    py_data = pandas.read_csv('sim_py.grav')
    f_data = pandas.read_csv('model1.out', delim_whitespace=True, header=0, names=['Station','kper','kstp','totim','g','d2'])
    for station in py_data['Station'].unique():
        dp = py_data[py_data['Station'] == station]['g'].reset_index(drop=True)
        df = f_data[f_data['Station'] == 'G' + station]['g'].reset_index(drop=True)
        mse = np.sqrt(((df - dp) ** 2).mean())
        print(f'mse, {station}: {mse:.3f}')
        assert mse < 1.0
    os.chdir(cd)

def test_abq_model():
    cd = os.getcwd()
    print (cd)
    print("Heavy test: Rio Grande Valley transient model")
    os.chdir(os.path.join('..','test','abq','2013y_parent'))
    os.system(r'../../../autotest/heavy -g 6 tran_hvy.nam')    
    assert os.path.exists("./tran_hvy.out")
    os.chdir(cd)
    assert True

def test_truxton_model():
    cd = os.getcwd()
    print("Heavy test: Truxton Basin Hydrologic Model Xhigh Qhigh scenario")
    os.chdir(os.path.join('..','test','truxton','tbhm_tr_Xhigh_Qhigh'))
    os.system(r'../../../autotest/heavy -g 6 tbhm_tr_Xhigh_Qhigh_hvy.nam')    
    assert os.path.exists("./tbhm_tr_Xhigh_Qhigh.out")
    os.chdir(cd)

def test_setup():
    tempdir = os.path.join('.', 'temp')
    if os.path.isdir(tempdir):
        shutil.rmtree(tempdir)
    os.mkdir(tempdir)
    return


if __name__ == "__main__":
    test_setup()
