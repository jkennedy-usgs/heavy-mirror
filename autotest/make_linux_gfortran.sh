gfortran --version
python make_gfortran.py -fc gfortran -ff='-w' -sd -mc ../src heavy
cp heavy ../test/heavy