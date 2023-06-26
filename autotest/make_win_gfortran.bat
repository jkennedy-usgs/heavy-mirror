set PATH=%PATH%;C:\MinGW\bin
gfortran --version
python make_gfortran.py -fc gfortran -sd -mc ..\src heavy.exe 
pause