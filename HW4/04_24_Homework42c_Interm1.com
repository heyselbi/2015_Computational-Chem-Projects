%NProcShared=12
%Chk=04_24_Homework42c_Interm1.chk
#n RHF/6-31G(d,p) pop=Reg Opt

 04_24_Homework42c_Interm1 

0 1
C  
C   1 B1
H   1 B2 2 A2
H   1 B3 2 A3 3 D3
H   1 B4 2 A4 3 D4
C   2 B5 1 A5 3 D5
H   6 B6 2 A6 1 D6
H   6 B7 2 A7 1 D7
H   6 B8 2 A8 1 D8
C   2 B9 1 A9 3 D9
H   10 B10 2 A10 1 D10
H   10 B11 2 A11 1 D11
H   10 B12 2 A12 1 D12
Cl  4 B13 1 A13 2 D13
Variables:
B1        1.49620
B2        1.09429
A2      110.68314
B3        1.09411
A3      110.48317
D3      239.74445
B4        1.09435
A4      110.84208
D4      119.24581
B5        1.49613
A5      117.22667
D5      196.46640
B6        1.09434
A6      110.83916
D6       44.43151
B7        1.09431
A7      110.68707
D7      163.67723
B8        1.09412
A8      110.48570
D8      283.93676
B9        1.49225
A9      118.27777
D9       44.58469
B10        1.09429
A10      110.79009
D10      315.48224
B11        1.09431
A11      110.78782
D11      196.15530
B12        1.09416
A12      110.47808
D12       75.81962
B13        3.90393
A13      117.99511
D13      270.69110

