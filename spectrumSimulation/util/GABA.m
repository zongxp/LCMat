function s=GABA(freq)
s.nspins=6;

s.pairs=[1,3;1,4; 2,3; 2,4;  3,5;3,6; 4,5; 4,6];

%s.pairs=[2,3;2,3';2',3;2',3';3,4;3,4';3',4;3',4'];

s.J=[5.372,7.127,10.578,6.982,7.755,7.432,6.173,7.933];
s.shift=[3.0128,3.0128,1.8890,1.8890,2.284,2.284];

s.f0=freq;
