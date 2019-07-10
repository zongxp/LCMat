[h,a]=read_mrui('Cr_3ppm');
[h,a2]=read_mrui('Cr_3p9ppm');

f=(1-exp(-3/0.88))/(1-exp(-3/1.34));

write_mrui2(a2*f+a,'Cr_T1corrected');


[h,a]=read_mrui('pCr_3ppm');
[h,a2]=read_mrui('pCr_3p9ppm');

f=(1-exp(-3/0.88))/(1-exp(-3/1.34));

write_mrui2(a2*f+a,'pCr_T1corrected');