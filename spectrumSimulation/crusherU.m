function sig=crusherU(sig,sz)
sig=get_commute(sig,sz);
% crusher gradient will remove components of sig that are anti-symmetry
% with 