function lcbg()

p = parameter('LCModel GUI for Bruker');

p=add(p,'directoryname','scan folder','');


%p=add(p,'filename','mrui file','');  % leave blank if use the fid file in scan folder.
p=add(p,'int','extra phase (deg)',90);

p = add(p,'button','Show Sp.','run_LCModel(params)');

p=add(p,'string','baseline range','');

p = add(p,'directoryname','alt. water folder','');  %leave blank if use fid.refscan in the fid scan folder as reference

p = add(p,'button','Calc Norm. Area','run_LCModel(params)');
p=add(p,'bool','Cr');
p=add(p,'bool','GABA');
p=add(p,'bool','Glu');
p=add(p,'bool','Gln');
p=add(p,'bool','Lac');
p=add(p,'bool','Ins');
p=add(p,'bool','NAA');
p=add(p,'bool','pCho');
p=add(p,'bool','Tau');
p=add(p,'bool','pCr');

p = add(p,'directoryname','model sp. Folder',''); 


p = add(p,'button','Go','run_LCModel(params)');

p = add(p,'button','Close','close');

p = parametergui(p);


