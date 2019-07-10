function [M,ra,offset]=get_model_spectra(a)
% [M,ra,offset]=get_model_spectra(a)
% offset frequency of the center of the model spectra relative to water in Hz.
if ~exist('a','var')
   a={'Cr','GABA','Glu','Gln','Lac','Ins','NAA','pCho','Tau','pCr','Mac'};
end
nm=length(a);
M=zeros(2043,nm);
ra=zeros(1,nm);
offset=zeros(1,nm);
for i=1:nm
     root=fileparts(fileparts(mfilename('fullpath')));
    
    switch a{i}
        case 'Ins'
         fid=read_mrui(fullfile(root,'MeasuredModelSpectra/Ins_Phantom_052612/sp1.fid/fid_model_sw4000.mrui'));
         offset(i)=0;
        case 'pCho'
            
         [hd,fid]=read_mrui(fullfile(root,'spectrum_simulation/sw4000',a{i}));
         offset(i)=-730;
        case 'Mac'
        % fid=read_mrui('c:/labhome/backup/matwork/LCModel/MeasuredModelSpectra/Mac.mrui');   %data from 1118_MCAO
        % offset(i)=-730;     
          fid=read_mrui(fullfile(root,'MeasuredModelSpectra/MacMol_112712/sum_mac234.mrui'));   %data from 1118_MCAO
          offset(i)=-760;     
        otherwise
         [hd,fid]=read_mrui(fullfile(root,'spectrum_simulation/sw4000',a{i}));
         offset(i)=-670;
    end
    
  fa=fft(fid);
  fa=fftshift(fa);
  ra(i)=abs(fid(1,:));
  ph=phase(fid(1));
  M(:,i)=fa.*exp(-1i*ph);
end

M=M./repmat(ra,[size(M,1),1]);  %normalize to avoid numerical errors in error calculation for Macromolecules.
ra=ones(1,nm);



