clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Extrapole and interpole biology data and write in climatology file
%
%  Data source : IRI/LDEO Climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA2001/
%
%  Pierrick Penven, IRD, 2005.
%  Olivier Aumont the master, IRD, 2006.
%  Patricio Marchesiello, chief, IRD, 2007.
%  Christophe Eugene Raoul Menkes, the slave, IRD, 2007.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
%  Switches for selecting what to process (1=ON)
%
romstools_param
%
% Climatology files
%
no3_seas_data  = '../WOAPISCES/no3_seas.cdf';
no3_ann_data    = '../WOAPISCES/no3_ann.cdf';
po4_seas_data  = '../WOAPISCES/po4_seas.cdf';
po4_ann_data    = '../WOAPISCES/po4_ann.cdf';
o2_seas_data   = '../WOAPISCES/o2_seas.cdf';
o2_ann_data     = '../WOAPISCES/o2_ann.cdf';
sio3_seas_data = '../WOAPISCES/sio3_seas.cdf';
sio3_ann_data   = '../WOAPISCES/sio3_ann.cdf';
dic_seas_data  = '../WOAPISCES/dic_seas.cdf';
dic_ann_data    = '../WOAPISCES/dic_ann.cdf';
talk_seas_data = '../WOAPISCES/talk_seas.cdf';
talk_ann_data   = '../WOAPISCES/talk_ann.cdf';
doc_seas_data  = '../WOAPISCES/doc_seas.cdf';
doc_ann_data    = '../WOAPISCES/doc_ann.cdf';
fer_seas_data  = '../WOAPISCES/fer_seas.cdf';
fer_ann_data    = '../WOAPISCES/fer_ann.cdf';
dust_seas_data = '../WOAPISCES/dust_seas.cdf';
dust_ann_data   = '../WOAPISCES/dust_ann.cdf';
%
cycle=woa_cycle;
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Add variables in the files
%
add_no3(oaname,clmname,ininame,grdname,no3_seas_data,...
        no3_ann_data,cycle,makeoa,makeclim)

add_po4(oaname,clmname,ininame,grdname,po4_seas_data,...
        po4_ann_data,cycle,makeoa,makeclim)

add_sio3(oaname,clmname,ininame,grdname,sio3_seas_data,...
        sio3_ann_data,cycle,makeoa,makeclim)

add_o2(oaname,clmname,ininame,grdname,o2_seas_data,...
        o2_ann_data,cycle,makeoa,makeclim)

add_dic(oaname,clmname,ininame,grdname,dic_seas_data,...
       dic_ann_data,cycle,makeoa,makeclim)

add_talk(oaname,clmname,ininame,grdname,talk_seas_data,...
       talk_ann_data,cycle,makeoa,makeclim)

add_doc(oaname,clmname,ininame,grdname,doc_seas_data,...
       doc_ann_data,cycle,makeoa,makeclim)
            
add_fer(oaname,clmname,ininame,grdname,fer_seas_data,...
       fer_ann_data,cycle,makeoa,makeclim)

%
% Horizontal extrapolation of NO3
%
if (makeoa)
  ext_tracers(oaname,no3_seas_data,no3_ann_data,...
              'nitrate','NO3','no3_time','Zno3',Roa);

  ext_tracers(oaname,po4_seas_data,po4_ann_data,...
              'phosphate','PO4','po4_time','Zpo4',Roa);

  ext_tracers(oaname,sio3_seas_data,sio3_ann_data,...
              'silicate','Si','si_time','Zsi',Roa);

  ext_tracers(oaname,o2_seas_data,o2_ann_data,...
              'oxygen','O2','o2_time','Zo2',Roa);

  ext_tracers(oaname,dic_seas_data,dic_ann_data,...
             'dic','DIC','dic_time','Zdic',Roa);

  ext_tracers(oaname,talk_seas_data,talk_ann_data,...
              'talk','TALK','talk_time','Ztalk',Roa);

  ext_tracers(oaname,doc_seas_data,doc_ann_data,...
             'doc','DOC','doc_time','Zdoc',Roa);

  ext_tracers(oaname,fer_seas_data,fer_ann_data,...
             'fer','FER','fer_time','Zfer',Roa);
end
%
% Vertical interpolations 
%
if (makeclim)
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' NO3...')
  vinterp_clm(clmname,grdname,oaname,'NO3','no3_time','Zno3',0,'r');
%
%  PO4
%
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' PO4...')
  vinterp_clm(clmname,grdname,oaname,'PO4','po4_time','Zpo4',0,'r');
 %
 %  Si
 %
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' Si...')
  vinterp_clm(clmname,grdname,oaname,'Si','si_time','Zsi',0,'r');
%
%  O2
%
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' O2...')
  vinterp_clm(clmname,grdname,oaname,'O2','o2_time','Zo2',0,'r');

%
% DIC 
%
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' DIC...')
  vinterp_clm(clmname,grdname,oaname,'DIC','dic_time','Zdic',0,'r');
%
% TALK
%
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' TALK...')
  vinterp_clm(clmname,grdname,oaname,'TALK','talk_time','Ztalk',0,'r');

%
% DOC
%
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' DOC...')
  vinterp_clm(clmname,grdname,oaname,'DOC','doc_time','Zdoc',0,'r');
%
% FER
%
  disp(' ')
  disp(' Vertical interpolations')
  disp(' ')
  disp(' FER...')
  vinterp_clm(clmname,grdname,oaname,'FER','fer_time','Zfer',0,'r');
end
%
% Initial file
%
if (makeini)
  disp(' ')
  disp(' Initial')
  disp(' ')
  disp(' NO3...')
  add_ini_no3(ininame,grdname,oaname,cycle);
%
%  PO4
%
  disp(' ')
  disp(' Initial')
  disp(' ')
  disp(' PO4...')
  add_ini_po4(ininame,grdname,oaname,cycle);
%
%  Si
%
  disp(' ')
  disp(' Initial')
  disp(' ')
  disp(' Si...')
  add_ini_sio3(ininame,grdname,oaname,cycle);
%
%  O2
%
 disp(' ')
 disp(' Initial')
 disp(' ')
 disp(' O2...')
 add_ini_o2(ininame,grdname,oaname,cycle);
%
%  O2
%
 disp(' ')
 disp(' Initial')
 disp(' ')
 disp(' DIC...')
 add_ini_dic(ininame,grdname,oaname,cycle);
%
%  TALK
%
 disp(' ')
 disp(' Initial')
 disp(' ')
 disp(' TALK...')
 add_ini_talk(ininame,grdname,oaname,cycle);
%
%  DOC
%
 disp(' ')
 disp(' Initial')
 disp(' ')
 disp(' DOC...')
 add_ini_doc(ininame,grdname,oaname,cycle);
%
%  TALK
%
 disp(' ')
 disp(' Initial')
 disp(' ')
 disp(' FER...')
 add_ini_fer(ininame,grdname,oaname,cycle);
end
%
% Make a few plots
%
disp(' ')
disp(' Make a few plots...')
test_clim(clmname,grdname,'NO3',1,coastfileplot)
figure
test_clim(clmname,grdname,'PO4',1,coastfileplot)
figure
test_clim(clmname,grdname,'Si',1,coastfileplot)
figure
test_clim(clmname,grdname,'O2',1,coastfileplot)
figure
test_clim(clmname,grdname,'DIC',1,coastfileplot)
figure
test_clim(clmname,grdname,'TALK',1,coastfileplot)
figure
test_clim(clmname,grdname,'DOC',1,coastfileplot)
figure
test_clim(clmname,grdname,'FER',1,coastfileplot)
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











