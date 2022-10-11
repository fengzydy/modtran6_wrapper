function [] = modtran6_airs_tape5_make_cld_profile(path_modtran,nlyr,z,p,t,q,co2,o3,resolution,fwhm,cloud,v1,v2,angle,iLoc,emis,zt,ts)
    %%%%%%% Jing Feng (jing.feng@noaa.gov), Oct 10, 2022
    % This program generates tape5 file for the modtran computation
    %% path_modtran: path of modtran6 excutable file
    %% nlyr: number of vertical layer - 1
    %% z: altitude, km
    %% p: pressure, hPa
    %% t: K
    %% q: mixing ratio, g/kg; change unit via card 2c
    %% co2: ppmv
    %% o3: mixing ratio, g/kg
    %% resolution: spectral interval for radiance output, cm-1
    %% fwhm: full width half maximum
    % cloud: structure; include: cloud.qi, ice particle density g/m^3
    %		             .ql, liquid particle density g/m^3
    %			           .z,  cloud height
    %% v1, v2: spectral ranges, cm-1
    %% angle: observer angle 
    % iLoc: string for specifying gas concentration
    % emis: emissivity
    % zt: observer height
    % ts: skin temperature
    % !!!!! please read through MODTRAN 6 instruction file !!!!!!!

    if nargin < 1
    path_modtran = '';
    end
    if nargin<18
      ts = t(1);                    % boundary temperature, > 0 / <= 0
    end
    if nargin>=17 & angle<90
    else
      zt=z(end);
    end
    %%% modroot.in
    modroot = ['mlw'];
    fid = fopen('modroot.in','w');
    fprintf(fid,'%s\n',modroot);
    fclose(fid);
    
    tape5_name = ['modroot, '.tp5''];
    zenith_angle = 180-angle; % 180 - 53; %

    %%% set the values
    % card 1
    card1.modtran = 'H';                 % T(M) / C(K) / F(L)
    card1.speed = 'M';                   % S(blank) / M
    card1.binary =' ';                   % F /T
    card1.lymolc =' ';                   % +
    card1.model = 7;                     % 0-7
    card1.t_best =' ';                   % T
    card1.itype = 2;                     % 1-3
    card1.iemsct =1;                    % 0-3
    card1.imult = 0;                     % 0, +-1
    card1.m1 = 0;                        % T and q
    card1.m2 = 0;                        % H2O
    card1.m3 = 0;                        % O3
    card1.m4 = 0;                        % CH4
    card1.m5 = 0;                        % N2O
    card1.m6 = 0;                        % CO
    card1.mdef = 0;                      % 0/1; for O2, NO, SO2, NO2, NH3, HNO3
    card1.i_rd2c= 1;                     % 0/1
    card1.noprnt = -1;                   % 0 / 1 / -1 /-2
    card1.tptemp = ts;                   % boundary temperature, > 0 / <= 0

    card1.surref = 1-emis;%'      0';            % B(b) / L(l) / >=0(blank) / < 0 
    % card 1a
    card1a.dis = 'F';                    % T(t) / F(f, or blank)
    card1a.disazm = ' ';                 % T(t) / F(f, or blank)
    card1a.disalb = ' ';                 % T/F
    card1a.nstr = 8;                     % 2/4/8/16
    card1a.sfwhm= 0;                      % 0/1/-1
    card1a.co2mx = 380;                  % 0(blank, 330), or another value
    card1a.h2ostr = '+4.0';                 % 0(blank) / gxx / axx / positive factor
    card1a.o3str = 1;                    % 0(blank) / gxx / axx / positive factor
    card1a.c_prof = 0;                   % 0 - 7
    card1a.lsunfl = 'f' ;                  % T/F
    card1a.lbmnam = 4 ;                  % F/4/2
    card1a.lfltnm = 'T' ;                  % T/F
    card1a.h2oaer = ' ';                 % 
    card1a.cdtdir = 'T';
    card1a.solcon = 0;
    card1a.cdastm = '';
    card1a.astmc = '';
    card1a.astmx = '';
    card1a.astmo = '';
    card1a.aerrh = '';
    card1a.nssalb = '';
    % card 1a1 1a2 1a3
    card1a1.usersun = '';               
    if resolution<1
    card1a2.bmname = 'LBL2013';                
    else
    card1a2.bmname = '       ';
    end
    card1a3.filtnm = 'DATA/airs.flt';  
    % card 2
    card2.aplus = '  ';                  % blank / A+
    card2.ihaze = -1;                     % -1, 0, 1-10
    card2.cnovam = ' ';                  % blank / N
    card2.iseasn = 0;                    % 0-2
    card2.aruss = '   ';                 % blank / USS
    card2.ivulcn = 0;                    % 0-8
    card2.icstl = 0;                     % 1-10
    
    if length(cloud.z)>=1
    card2.icld = 4;                      % 0-11 / 18 / 19
    else
    card2.icld = 0;
    end
    
    card2.ivsa = 0;                      % 0 / 1
    card2.vis = 0.0;                     % 0 / positive value
    card2.wss = 0.0;                     % wind speed [ m/s ]
    card2.whh = 0.0;                     % 24-h avg [ m/s ]
    card2.rainrt = 0.0;                  % [ mm/hr ]
    card2.gndalt = 0.0;                  % [ km ]
    % card 2a
    card2a.cthik = -9.0;                 % -9.0 (default)
    card2a.calt = -9.0;
    card2a.cext = -9.0;
    card2a.ncralt = length(cloud.z);     %%% number of layers input
    card2a.ncrspc = 223;                 %%% number of bands
    card2a.cwavln = -9.0;
    card2a.ccolwd = -9.0;
    card2a.ccolip = -9.0;
    card2a.chumid = 500;
    card2a.asymwd = -9.0;
    card2a.asymip = -9.0;
    % card 2c
    card2c.ml = nlyr+1;                        
    card2c.ird1 = 0;%1;%%%%%%%%%%%%%%%%%%%%%%%%%change it after to 0;%%0;                     % 0 / 1
    card2c.ird2 = 0;                     % 0-2
    card2c.hmodel = 'tape5_maker.m';
    card2c.ree = ' ';
    card2c.nmolyc = 0;
    card2c.e_mass = 0;
    card2c.airmwt = 0;
    % card 2c1
    card2c1.wmol = zeros(3,card2c.ml);
    card2c1.jchar = repmat('', [card2c.ml 14]);
    card2c1.jcharx = ' ';
    card2c1.jchary = ' ';
    for i = 1:card2c.ml
      card2c1.zm(i) = z(i);
      card2c1.p(i) = p(i);
      card2c1.t(i) = t(i);
      card2c1.wmol(1,i) = q(i)/4;            % H2O
      card2c1.wmol(2,i) = co2(i);          % CO2
      card2c1.wmol(3,i) = o3(i);           % O3
      
      if nansum(o3)>0
        card2c1.jchar(i,:)=['AACAC',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc]; % 1: A/B/C/1-6/blankAC',iLoc,'C',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc];AC',iLoc,'C',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc]; %
                                          % 2: A/B/1-6/blank
                                          % 3: A-H/1-6/blank 
      else
      card2c1.jchar(i,:)=['AACA',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc]; % 1: A/B/C/1-6/blankAC',iLoc,'C',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc];AC',iLoc,'C',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc]; %
                                          % 2: A/B/1-6/blank
                                          % 3: A-H/1-6/blank
      end 
      if nansum(t)==0
        card2c1.jchar(i,:)=['A',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc]; % 1: A/B/C/1-6/blankAC',iLoc,'C',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc];AC',iLoc,'C',iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc,iLoc]; %
                                          % 2: A/B/1-6/blank
                                          % 3: A-H/1-6/blank 
      end
      card2c1.jcharx(i) = ' ';
    end
    % card 2c2 2c2x
    card2c2.wmol = zeros(9,card2c.ml);   % N2O, CO, CH4, O2, NO, SO2, NO2, NH3, HNO3
    card2c2x.wmolx = zeros(13,card2c.ml);% CFCs (heavy molecular gases)
    card2c2y.wmoly = zeros(card2c.nmolyc);
    card2c2.wmol(1,:)=0;%0.316;
    card2c2.wmol(2,:)=0;%0.15;
    card2c2.wmol(3,:)=0;%1.797;
    card2c2.wmol(4,:)=0;%2.09*1e5;
    card2c2.wmol(5,:)=0;
    card2c2.wmol(6,:)=0;
    % card 2d
    card2d.irge = zeros(4,1); 
    % card 2e1
    if card2.icld >0
    zcld = zeros(card2a.ncralt,1);
    ql   = zeros(card2a.ncralt,1);
    qi   = zeros(card2a.ncralt,1);
    %[zcld, ql, qi] = read_cld;
    %zcld=[z(max(i-1,1)) z(icloud)-0.06 z(icloud)];
    %ql=[0 0 0];
    %qi=[1.5 1.5 0];
    zcld =cloud.z;
    ql   =cloud.ql;
    qi   =cloud.qi;
    icesize = zeros(card2a.ncralt,1);
    for i = 1:card2a.ncralt
      card2e1.zcld(i) = zcld(i);           % 
      card2e1.cld(i) = ql(i);
      card2e1.cldice(i) = qi(i);
      card2e1.rr(i) = 0.0; 
    end
    % card 2e2
    total_cl=sum((zcld(2:end)-zcld(1:end-1)).*(ql(2:end)+ql(1:end-1))/2)
    total_ci=sum((zcld(2:end)-zcld(1:end-1)).*(qi(2:end)+qi(1:end-1))/2)
          if ~isfield(cloud,'re')
            cloud.re=10;
          end
          if  ~isfield(cloud,'habit')
            cloud.habit='solid_column_r050';
          end
          load('/project/rrg-yihuang-ad/jfeng93/model/modtran/mat/cld_v2.mat','freq','sca_liq','abs_liq','band_index','bpr_liq');
          load(['/lustre03/project/6003571/jfeng93/data/ice_optics_yang/',cloud.habit,'.mat'])
          reff_ice=cloud.re;    %% um
          re_ice_id=round(reff_ice);
            %% include: reff, sca, abs, wavlen, band_index
            cloud_spectrum.nbands = length(freq);
            card2e2.extc = zeros(7,cloud_spectrum.nbands);
            card2e2.absc = zeros(7,cloud_spectrum.nbands);
            card2e2.asym = zeros(7,cloud_spectrum.nbands);
            cloud_spectrum.wavelength=freq;
            V_ice= v(re_ice_id);
            reff_liq=10*10^(-6);
            V_liq= reff_liq^3*4/3*pi;
            rau_ice=917;
            rau_liq=1000;
            band_index=1:length(freq);
            cloud_spectrum.liquid_water_extinction=(abs_liq(band_index,2)+sca_liq(band_index,2))/V_liq/rau_liq;
            cloud_spectrum.liquid_water_absorption=abs_liq(band_index,2)/V_liq/rau_liq;
            cloud_spectrum.ice_extinction=(ext_ice(band_index,re_ice_id))/V_ice/rau_ice;%*1000;
            cloud_spectrum.ice_absorption=abs_ice(band_index,re_ice_id)/V_ice/rau_ice;
            
            cloud_spectrum.liquid_water_asymmetry=bpr_liq(band_index,2);%(((1-bpr_liq(band_index))./bpr_liq(band_index)).^(1/3)-1)./(((1-bpr_liq(band_index))./bpr_liq(band_index)).^(1/3)+1);%bpr_liq(band_index);
            cloud_spectrum.ice_asymmetry=asy_ice(band_index,re_ice_id);
        for i = 1:cloud_spectrum.nbands
          card2e2.wavlen(i) = 10^4/cloud_spectrum.wavelength(cloud_spectrum.nbands+1-i);
          card2e2.extc(6,i) = cloud_spectrum.liquid_water_extinction(cloud_spectrum.nbands+1-i);
          card2e2.absc(6,i) = cloud_spectrum.liquid_water_absorption(cloud_spectrum.nbands+1-i);
          card2e2.asym(6,i) = cloud_spectrum.liquid_water_asymmetry(cloud_spectrum.nbands+1-i);
          card2e2.extc(7,i) = cloud_spectrum.ice_extinction(cloud_spectrum.nbands+1-i);
          card2e2.absc(7,i) = cloud_spectrum.ice_absorption(cloud_spectrum.nbands+1-i);
          card2e2.asym(7,i) = cloud_spectrum.ice_asymmetry(cloud_spectrum.nbands+1-i);
        end
        card2a.ncrspc=length(band_index);
end
    % card 3
    card3.h1 = zt;                    % [ km ]
    card3.h2 = z(1);       
    re=6371.23;
    if str2num(iLoc)==4 || str2num(iLoc)==5
        re=6356.91;
    end
    if angle>=80
    card3.h2=max([0 (re+card3.h1)*sin(angle/180*pi)-re]);
    end
    if angle>90
      angle
      card3.h2=z(end);
      card3.h1=z(1);
    end             
    card3.angle = zenith_angle;          % measured @ h1
    card3.range = 0.0;                   % [ km ]
    card3.beta = 0.0;                    % 
    card3.ro = 0.0;                      % [ km ]
    card3.lenn = 0;                      % 0 / 1 
    card3.phi = 0.0;                     % 
    % card 3a13
    card3.IPARM=0;
    card3.IPH=0;
    card3.IDAY=63;
    card3.ISOURC=0;
    % card 4
    card4.v1 = v1;                      % [ cm-1 ]
    card4.v2 = v2;                      % 
    card4.dv = resolution;               %
    card4.fwhm = fwhm;             %
    card4.yflag = 'R';                   % T / R
    card4.xflag = 'W';                   % W / M / N
    card4.dlimit = '        ';
    card4.flags = '       ';
    card4.mlflx = '   ';
    card4.vrfrac = '          ';
    % card 5
    card5.irpt = 0;                      % 0 / +-1 / +-3 /+-4
    
    %%% write to file
    filename = [path_modtran, tape5_name];
    fid = fopen(filename,'w');
    % card 1
    fprintf(fid,'%1c%1c%1c%1c%1d%1c%4d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d %4d%8.3f%7f\n',...
      card1.modtran, card1.speed, card1.binary,card1.lymolc,card1.model,card1.t_best,...
      card1.itype, card1.iemsct, card1.imult, ...
      card1.m1, card1.m2, card1.m3, card1.m4, card1.m5, card1.m6, card1.mdef, ...
      card1.i_rd2c, card1.noprnt, card1.tptemp, card1.surref);
    % card 1a
    fprintf(fid,'%1c%1c%1c%2d %4.0f%10.3f%10s%10f%1d%1c %1d %1c %1c %1c%10.3f%1c%9.0f%10.3f%10.3f%10s%10s\n',...
      card1a.dis, card1a.disazm, card1a.disalb,card1a.nstr, card1a.sfwhm, card1a.co2mx, ...
      card1a.h2ostr, card1a.o3str, card1a.c_prof, card1a.lsunfl, card1a.lbmnam, card1a.lfltnm,...
      card1a.h2oaer, card1a.cdtdir, card1a.solcon, card1a.cdastm, card1a.astmc, card1a.astmx,...
      card1a.astmo, card1a.aerrh, card1a.nssalb); 
    % card 1a2
    if (card1a.lbmnam == 4)
    fprintf(fid,'%40s\n',card1a2.bmname);
    end
    fprintf(fid,'%40s\n',card1a3.filtnm);
    fprintf(fid,'%40s\n','/project/rrg-yihuang-ad/jfeng93/model/MODTRAN6.0/DATA/');
    % card 2
    fprintf(fid,'%2s%3d%1c%4d%3s%2d%5d%5d%5d%10.5f%10.5f%10.5f%10.5f%10.5f\n',...
      card2.aplus, card2.ihaze, card2.cnovam, card2.iseasn, card2.aruss, card2.ivulcn, ...
      card2.icstl, card2.icld, card2.ivsa, card2.vis, card2.wss, card2.whh, card2.rainrt, card2.gndalt);
    % card 2a
    if (card2.icld > 0) && (card2.icld <= 10 )
    fprintf(fid,'%8.3f%8.3f%8.3f%4d%4d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n',...
      card2a.cthik, card2a.calt, card2a.cext, card2a.ncralt, card2a.ncrspc, card2a.cwavln, ...
      card2a.ccolwd, card2a.ccolip, card2a.chumid, card2a.asymwd, card2a.asymip);
    end 
    if (card2.icld ==18) || (card2.icld == 19 )
    fprintf(fid,'%8.3f%8.3f%8.3f\n',...
      card2a.cthik, card2a.calt, card2a.cext);
    end     
    % card 2c, 2c1, 2c2, 2c2x
    if ((card1.model == 7) || (card1.model == 0) ||(card1.model == 8)  && (card1.i_rd2c == 1))
    fprintf(fid,'%5d%5d%5d%20s%10s%5d%10.3f%10.3f\n',...
      card2c.ml, card2c.ird1, card2c.ird2, card2c.hmodel,card2c.ree,card2c.nmolyc,card2c.e_mass,card2c.airmwt);
    for i = 1:card2c.ml
    fprintf(fid,'%10.3f%10.5f%10.3f%10.6f%10.3f%10.6f%14s %1c%1c\n',...
      card2c1.zm(i), card2c1.p(i), card2c1.t(i), card2c1.wmol(1,i), card2c1.wmol(2,i), card2c1.wmol(3,i), ...
      card2c1.jchar(i,:), card2c1.jcharx(i),card2c1.jchary );
    if (card2c.ird1 == 1)
    fprintf(fid,'%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n%10.3f\n',...
      card2c2.wmol(1,i), card2c2.wmol(2,i), card2c2.wmol(3,i), ...
      card2c2.wmol(4,i), card2c2.wmol(5,i), card2c2.wmol(6,i), ...
      card2c2.wmol(7,i), card2c2.wmol(8,i), card2c2.wmol(9,i) );
    end
    end
    if (card2c.ird1 == 1) && (card1.mdef == 2)
    fprintf(fid,'%10.0e%10.0e%10.0e%10.0e%10.0e%10.0e%10.0e%10.0e\n%10.0e%10.0e%10.0e%10.0e%10.0e\n',...
      card2c2x.wmolx(1,i), card2c2x.wmol(2,i), card2c2x.wmol(3,i), ...
      card2c2x.wmolx(4,i), card2c2x.wmol(5,i), card2c2x.wmol(6,i), ...
      card2c2x.wmolx(7,i), card2c2x.wmol(8,i), card2c2x.wmol(9,i), ...
      card2c2x.wmolx(10,i), card2c2x.wmol(11,i), card2c2x.wmol(12,i), ...
      card2c2x.wmolx(13,i) );
    end 
    if (card2c.nmolyc >0 && card2c.ird1==1)
    fprintf(fid,'%10.0e%10.0e%10.0e%10.0e%10.0e%10.0e%10.0e%10.0e\n',...
      card2c2y.wmoly(1:card2c.nmolyc,i));
    end
    end % 2c, 2c1, ...
    % card 2d
    if (card2.ihaze == 7) || (card2.icld == 11)
     clofprintf(fid,'%5d%5d%5d%5d\n', ...
      card2d.ireg(1), card2d.ireg(2), card2d.ireg(3), card2d.ireg(4) );
    end
    % card 2e1
    
    if ( (card2.icld >= 1) && (card2.icld <= 10) ) && (card2a.ncralt >= 3) 
    for i = 1:card2a.ncralt
    fprintf(fid,'%10.3f%10.5f%10.5f%10.5f\n',...
      card2e1.zcld(i), card2e1.cld(i), card2e1.cldice(i), card2e1.rr(i) );
    end
    end
    % card 2e2
    if ( (card2.icld >= 1) && (card2.icld <= 10) ) && (card2a.ncrspc >= 2)
    for i = 1:card2a.ncrspc
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',...
      card2e2.wavlen(i), card2e2.extc(6,i), card2e2.absc(6,i), card2e2.asym(6,i), ...
      card2e2.extc(7,i), card2e2.absc(7,i), card2e2.asym(7,i) );
    end
    end
    % card 3
    fprintf(fid,'%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f     %10.3f%10.3f\n',...
      card3.h1, card3.h2, card3.angle, card3.range, card3.beta, card3.ro, ...
      card3.lenn, card3.phi);
    % card 4
    fprintf(fid,'%10.0f%10.0f%10.3f%10.3f%1c%1c%8s%7s%3s%10s\n',...
      card4.v1, card4.v2, card4.dv, card4.fwhm, card4.yflag, card4.xflag, ...
      card4.dlimit, card4.flags,card4.mlflx,card4.vrfrac);
    % card 5
    fprintf(fid,'%5d',...
      card5.irpt);
    fclose(fid);
    