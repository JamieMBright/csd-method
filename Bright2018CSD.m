function csd = Bright2018CSD(ghi,ghics,dni,dnics,dif,difcs,zen,plot_figure)

%% thresholds
M = 5; %window size
r_ghi_max = 5;
r_ghi_min = 0.6;
r_dni_max = 10;
r_dni_min = 2;
r_dif_max = 40;
r_dif_min = 5;
p_ghi_max = 0.4;
p_ghi_min = 0.05;
p_dni_max = 0.5;
p_dni_min = 0.05;
p_dif_max = 0.1;

%increase the ghi and ghics in size to facilitate m timeteps ahead
ghim=[ghi;NaN(M,1)];
ghicsm=[ghics;NaN(M,1)];
dnim=[dni;NaN(M,1)];
dnicsm=[dnics;NaN(M,1)];
difm=[dif;NaN(M,1)];
difcsm=[difcs;NaN(M,1)];

% % rate change limits
h = 90-zen;
x = (0:90);
y_ghi = linspace(r_ghi_min,r_ghi_max,length(x));
y_dni = linspace(r_dni_min,r_dni_max,length(x));
y_dif = linspace(r_dif_min,r_dif_max,length(x));
idxs_z = knnsearch(x',h);
R_ghi = [y_ghi(idxs_z),NaN(1,M)]';
R_dni = [y_dni(idxs_z),NaN(1,M)]';
R_dif = [y_dif(idxs_z),NaN(1,M)]';

% abs fractional difference to clear-sky curve
y_ghi = linspace(p_ghi_max,p_ghi_min,length(x));
y_dni = linspace(p_dni_max,p_dni_min,length(x));
P_ghi = [y_ghi(idxs_z),NaN(1,M)]';
P_dni = [y_dni(idxs_z),NaN(1,M)]';
% real fractional diff upper limit
P_dif = p_dif_max;

% pre allocate the rate change variable
rate_ghi = NaN(length(ghim),M);
rate_dni = NaN(length(ghim),M);
rate_dif = NaN(length(ghim),M);
% for each differential...
% rate change difference
for m=1:M
    % get the indexing
    idx1=(m:length(ghim)-(M-m+1));
    idx2 =(m+1:length(ghim)-(M-m));
    % calculate the rate change.
    rate_ghi(1:length(idx1),m) = (ghim(idx1)-ghim(idx2))./(ghicsm(idx1)-ghicsm(idx2));
    rate_dni(1:length(idx1),m) = (dnim(idx1)-dnim(idx2))./(dnicsm(idx1)-dnicsm(idx2));
    rate_dif(1:length(idx1),m) = (difm(idx1)-difm(idx2))./(difcsm(idx1)-difcsm(idx2));
end
% take a mean across all m time steps.
rate_ghi = nanmean(rate_ghi,2);
rate_dni = nanmean(rate_dni,2);
rate_dif = nanmean(rate_dif,2);


rate_ghi(isinf(rate_ghi) | isnan(rate_ghi))=1;
rate_dni(isinf(rate_dni) | isnan(rate_dni))=1;
rate_dif(isinf(rate_dif) | isnan(rate_dif))=1;

% perentage difference
pdiff_ghi = abs(ghim-ghicsm)./ghicsm;
pdiff_dni = abs(dnim-dnicsm)./dnicsm;
pdiff_dif = abs(difm-difcsm)./difcsm;

% preallocate
csd_ghi = ones(size(ghim));
csd_dni = ones(size(ghim));
csd_dif = ones(size(ghim));

% criteria
csd_ghi((rate_ghi<1+R_ghi & rate_ghi>1-R_ghi) & pdiff_ghi<P_ghi) = 0;
csd_dni((rate_dni<1+R_dni & rate_dni>1-R_dni) & pdiff_dni<P_dni) = 0;
csd_dif((rate_dif<1+R_dif & rate_dif>1-R_dif) & pdiff_dif<P_dif) = 0;

csd_inverse=zeros(size(ghi));
csd_inverse(csd_ghi==0 & csd_dni==0 & csd_dif==0)=1;

% final duration criteria
duration_threshold = 15; %hour
durations = sum(hankel(csd_inverse,[csd_inverse(end),NaN(1,duration_threshold-1)]),2);

% final test if the duration of the clear-sky period has been at least 15
% mins
csd=ones(size(ghi));
csd(durations==duration_threshold)=0;


% figure
if exist('plot_figure','var')
    t1=1;
    t2= 6000;
    f=figure('name','Long Ackerman 2006 CSD example','color','w');
    f.Position = [250,250,1200,500];
    hold on
    plot(ghim)
    plot(ghicsm,'k:')
    CSD = ghim;
    CSD(csd_ghi==1)=NaN;
    plot(CSD,'linewidth',2,'color','k')
    CSD = ghim;
    CSD(csd==1)=NaN;
    plot(CSD,'linewidth',2,'color','g')
    plot(dnim)
    plot(dnicsm,'k:')
    CSD = dnim;
    CSD(csd_dni==1)=NaN;
    plot(CSD,'linewidth',2,'color','k')
    CSD = dnim;
    CSD(csd==1)=NaN;
    plot(CSD,'linewidth',2,'color','g')
    plot(difm)
    plot(difcsm,'k:')
    CSD = difm;
    CSD(csd_dif==1)=NaN;
    plot(CSD,'linewidth',2,'color','k')
    CSD = difm;
    CSD(csd==1)=NaN;
    plot(CSD,'linewidth',2,'color','g')
    
    hold off
    legend('GHI','Clear detected')
    ylabel(gca,'Irradiance [Wm^{-2}]')
    set(gca,'xlim',[t1,t2])
    set(gca,'xtick',[]);
    xlabel('Time')
end

end
