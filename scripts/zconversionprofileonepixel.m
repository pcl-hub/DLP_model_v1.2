function conversionprofile = zconversionprofileonepixel(formulation,printer,numberofLayers,...
                                layerThicknesses,exposureTimes,grayscaleValues,voidStructure,...
                                depthPenetration,criticalEnergy,criticalConversion,oxygenInhibition,...
                                ultimateConversion,rateConstant,reactionOrder,doseRateFactor,...
                                zinterval,stepplot,finalplot)
%% this function calculates z conversion profile for one pixel
                         
%% define all the parameters
nl = numberofLayers; % number of layers

dz = layerThicknesses; % layer thicknesses
dt = exposureTimes; % exposure times

gsv = grayscaleValues; % grayscale values;

vstr = voidStructure;

Dp = depthPenetration; % resin depth of penetration
Ec = criticalEnergy; % resin critical energy

ac = criticalConversion; % resin critical conversion

Eoi = oxygenInhibition;
woi = 1; % woi = oxygenInhibitFactor;

au = ultimateConversion;
k = rateConstant;
n = reactionOrder;
w = doseRateFactor;

% make t, dz, and gsv n-by-1 vectors
if ~(length(dt)==nl); dt = dt(1)*ones(nl,1);end
if ~(length(dz)==nl); dz = dz(1)*ones(nl,1);end
if ~(length(gsv)==nl); gsv = gsv(1)*ones(nl,1);end
if ~(length(vstr)==nl); vstr = vstr(1)*ones(nl,1);end

% default formulation and printer if not given
if isempty(formulation)
    formulation = struct('material','DA-2 Resin',...
                         'photoinitiator','PPO',...
                         'piConcentration',0.007,...
                         'photoabsorber',[],...
                         'paConcentration',[]);
end
if isempty(printer)
    printername = 'Anycubic';
    unitpd = 0.02; darkpd = 0.03;
    printer = struct('name',printername,...
                     'unitPowerDensity',unitpd,...
                     'darkPowerDensity',darkpd);
else
    unitpd = printer.unitPowerDensity;
    darkpd = printer.darkPowerDensity; 
    % power density I = unitpd*gsv+darkpd
end

%% calculate conversion profile
z = (sum(dz(:)):-zinterval:zinterval)'/1000;

effectiveDoseProfile = zeros(size(z));
conversionProfile = zeros(size(z));

ptInds = zeros(nl,2);
layerStatuses = zeros(nl,nl);
conversionProfiles = zeros(length(z),nl);

layers(1:nl,1) = struct('z',[],'effectiveDose',[],'conversion',[]);

tic
for i = nl:-1:1 % this loop goes through the layers
    %% initialize the points in layer i
    np = dz(i)/zinterval;
    d = (0:-1:(1-np))'*zinterval; % number of points in layer i set by dz(i)/step
    dp = zeros(np,1); % dose profile, initial values: 0
    cp = zeros(np,1); % conversion profile, initial values: 0
    rEoi = Eoi*ones(np,1); % residual oxygen
    resinRefresh = true;
    
    ptInds(i,1) = sum(dz(1:i-1))/zinterval+1;
    ptInds(i,2) = sum(dz(1:i))/zinterval;
    
    
    %% the loop below calculates the conversions for layer i when printing layer j
    for j = i:nl
        %% increment in depth and recalculates the light intensity
        d = d+dz(j); % increase the depth of the points in layer i by the thickness of layer j, dz(j)
        I0 = unitpd*gsv(j)+darkpd; % the printer light intensity set by the grayscale values for layer j
        I = I0*exp(-d/Dp); % the light intensity that reaches layer i,
                           % whose value is controlled by the depth of penetration
        
        %% accumulated dose profile
        dp = dp+I.^w*dt(j);
        
        %% equivalent time and effective time
        teq = (au^(1-n)-(au-cp).^(1-n))./((1-n)*k.*I.^w); % calculate the equivalent start time 
                                                          % from the start conversion of layer i,
                                                          % before the printing of layer j
                                                          % based on the cure kinetics (I) for layer i
        
        %% oxygen inhibition
        toi = rEoi./(I.^woi); % oxygen inhibition time length; woi is set to 1, see line 25
        
        %% effective time for cure kinetics
        teff = teq+dt(j)-toi;
        oi = teff<0;
        teff(oi) = 0;
        
        %% residual Eoi
        rEoi(~oi) = 0;
        rEoi(oi) = rEoi(oi)-I(oi).^woi*dt(j);

        %% calculate conversions from the kinetic model based on effective time
        cp = au-(au^(1-n)-(1-n)*k*(I.^w).*teff).^(1/(1-n));
        %%
        [cp,dp,rEoi,layerStatuses,resinRefresh] = structureproperty(cp,dp,rEoi,layerStatuses,resinRefresh,vstr,Eoi,ac,i,j);
        
        conversionProfiles(ptInds(i,1):ptInds(i,2),j) = cp;
    end
    layers(i).z = d/1000;
    layers(i).effectiveDose = dp;
    layers(i).conversion = cp;
    
    effectiveDoseProfile(ptInds(i,1):ptInds(i,2)) = dp;
    conversionProfile(ptInds(i,1):ptInds(i,2)) = cp;
end
toc

%% assign results
conversionprofile = struct('formulation',formulation,'printer',printer,...
        'resinProperties',struct('depthPenetration',Dp,'criticalEnergy',...
        Ec,'criticalConversion',ac,'oxygenInhibition',Eoi),'cureKinetics',...
        struct('ultimateConversion',au,'rateConstant',k,'reactionOrder',n,...
        'doseRateFactor',w),'printingParams', struct('numberofLayers',nl,...
        'layerThicknesses',dz,'exposureTimes',dt,'grayscaleValues',gsv,...
        'voidStructure',vstr),'layerStatuses',layerStatuses,'layers',layers,...
        'prediction',struct('z',z,'effectiveDoseProfile',effectiveDoseProfile,...
        'conversionProfile',conversionProfile),'units',struct('depthPenetration',...
        '\mum','criticalEnergy','mJ/cm^{2}','criticalConversion',[],...
        'oxygenInhibition','mJ/cm^{2}','ultimateConversion',[],'rateConstant',[],...
        'reactionOrder',[],'doseRateFactor',[],'numberLayers',[],...
        'layerThicknesses','\mum','exposureTimes','seconds',...
        'grayscaleValues',[],'z','mm','doseProfile',[],'conversionProfile',[]));


%% plot the conversion profile

if stepplot

afigure; hold on; grid off

set(gcf,'Position',[240 120 360 500]);

set(gca,'fontsize',12);

xlabel('conversion, \alpha');
ylabel('depth, {\itz} (mm)');

xlim([0 0.8]);
ylim([0 max(z)]);

tickspacings = [0.05 0.1 0.2 0.4 1 2 4];
ytkmaxs = [0 0.2 0.5 1 3 6 10 inf];
yi = (max(z)-ytkmaxs(1:end-1))>0 & ...
    ~((max(z)-ytkmaxs(2:end))>0);
tickspacing = tickspacings(yi);
if tickspacing<max(dz/1000)
    tickspacing = min(dz/1000);
end
yticks(0:tickspacing:max(z));

ax = gca;

minortickvalues = zeros(1,nl);
for i = 1:nl
    minortickvalues(i) = sum(dz(i:nl))/1000;
end
minortickvalues = flip(minortickvalues);

ax.YAxis.MinorTickValues = minortickvalues;
set(ax,'GridAlpha',0.8,'GridLineStyle',':',...
    'MinorGridAlpha',0.9, 'YMinorGrid','on');


mc = distinguishable_colors(nl+6);
mc = mc(7:end,:);

zs = z*ones(1,nl);

for j=1:nl
    zs(:,j) = zs(:,j)-sum(dz(j+1:nl))/1000;
    if j == 1
        title(sprintf(['Printer initialized.\nStarting to print Layer ' ...
              num2str(j) ' ...\nPress any key to continue...']),...
             'fontsize',10); pause;
         gcap = get(gca,'position');
         origin0 = gcap([1 2]);
         fullbox = gcap([3 4]);
    else
        title(sprintf(['Layer ' num2str(j-1) ' done.' ...
             '\nPrinting Layer ' num2str(j) ' ...' ...
             '\nPress any key to continue...']),...
             'fontsize',10); pause; 
    end


for i = 1:j
    zp = zs(ptInds(i,1):ptInds(i,2),j);
    for m = i:j
    cp = conversionProfiles(ptInds(i,1):ptInds(i,2),m);
    nz = cp>0;
    if m==1; cla; end
    if m==j
        plot(cp(nz),zp(nz),'color',mc(m,:),...
            'Linewidth',2,'Linestyle','-');
        plot(cp(~nz),zp(~nz),'color',mc(m,:),...
            'Linewidth',2,'Linestyle','-');
    else
        plot(cp(nz),zp(nz),'color',mc(m,:),...
            'Linewidth',1.5,'Linestyle',':');
        plot(cp(~nz),zp(~nz),'color',mc(m,:),...
            'Linewidth',1.5,'Linestyle',':');
    end
    end
end

plot(ac*ones(1,2),max(zs(:,j))*[-0.1 1.1],...
    'color','b','linewidth',1.5,'linestyle','-');
% scatter(ac*ones(size(z)),z,'sizedata',3,'marker','o',...
%     'markerfacecolor','b','markeredgecolor','none');

set(gca,'fontsize',12);

xlabel('conversion, \alpha');
ylabel('depth, {\itz} (mm)');
xlim([0 0.8]);
ylim([0 max(zs(:,j))]);

boxsize = [fullbox(1) fullbox(2)*max(zs(:,j))/max(z)];
cornerp = [origin0(1) origin0(2)+fullbox(2)-boxsize(2)];

set(gca,'position',[cornerp boxsize]);

tickspacings = [0.05 0.1 0.2 0.4 1 2 4];
ytkmaxs = [0 0.2 0.5 1 3 6 10 inf];
yi = (max(zs(:,j))-ytkmaxs(1:end-1))>0 & ...
    ~((max(zs(:,j))-ytkmaxs(2:end))>0);
tickspacing = tickspacings(yi);
if tickspacing<max(dz(1:j)/1000)
    tickspacing = max(dz(1:j)/1000);
end
yticks(0:tickspacing:max(zs(:,j)));

ax = gca;

minortickvalues = zeros(1,j);
for i = 1:j
    minortickvalues(i) = sum(dz(i:j))/1000;
end
minortickvalues = flip(minortickvalues);

ax.YAxis.MinorTickValues = minortickvalues;
set(ax,'GridAlpha',0.8,'GridLineStyle',':',...
    'MinorGridAlpha',0.9, 'YMinorGrid','on');

end

title(sprintf('\nPrinting Completed\n'),'fontsize',10);
end

if finalplot
afigure; 
hold on
set(gca,'fontsize',12);
set(gcf,'Position',[720 120 360 480]);

xlabel('conversion, \alpha');
ylabel('depth, {\itz} (mm)');

xlim([0 0.8]);
ylim([0 max(z)]);

ax = gca;
gcap = get(ax,'position');
gcap = gcap-[0 0 0.01 0.01];
set(ax,'position',gcap);

tickspacings = [0.05 0.1 0.2 0.4 1 2 4];
ytkmaxs = [0 0.2 0.5 1 3 6 10 inf];
yi = (max(z)-ytkmaxs(1:end-1))>0 & ...
    ~((max(z)-ytkmaxs(2:end))>0);
tickspacing = tickspacings(yi);
if tickspacing<max(dz/1000)
    tickspacing = max(dz/1000);
end
yticks(0:tickspacing:max(z));

minortickvalues = zeros(1,nl);
for i = 1:nl
    minortickvalues(i) = sum(dz(i:nl))/1000;
end
minortickvalues = flip(minortickvalues);
ax.YAxis.MinorTickValues = minortickvalues;
set(ax,'GridAlpha',0.8,'GridLineStyle',':',...
    'MinorGridAlpha',0.9, 'YMinorGrid','on');

h1 = plot(ac*ones(1,2),max(z)*[-0.1 1.1],...
    'color','b','linewidth',1.5,'linestyle','-');

% scatter(ac*ones(size(z)),z,'sizedata',3,...
%     'marker','o','markerfacecolor','b',...
%     'markeredgecolor','none'); % critical conversin line

for i = 1:nl
    zp = z(ptInds(i,1):ptInds(i,2));
    cp = conversionProfile(ptInds(i,1):ptInds(i,2));
    nz = cp>0;
   
    plot(cp(nz),zp(nz),'-r');
    plot(cp(~nz),zp(~nz),'-r');
end


% legend(h1,'{\alpha}_{c}',...
%     'location','northeast','fontsize',11);
end


%% Do not edit the script below - JT %%

function [cp,dp,rEoi,layerStatuses,resinRefresh] = structureproperty...
             (cp,dp,rEoi,layerStatuses,resinRefresh,vstr,Eoi,ac,i,j)
ii = cp<ac; % indices for points below critical conversion

switch j
    case i
        switch vstr(i)
            case 0
                cp(ii) = 0;
                dp(ii) = 0;
                rEoi(ii) = Eoi;
                if all(ii)
                    resinRefresh = true;
                    layerStatuses(i,j) = 0;
                elseif all(~ii)
                    resinRefresh = false;
                    layerStatuses(i,j) = 1;
                else
                    resinRefresh = true;
                    layerStatuses(i,j) = 0.5;
                end

            case 1
                if all(ii)
                    cp(ii) = 0;
                    dp(ii) = 0;
                    rEoi(ii) = Eoi;
                    layerStatuses(i,j) = 0; 
                    resinRefresh  = true;
                elseif all(~ii)
                    layerStatuses(i,j) = 1; 
                    resinRefresh  = false;
                else
                    layerStatuses(i,j) = 0.5; 
                    resinRefresh  = false;
                end
        end
    case i+1
        if resinRefresh

        switch vstr(i)
            case 0
                cp(ii) = 0;
                dp(ii) = 0;
                rEoi(ii) = Eoi;
                if all(ii)
                    resinRefresh = true;
                    layerStatuses(i,j) = 0;
                elseif all(~ii)
                    resinRefresh = false; 
                    layerStatuses(i,j) = 1;
                else
                    resinRefresh = true;
                    layerStatuses(i,j) = 0.5;
                end
            case 1
                switch vstr(j)
                    case 0
                        switch layerStatuses(j,j)
                            case {0, 0.5}
                                cp(ii) = 0;
                                dp(ii) = 0;
                                rEoi(ii) = Eoi;
                                resinRefresh = true;
                                layerStatuses(i,j) = 0;
                            case 1
                                resinRefresh = false;
                                if all(ii)
                                    layerStatuses(i,j) = 0;
                                elseif all(~ii)
                                    layerStatuses(i,j) = 1;
                                else
                                    layerStatuses(i,j) = 0.5;
                                end
                        end
                    case 1
                        switch layerStatuses(j,j)
                            case 0
                                cp(ii) = 0;
                                dp(ii) = 0;
                                rEoi(ii) = Eoi;
                                resinRefresh = true;
                                layerStatuses(i,j) = 0;
                            case {0.5, 1}
                                resinRefresh = false;
                                if all(ii)
                                    layerStatuses(i,j) = 0;
                                elseif all(~ii)
                                    layerStatuses(i,j) = 1;
                                else
                                    layerStatuses(i,j) = 0.5;
                                end
                        end
                end
        end
        else
            layerStatuses(i,j) = layerStatuses(i,j-1);
        end

    otherwise
        if resinRefresh

        switch vstr(i)
            case 0
                cp(ii) = 0;
                dp(ii) = 0;
                rEoi(ii) = Eoi;
                if all(ii)
                    resinRefresh = true;
                    layerStatuses(i,j) = 0;
                elseif all(~ii)
                    resinRefresh = false; 
                    layerStatuses(i,j) = 1;
                else
                    resinRefresh = true;
                    layerStatuses(i,j) = 0.5;
                end
            case 1
                switch vstr(i+1)
                    case 0
                        cp(ii) = 0;
                        dp(ii) = 0;
                        rEoi(ii) = Eoi;
                        resinRefresh = true;
                        layerStatuses(i,j) = 0;

                    case 1
                        switch layerStatuses(j,j)
                            case 0
                                cp(ii) = 0;
                                dp(ii) = 0;
                                rEoi(ii) = Eoi;
                                resinRefresh = true;
                                layerStatuses(i,j) = 0;
                            case {0.5, 1}
                                resinRefresh = false;
                                if all(ii)
                                    layerStatuses(i,j) = 0;
                                elseif all(~ii)
                                    layerStatuses(i,j) = 1;
                                else
                                    layerStatuses(i,j) = 0.5;
                                end
                        end
                end
        end
        else
            layerStatuses(i,j) = layerStatuses(i,j-1);
        end
end