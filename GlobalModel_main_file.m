clc
parametermapplot=true;

figures = get(groot, 'Children');
if parametermapplot
    if ~isempty(figures)
        fig20 = find([figures.Number] == 20, 1);
        if isempty(fig20)
            clear all
            parametermapplot=true;
            numSurfaces=0;
            surfhandles={};
            handles2d={};
            surflabels={};
            figure(20)
        else
            clearvars -except R0_PS numSurfaces surfhandles surflabels parametermapplot handles2d
        end
    else %if empty(figures)
        clear all
        parametermapplot=true;
        numSurfaces=0;
        surfhandles={};
        surflabels={};
        figure(20)
    end
else 
    clear all
    parametermapplot=false;
end
plot2dslices=true;
%legendstr2d={'Original cell parameters (OCP)','OCP + max. Li. concentration','OCP + diffusivity','OCP + porosity','OCP + thickness','OCP + radius','Target cell'};
legendstr2d={'Original cell parameters (OCP)','Target cell parameters'};

NoSOC=false; %set false to include SOC dependence, true to reduce model to SOC independence
Nocurrent=false;
noTemp=true;
static=false; %Neither SOC, current nor Temperature
pybamm=true; %use the orginal parameters or stuff from pybamm

smooth=false; %smooth parameter maps for visual representation
pybamm_foldername='Marquis';
pybamm_name='Marquis';
%pybamm_foldername='Chen';
%pybamm_name='Chen';
sol_data_name='sol_data_Chen_Verification.mat';
OCV_filename='OCV_SOC_Chen.mat';

%Chen_Marquis_MaxConcentrationPosEl Chen_Marquis_PosElDiffusivity
%Chen_Marquis_PosElPorosity Chen_Marquis_PosElThickness Chen_Marquis_PosParticleRadius
titletext='';
k_PS=1;
FaceAlpha_value=0.5; %transparancy value

facecolorvalue='r';
V_ECM_color='b';
legendstrV_ECM={'Chen ECM','Chen DFN'};

pybamm_shift=false; %vertical shift for plotting multiple surfaces
pybamm_normalize=false;
start=1; %Start time 
endglobalmodel=1;%*5000; number of chunks
fs = 20;
run param/params_LCO
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
OneC=p.OneC;
filesofinterest=5;
if pybamm
    filesofinterest=1; %for later loop over files
    load(sprintf('./Results/HPPC pybamm %s/Parameters_pybamm_%s.mat', pybamm_foldername, pybamm_name))
end



mat_files = dir('results_*.mat'); %Merge data
all_SOC = [];
all_I = [];
all_Volt= [];
all_T=[];
for i = 1:length(mat_files)
    load(mat_files(i).name, 'SOC_save','Volt_save', 'I_save','T_save');
    if iscolumn(I_save)
        I_save=I_save';
    end
    if iscolumn(SOC_save)
        SOC_save=SOC_save';
    end
    if iscolumn(T_save)
        T_save=T_save';
    end
    if iscolumn(Volt_save)
        Volt_save=Volt_save';
    end
    all_SOC = [all_SOC SOC_save];
    all_I = [all_I I_save]; %ATT sometimes needs a '
    all_Volt = [all_Volt Volt_save];
    all_T = [all_T T_save];
end
t=0:length(all_I)-1;
if pybamm && ~parametermapplot
    
    load(sol_data_name, 'I','OneC', 'SoC','T','V');
    t = 1:max(T);
    all_I = interp1(T, I, t, 'linear');
    all_SOC = interp1(T, SoC, t, 'linear');
    all_Volt = interp1(T, V, t, 'linear');
end

if exist("endglobalmodel",'var')==0
    endglobalmodel=ceil((length(all_I)-start)/5000);
    endglobalmodelplot=length(t);
else
    endglobalmodelplot=endglobalmodel*5000+start-1;
end




TempRange = './Parameters_Temp/pybamm'; %ATT: noSOC: './Parameters_Temp_noSOC'
files1 = dir(fullfile(TempRange, 'Parameters1RC_*.mat')); % ATT: for noSOC: 'Parameters1RC_noSOC_*.mat'
if pybamm && parametermapplot
    TempRange = sprintf('./Results/HPPC pybamm %s', pybamm_foldername);   
    files1 = dir(fullfile(TempRange, sprintf("Parameters_pybamm_%s.mat", pybamm_name)));
end
if NoSOC
    TempRange = './Parameters_Temp_noSOC'; 
    files1 = dir(fullfile(TempRange, 'Parameters1RC_noSOC_*.mat'));
elseif Nocurrent
    TempRange = './Parameters_Temp_nocurrent'; 
    files1 = dir(fullfile(TempRange, 'Parameters1RC_nocurrent_*.mat'));
elseif static
    TempRange = './Parameters_static'; 
    files1 = dir(fullfile(TempRange, 'Parameters1RC_nocurrent_*.mat'));
% elseif noTemp
% 
end

OCV_TempRange= './OCV_SOC_Temp/pybamm'; 
%OCV_TempRange= './Results_CHEN/OCV_SOC_Temp'; %ATT

files2 = dir(fullfile(OCV_TempRange, OCV_filename));


V_ECM_Temp=cell(length(files1),1);

Temperatures=cell(length(files1),1);
Params_c=cell(length(files1),1);
Params_d=cell(length(files1),1);

for i = filesofinterest%:7%1:length(files1) %5:25°C  
    filePath1 = fullfile(TempRange, files1(i).name) 
    filePath2 = fullfile(OCV_TempRange, files2(1).name) %always use the
    %OCV_SOC curve for 25 degrees ATT: i=5
    %filePath2 = fullfile(OCV_TempRange, files2(3).name); %ATT always use the OCV_SOC curve for 25 degrees

    
    data1 = load(filePath1);
    data2 = load(filePath2);
    tempPattern = 'Parameters1RC_(.*?)deg'; %ATT: for noSOC: 'Parameters1RC_noSOC_(.*?)deg'
    if NoSOC
        tempPattern = 'Parameters1RC_noSOC_(.*?)deg'; 
    elseif Nocurrent
        tempPattern = 'Parameters1RC_nocurrent_(.*?)deg'; 
    end
    tokens = regexp(files1(i).name, tempPattern, 'tokens');
    if ~isempty(tokens)
        temperature = tokens{1}{1}; % Extracted temperature value as a string
    else
        temperature = 'Unknown'; % Default value if pattern not found
    end
    Temperatures{i}=num2str(temperature);
    if static
        temperature='25';
    end
    
    % run global model
    if ~static
        [V_ECM, Params_c{i},Params_d{i}] = GlobalModel_piecewise_function(data1,data2,all_I(start:end),all_SOC(start:end),t(start:end),NoSOC,Nocurrent,static,endglobalmodel);
    else
        V_ECM = GlobalModel_piecewise_function(data1,data2,all_I(start:end),all_SOC(start:end),t(start:end),NoSOC,Nocurrent,static,endglobalmodel);
    end
    %V_ECM = GlobalModel1RC_piecewise_function(I_RCfit,SOC_RCfit,t);
    V_ECM_Temp{i}=V_ECM;
    
    if strcmp(temperature, '25') && noTemp
        V_ECM_Temperature_interp=V_ECM;
        save('ECM_noTemp.mat', 't', 'V_ECM_Temperature_interp','all_I','all_T','all_SOC','start')
    end

    if pybamm && ~parametermapplot
        figure(8)
        errors = all_Volt(4900:5420) - V_ECM(4900:5420);
        rmse = sqrt(mean(errors.^2))*1000;
            
        title(['Chen (OCV Marquis), standard. RMSE ', num2str(rmse), ' mV'],'interpreter','latex')
        hold on
        plot(t(start:endglobalmodelplot)+1,V_ECM,'LineWidth',2,'Color',V_ECM_color)
        plot(t(start:endglobalmodelplot)+1,all_Volt,'LineWidth',2,'Color',V_ECM_color,'LineStyle','--')

        legend(legendstrV_ECM,'interpreter','latex','Location','southeast')
        %legend({'DFN Chen mod','ECM Chen mod', 'DFN Chen unmod + OCV Marquis','ECM Chen unmod + OCV Marquis','DFN Chen unmod','ECM Chen unmod','DFN Marquis','DFN Chen mod + OCV Marquis','ECM Chen mod + OCV Marquis'},'interpreter','latex','Location','southeast')
        ylabel('Voltage','FontSize' , fs)
        xlabel('time','FontSize', fs)
        figure(9)
        plot(all_Volt(4900:5420),'LineWidth',2);
        hold on
        plot(V_ECM(4900:5420),'LineWidth',2);   
        title(['RMSE ', num2str(rmse), ' mV'],'interpreter','latex')
        legend({'Marquis (DFN)','Marquis (ECM)'},'interpreter','latex','Location','southeast')

        %remove second plot
        % ax = findobj(gcf, 'type', 'axes');
        % delete(ax(1)); % ax(1) refers to the most recently created subplot, which is 'Second Subplot'
        %

        % ax2=subplot(2,1,2);
        % plot(t(start:endglobalmodelplot),all_I(start:endglobalmodelplot),'LineWidth',2,'Color','b')
        % ylabel('Current','FontSize' , fs)
        % xlabel('time','FontSize', fs)
        % set(gca,'FontSize', fs)
        % linkaxes([ax1 ax2],'x')
        % Error{i}=all_Volt(start:endglobalmodelplot) - V_ECM;
        % RMSE{i} = sqrt(mean((all_Volt(start:endglobalmodelplot) - V_ECM).^2))*1000;
        % fprintf('Temperature: % 3.2f RMSE : %3.2f mV',str2double(Temperatures{i}) ,RMSE{i})
    end
    

end

if ~pybamm
    figure()
    hold on
    ylabel('Voltage','FontSize', fs)
    xlabel('time','FontSize', fs)
    legendStr = {};
    V_ECM_Temp_nonempty=find(~cellfun(@isempty,V_ECM_Temp));
    for V_ECM_Temp_idx=V_ECM_Temp_nonempty'
        plot(t(start:endglobalmodelplot), V_ECM_Temp{V_ECM_Temp_idx}, 'LineWidth', 2,'Color','b')
        legendStr{end + 1} = ['Temperature (' Temperatures{i} 'C)'];
        legend(legendStr, 'interpreter', 'latex', 'Location', 'southeast');
    end
    legendStr{end + 1} = ['DFN'];
    plot(t(start:endglobalmodelplot), all_Volt(start:endglobalmodelplot), 'LineWidth', 2,'Color','r')
    legend(legendStr, 'interpreter', 'latex', 'Location', 'southeast');
end

%3D-Plot: R0 over SOC, I


if smooth
    '3D Maps are smoothed'
end
azimuth=45;
elevation=45;
SOC_of_interest=0.1;
current_of_interest=50; 
R0_SOC=cell(size(Temperatures,1),1);
R0_I=cell(1,size(Temperatures,1));
i_start=1; %att: 3 usually 
for i=filesofinterest%:7%i_start:length(Params)

    ParamStruct_c=Params_c{i};
    ParamStruct_d=Params_d{i};
    current_c=ParamStruct_c.current;
    current_d=ParamStruct_d.current;
    SOC=ParamStruct_c.SOC; %SOC same for charge and discharge
    R0_c=ParamStruct_c.R0;
    R0_d=ParamStruct_d.R0;
    R0_cd=zeros([size(R0_d,1),size(R0_d,2)*2-1]); %0 exists twice
    current_cd=zeros([length(current_c)*2-1,1]);
    current_cd(1:length(current_d)-1)=current_d(1:end-1);
    current_cd(length(current_c):end)=current_c(1:end);
    current_cd=flip(current_cd);
    R0_cd(:,1:length(current_d)-1)=R0_d(:,1:end-1);
    R0_cd(:,length(current_c):end)=R0_c(:,1:end);
    
    %particle size start (uncomment below)

    %R0_PS=zeros(8,size(R0_cd,1));
    % if size(R0_cd,1) ~= size(R0_PS,2)
    %     lenDiff = size(R0_PS,2) - size(R0_cd,1);
    %     R0_cd15=R0_cd(:,15);
    %     lastValue = R0_cd15(end-1);
    %     closestHigherValues = repmat(lastValue, lenDiff, 1);
    %     R0_cd15 = [R0_cd15; closestHigherValues];
    %     R0_PS(k_PS,:)=R0_cd15;
    % else
    %     R0_PS(k_PS,:)=R0_cd(:,15); %15 = -1C
    % end

    % figure
    % PartSize=[25,50,75,100,125,150,175,200]';
    % standard_SOC=1:-0.025:0;
    % surf(PartSize,standard_SOC,R0_PS'*1000)
    % title('R0 at -1C over particle size and SOC')
    % xlabel('particle size in % of original')
    % ylabel('SOC')
    % zlabel('R0 [m$\Omega$]', 'Interpreter','latex')
    % ax=gca;
    % ax.FontSize = fs;
    % ax.XTick = linspace(min(PartSize), max(PartSize), 8); % Adjust 'numberOfTicks' as needed
    
    %particle size end

    %ATT not yet corrected
    tau_cell{i}=ParamStruct_c.tau; %ATT
    SOC_idx_of_interest=find(round(SOC,2)==SOC_of_interest);
    current_idx_of_interest=find(round(current_c,0)==current_of_interest);
    Temperature_now_cell{i-i_start+1}=str2num(Temperatures{i});
    R0_SOC{i}=R0_c(SOC_idx_of_interest,:);
    R0_I{i}=R0_c(:,current_idx_of_interest);

    % if i==1  
    %     current_standard=current;
    %     SOC_standard=SOC;
    % end
    % if i==3 %To reduce the larger R0 matrix at 25 degrees to the smaller size at the other Temp. Instead of interpolating the smaller matrices up.
    %     [~,current_idx_standard_at25]=ismember(round(current_standard),round(current));
    %     current_idx_standard_at25=nonzeros(current_idx_standard_at25);
    %     [~,SOC_idx_standard_at25]=ismember(round(SOC_standard,2),round(SOC,2));
    %     SOC_idx_standard_at25=nonzeros(SOC_idx_standard_at25);
    %     R0_SOC{i}=R0_SOC{i}(current_idx_standard_at25);
    %     R0_I{i}=R0_I{i}(SOC_idx_standard_at25);
    %     tau_cell{i}=tau_cell{i}(SOC_idx_standard_at25);
    % end
    windowSize=7.5;
    if smooth
        
        R0_cd = smoothdata(R0_cd, 1, 'gaussian', windowSize); 
        R0_cd = smoothdata(R0_cd, 2, 'gaussian', windowSize); 
    end
    if i==5
        R0_nonoise=R0_cd;
    elseif i==6
        R0_noise=R0_cd;
    elseif i==7
        R0_noise_prop=R0_cd;
    end

    maxR0=max(R0_cd(:));
    hold on
    figure(20)
    if pybamm_shift
        fig = gcf;    
        ax = get(fig, 'CurrentAxes');
        totalNumberOfPlots = numel(get(ax, 'Children'));
        pybamm_shift_num=0.2*totalNumberOfPlots-0.2;
        if pybamm_normalize
            surf(current_cd/OneC,SOC,1000*R0_cd/maxR0+pybamm_shift_num,'FaceAlpha',FaceAlpha_value);
       else
            surf(current_cd/OneC,SOC,1000*R0_cd+pybamm_shift_num,'FaceAlpha',FaceAlpha_value);
        end
    else
        if pybamm_normalize
            surf(current_cd/OneC,SOC,1000*R0_cd/maxR0,'FaceAlpha',FaceAlpha_value);
        else
            surfhandles{numSurfaces+1}=surf(current_cd/OneC,SOC,R0_cd*1000,'FaceColor',facecolorvalue,'FaceAlpha',FaceAlpha_value);      
            surflabels{numSurfaces+1} = sprintf('ECM %s',  pybamm_name);
            %set(gca, 'YDir', 'reverse'); % Reverse the direction of the y-axis

        end
    end
    ax = gca;
    allChildren = ax.Children;
    surfObjects = allChildren(strcmp('surface', get(allChildren, 'Type')));
    numSurfaces = numel(surfObjects);
    legend([surfhandles{:}], surflabels);
    %legend([h1, h2], {'First Surface', 'Second Surface'});
    %save parameter map
    %save('ParameterMap_Chen_Marquis_diff_conc_radius_mod_OCV_Marquis.mat', 'current_cd','OneC','SOC','R0_cd');
    
    
    title(titletext,'Fontsize', fs, 'interpreter', 'latex')
    %title(sprintf('Parameter map $R_0$(current, SOC): %s: scaled by total capacity', pybamm_name),'FontSize', fs, 'interpreter', 'latex');
    %title(sprintf('Parameter maps $R_0$(current, SOC): scaled by max($R_0$)', pybamm_name),'FontSize', fs, 'interpreter', 'latex');
    xlabel('C-rate', 'FontSize', fs, 'interpreter', 'latex');
    ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
    zlabel('$R_0$ (m$\Omega$)', 'FontSize', fs, 'interpreter', 'latex');
    if pybamm_normalize
        zlabel('$R_0$/max($R_0$)', 'FontSize', fs, 'interpreter', 'latex');
    end
    %zlim([0,60])
    % hLeg=legend(['1x','2x','4x','0.5x'], 'Interpreter','latex')
    % set(hLeg,'visible','off')
    hold on
    view(-azimuth-90, elevation);
    %axis off
    grid on
    ax=gca;
    ax.FontSize = fs;
    ax.XTick = round(linspace(min(current_cd/OneC), max(current_cd/OneC), 6),1); 
    ax.YTick = round(linspace(min(SOC), max(SOC), 6),1); 
    figure20=figure(20);
    set(20, 'units', 'centimeters', 'pos', [0 0 35 20])
    %exportgraphics(tiles,'ECM visual2.pdf')
    %exportgraphics(figure20,'Chen_default_noMarquisOCV.pdf')
    exportgraphics(figure20,'PM_Chen_Marquis_pre.pdf')

    if plot2dslices
        ax=figure(21);
        hold on
        h=plot(SOC,R0_cd(:,15)*1000,'LineWidth', 2, 'Color','r','LineStyle','-.');      
        xlabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
        ylabel('$R_0$ (m$\Omega$)', 'FontSize', fs, 'interpreter', 'latex');
        ax = gca;
        ax.FontSize = fs;  
        legend(legendstr2d,'Interpreter','latex');
        % ax = gca;  % Get the current axes
        % lastPlot = ax.Children(end);  % Get the handle of the last added plot
        % delete(lastPlot);  % Delete the last plot
        ylim([20,40])
        ylim([0,180])
        ax=gca;
        ax.FontSize = fs;
        set(21, 'units', 'centimeters', 'pos', [0 0 35 15])
        %set(21, 'units', 'centimeters', 'pos', [0 0 35 20])
        textbox = annotation("textbox",'interpreter','latex','FontSize',fs);
        textbox.String='\textbf{Highlighted section}';
        box on
        xticks([0:0.1:1])
        %yticks([20:5:40])
        yticks([0:50:150])
        %exportgraphics(ax,'EC Influence_V3.pdf')
        exportgraphics(ax,'EC Influence Marquis_V4.pdf')
        
    end
end

%plot tau

figure(22)
hold on
plot(SOC(1:end),tau1d(1:end),'LineWidth', 2,'Color','r','LineStyle','-.')
xlabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
ylabel('$\tau_1$ (s)', 'FontSize', fs, 'interpreter', 'latex');
ylim([0,30])
ax=gca;
ax.FontSize = fs;
legend({'Original cell - NMC', 'Target cell - LCO'},'Interpreter','latex','Location','southeast')
set(22, 'units', 'centimeters', 'pos', [0 0 35 10])
box on
%exportgraphics(ax,'tau_Chen.pdf')
exportgraphics(ax,'tau chen marquis.pdf')

%% plot parameter map difference
figure
h1=surf(current_cd/OneC,SOC,abs(R0_nonoise-R0_noise_5)*1000,'FaceAlpha',FaceAlpha_value, 'FaceColor','b');
hold on
h2=surf(current_cd/OneC,SOC,abs(R0_nonoise-R0_noise_10)*1000,'FaceAlpha',FaceAlpha_value, 'FaceColor','r');
h3=surf(current_cd/OneC,SOC,abs(R0_nonoise)*1000,'FaceAlpha',FaceAlpha_value, 'FaceColor','g');

%title('Absolute Parameter Map Difference','FontSize', fs, 'interpreter', 'latex');
xlabel('C-Rate', 'FontSize', fs, 'interpreter', 'latex');
ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
zlabel('$\Delta R_0$ [$m\Omega$]', 'FontSize', fs, 'interpreter', 'latex');
view(azimuth, elevation);
ax=gca;
ax.FontSize = fs;
ax.XTick = round(linspace(min(current_cd/OneC), max(current_cd/OneC), 6),1); 
ax.YTick = round(linspace(min(SOC), max(SOC), 6),1); 
%zlim([0,round(max(R0*1000,[],'all'),-1)]);

legend([h1, h2, h3], {'Noise: 5\% standard deviation of voltage', 'Noise: 10\% standard deviation of voltage','Default parameter map'},'FontSize',fs,'Interpreter','latex');



%%
R0_SOC{1}(16)=R0_SOC{1}(15); %ATT a hack
R0_SOC=cell2mat(R0_SOC); %R0 as a function of T, I
R0_I=cell2mat(R0_I); %R0 as a function of T, SOC
Temperature_now=cell2mat(Temperature_now_cell);
tau_matrix=cell2mat(tau_cell);


figure
%surf(current_standard,Temperature_now,R0_SOC*1000)
surf(current,Temperature_now,R0_SOC*1000)
xlabel('Current [A]', 'FontSize', fs, 'interpreter', 'latex');
ylabel('Temperature [$^{\circ}$C]', 'FontSize', fs, 'interpreter', 'latex');
zlabel('$R_0$ [$m\Omega$]', 'FontSize', fs, 'interpreter', 'latex');
title(['$R_0$(Current, Temperature) at SOC: ' num2str(SOC_of_interest*100) '$\%$'],'FontSize', 0.8*fs, 'interpreter', 'latex');
view(azimuth, elevation);
figure
%surf(Temperature_now,SOC_standard,R0_I*1000)
surf(Temperature_now,SOC,R0_I*1000)
xlabel('Temperature [$^{\circ}$C]', 'FontSize', fs, 'interpreter', 'latex');
ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
zlabel('$R_0$ [$m\Omega$]', 'FontSize', fs, 'interpreter', 'latex');
title(['$R_0$(Temperature, SOC) at current: ' num2str(current_of_interest) '[A]'],'FontSize', 0.8*fs, 'interpreter', 'latex');
view(azimuth, elevation);
figure
%surf(Temperature_now,SOC_standard,tau_matrix)
surf(Temperature_now,SOC,tau_matrix)
xlabel('Temperature [$^{\circ}$C]', 'FontSize', fs, 'interpreter', 'latex');
ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
zlabel('$\tau$ [$s$]', 'FontSize', fs, 'interpreter', 'latex');
title(['$\tau$(Temperature, SOC)'],'FontSize', 0.8*fs, 'interpreter', 'latex');
view(azimuth, elevation);
%%
%clc
%close all
figure
hold on
grid on
legend_str={};
for i=[2,1,3:length(Temperatures)]
    plot(SOC(2:end-1), R0_I(2:end-1,i)*1000, 'LineWidth', 2);
    legend_str{end+1}=[strcat(string(Temperatures{i}), '$ ^\circ$ C')];
    %legendStr{end + 1} = ['Temperature (' Temperatures{i} 'C)'];
    %legend(legendStr, 'interpreter', 'latex', 'Location', 'southeast');
end
xlim([0,1])
xlabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
ylabel('$R_0 (m\Omega)$', 'FontSize', fs, 'interpreter', 'latex');
legend(legend_str, 'Interpreter','latex')
title(['$R_0$ at current: ' num2str(current_of_interest) '[A]'], 'Interpreter','latex')
%% Interpolate over Temperature
close all
T=all_T(start:endglobalmodelplot);
t_used=t(start:endglobalmodelplot)';
V_ECM_Temperature=cell2mat(V_ECM_Temp);
NumericTemperatures = cellfun(@str2double, Temperatures);
if NumericTemperatures(1)>NumericTemperatures(2)
    old=NumericTemperatures(1);
    NumericTemperatures(1)=NumericTemperatures(2);
    NumericTemperatures(2)=old;
end
NumericTemperatures = NumericTemperatures(~isnan(NumericTemperatures))+273.15;
%weights_temp=interp1(NumericTemperatures,flip(eye(5)),T,'linear','extrap');
if ~static
    V_ECM_Temperature_interp=interp2(NumericTemperatures,t_used,V_ECM_Temperature',T,t_used');
else
    V_ECM_Temperature_interp=V_ECM_Temperature(1,:); %not interpolated obviously.
    "Att: no Temp interpolation. constant parameter model"
end
%temp=sum(V_ECM_Temperature'.*weights_temp,2);
%%
if noTemp
    V_ECM_Temperature_interp=V_ECM_Temp{1};
    "Att: no Temp interpolation. T=const=25 deg C"
end

% Main Temperature adjusted plot
include_voltage_plot = true;
include_current_plot = true;
include_temperature_plot = false;
include_SOC_plot = true;
figure()
tiledlayout(3,2); 

if exist("endglobalmodel",'var')==1
    t=t(1:endglobalmodelplot);
    all_I=all_I(1:endglobalmodelplot);
    %all_T=all_T(1:endglobalmodelplot);
    all_Volt=all_Volt(1:endglobalmodelplot);
    all_SOC=all_SOC(1:endglobalmodelplot);

end
save('ECM_NoTemp.mat', 't', 'all_Volt', 'V_ECM_Temperature_interp','all_I','all_T','all_SOC','start')
if include_voltage_plot
    nexttile(1, [2 2]);
    plot(t(start:end), all_Volt(start:end), 'LineWidth', 2, 'Color', 'r');
    hold on
    %plot(t(start+1:end), V_ECM_Temp{5}(1:end-1), 'LineWidth', 2,'Color','b');
    plot(t(start+1:end), V_ECM_Temperature_interp(1:end-1), 'LineWidth', 2,'Color','b');
    ylabel('Voltage [V]', 'FontSize', fs, 'interpreter', 'latex');
    xlabel('Time [s]', 'FontSize', fs, 'interpreter', 'latex');
    legend('DFN','ECM', 'interpreter', 'latex', 'Location', 'southeast');
    %title('Temperature adjusted ECM ','FontSize', fs/2, 'interpreter', 'latex')
end

if include_current_plot
    nexttile(5);
    plot(t(start:end), all_I(start:end)/OneC, 'LineWidth', 2, 'Color', 'b');
    ylabel('C-Rate', 'FontSize', fs, 'interpreter', 'latex');
    xlabel('Time [s]', 'FontSize', fs, 'interpreter', 'latex');
end

% if include_temperature_plot
%     nexttile(8);
%     plot(t(start:end), all_T(start:end)-273.15, 'LineWidth', 2, 'Color', 'b');
%     ylabel('Temperature [$^{\circ}$C] ', 'FontSize', fs, 'interpreter', 'latex');
%     xlabel('Time [s]', 'FontSize', fs, 'interpreter', 'latex');
% end

if include_SOC_plot
    nexttile(6);
    plot(t(start:end), all_SOC(start:end), 'LineWidth', 2, 'Color', 'b');
    ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
    xlabel('Time [s]', 'FontSize', fs, 'interpreter', 'latex');
end


ax = findobj(gcf,'Type','Axes');
linkaxes(ax, 'x');

%%


MSE = mean((all_Volt(start:end) - V_ECM_Temperature_interp).^2);
RMSE{i+1} = sqrt(MSE)*1000; %Error in millivolt
Error{i+1} = all_Volt(start:end)-V_ECM_Temperature_interp;
figure()

centered_error=Error{i+1}-mean(Error{i+1});


plot(t(start:end),centered_error*1000,'LineWidth',2,'Color','b')
ylabel('Error [mV]','FontSize' , fs, 'interpreter', 'latex')
xlabel('Time [s]','FontSize', fs,'interpreter', 'latex')
title('DFN-ECM Error','interpreter', 'latex')
figure
histogram(centered_error)
histfit(centered_error)
skew=skewness(centered_error);
kurt = kurtosis(centered_error);
title('Histogram fit','interpreter', 'latex')
xlim([-0.25,0.25])
ylabel('Counts','FontSize' , fs, 'interpreter', 'latex')
xlabel('Error (mV)','FontSize', fs,'interpreter', 'latex')
%%
all_SOC_s=all_SOC(start:end);
all_I_s=all_I(start:end);
all_T_s=all_T(start:end);
sz=3;

include_SOC_plot2 = true;
include_current_plot2 = true;
include_temperature_plot2 = true;

figure()
tiledlayout('vertical'); 



if include_SOC_plot2
    nexttile;
    scatter(all_SOC_s,centered_error,sz,'filled')
    p = polyfit(all_SOC_s, centered_error, 1);
    fitted_values = polyval(p, all_SOC_s);
    hold on;
    grid on;
    %lfit=plot(all_SOC_s, fitted_values, 'r-', 'LineWidth', 2); 
    %legend(lfit,'Linear Fit', 'interpreter', 'latex', 'Location', 'southeast');
    xlabel('SOC','FontSize' ,fs, 'interpreter', 'latex');
    ylabel('Error (mV)','FontSize' ,fs, 'interpreter', 'latex');
    title('Voltage Error over SOC','FontSize' ,fs, 'interpreter', 'latex');
    xlim([min(all_SOC_s) max(all_SOC_s)]);
    %corrcoef(all_SOC,centered_error)
end

if include_current_plot2
    nexttile;
    scatter(all_I_s,centered_error,sz,'filled')
    p = polyfit(all_I_s, centered_error, 1);
    fitted_values = polyval(p, all_I_s);
    hold on;
    grid on
    %lfit=plot(all_I_s, fitted_values, 'r-', 'LineWidth', 2); 
    mask = abs(all_I_s) < 15;
    filtered_I = all_I_s(mask);
    filtered_error = centered_error(mask);
    p_filtered = polyfit(filtered_I, filtered_error, 1);
    fitted_values_filtered = polyval(p_filtered, filtered_I);
    %lfit_small=plot(filtered_I, fitted_values_filtered, 'g-', 'LineWidth', 2); 
    %legend(lfit, 'Linear Fit', 'interpreter', 'latex', 'Location', 'southeast');
    
    %legend([lfit, lfit_small], {'Linear Fit (all data)', 'Linear Fit $|I| < 15 A$'}, 'interpreter', 'latex', 'Location', 'southeast');
    
    xlabel('Current [A]','FontSize' ,fs, 'interpreter', 'latex');
    ylabel('Error [mV]','FontSize' ,fs, 'interpreter', 'latex');
    title('Voltage Error over current','FontSize' ,fs, 'interpreter', 'latex');
    xlim([min(all_I_s) max(all_I_s)]);
end

if include_temperature_plot2
    nexttile;
    scatter(all_T_s(1:2600)-273.15,centered_error(1:2600),sz,'filled')
    p = polyfit(all_T_s, centered_error, 1);
    fitted_values = polyval(p, all_T_s);
    hold on;
    grid on
    lfit=plot(all_T_s-273.15, fitted_values, 'r-', 'LineWidth', 2); 
    legend(lfit,'Linear Fit', 'interpreter', 'latex', 'Location', 'southeast');
    xlabel('Temperature [$^{\circ}$C]','FontSize' ,fs, 'interpreter', 'latex');
    ylabel('Error [mV]','FontSize' ,fs, 'interpreter', 'latex');
    title('Voltage Error over Temperature','FontSize' ,fs, 'interpreter', 'latex');
    xlim([min(all_T_s-273.15) max(all_T_s-273.15)]);
    ylim([-0.1,0.1])
end
%%



fprintf('\nTemperature adjusted ECM: RMSE : %3.2f mV',RMSE{i+1})
if ~static
    RMSE_notempty=cell2mat(RMSE);
    figure()
    plot(NumericTemperatures-273.15,RMSE_notempty(1:end-1),'o');
    yline(RMSE_notempty(end),'-','Temperature-adjusted RMSE','FontSize' ,fs/1.5, 'interpreter', 'latex')
    ylabel('RMSE [mV]','FontSize' ,fs, 'interpreter', 'latex')
    xlabel('Temperature [$^{\circ}$C]','FontSize' ,fs, 'interpreter', 'latex')
    ylim([0 inf]);
    title('RMSE','FontSize' ,fs, 'interpreter', 'latex')
end

%%


% whiteness_test(centered_error)
% %%
% [acf, lags, bounds] = autocorr(centered_error, 'NumLags', 20); % You can adjust the number of lags
% figure;
% stem(lags, acf);
% hold on;
% plot(lags, bounds(1)*ones(size(lags)), 'r--'); % Upper confidence bound
% plot(lags, bounds(2)*ones(size(lags)), 'r--'); % Lower confidence bound
% title('Autocorrelation of Residuals');
% xlabel('Lags');
% ylabel('Autocorrelation');
% hold off;
