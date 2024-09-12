function error = cost_function_PM(params, MapCell, Scale_limit, Translation_limit,t_measured,V_measured,I_measured,SOC_measured,OneC_measured,data2,initialScale)
    
    scale= (params(1) * (Scale_limit(2) - Scale_limit(1)) + Scale_limit(1));
    translation= params(2) * (Translation_limit(2) - Translation_limit(1)) + Translation_limit(1);
    %scaleR1= (params(3) * (Scale_limitR1(2) - Scale_limitR1(1)) + Scale_limitR1(1));
    %translationR1= params(4) * (Translation_limitR1(2) - Translation_limitR1(1)) + Translation_limitR1(1);


    
    R0=MapCell{1};
    R0_scaled=R0*scale+translation;

    %for SOC and I Range 
    FilepathR0 = '.\ParameterMaps';
    FilesR0 = dir(fullfile(FilepathR0, 'ParameterMap_Chen_Marquis_diff_conc_radius_mod_OCV_Marquis.mat'));
    fullFileNameR0 = fullfile(FilepathR0, FilesR0.name);  
    load(fullFileNameR0); 

    %For R1, C1
    FilepathR1 = '.\Results\HPPC pybamm Chen';
    FilesR1 = dir(fullfile(FilepathR1, 'Parameters_pybamm_Chen_OCV_Marquis_combined.mat'));
    fullFileNameR1 = fullfile(FilepathR1, FilesR1.name);  
    load(fullFileNameR1, 'R1d', 'C1d'); 

    R1d_scaled=R1d;%*scaleR1+translationR1;

    V_ECM = GlobalModel_piecewise_function_PM_scaling(R0_scaled,data2,I_measured,SOC_measured,t_measured,OneC_measured,current_cd,SOC, R1d_scaled, C1d);
    save('V_ECM_optimized.mat',"V_ECM")
    if int32(scale)==initialScale %a way to save base V_ECM
        save('V_ECM_0.mat',"V_ECM")
    end
    % figure(10)
    % hold on
    % plot(t_measured,V_ECM)
    % plot(t_measured,V_measured)
    % legend({'Chen ECM','Marquis DFN'})
    
    % figure(11)
    % hold on
    error = V_measured-V_ECM;
    meanerror=mean(error)
    % plot(t_measured,error+0.1) 
    
    %error = abs( (scale * MapCell{1} + translation) - MapCell{2} );
end