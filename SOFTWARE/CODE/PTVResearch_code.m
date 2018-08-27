classdef PTVResearch < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PTVResearchUIFigure            matlab.ui.Figure
        LoadrawvideoButton             matlab.ui.control.Button
        RUNButton                      matlab.ui.control.Button
        TestButton                     matlab.ui.control.Button
        Slider                         matlab.ui.control.Slider
        CurrentFrameEditFieldLabel     matlab.ui.control.Label
        CurrentFrameEditField          matlab.ui.control.NumericEditField
        DisplayResultsSwitchLabel      matlab.ui.control.Label
        DisplayResultsSwitch           matlab.ui.control.RockerSwitch
        DisplayVariableDropDown        matlab.ui.control.DropDown
        STOPButton                     matlab.ui.control.Button
        PROLabel                       matlab.ui.control.Label
        TabGroup                       matlab.ui.container.TabGroup
        PreTab                         matlab.ui.container.Tab
        BackgroundsubtractionButton    matlab.ui.control.StateButton
        ofimagesformeanEditFieldLabel  matlab.ui.control.Label
        ofimagesformeanEditField       matlab.ui.control.NumericEditField
        SharpenButton                  matlab.ui.control.StateButton
        amountEditFieldLabel           matlab.ui.control.Label
        amountEditField                matlab.ui.control.NumericEditField
        MaskButton                     matlab.ui.control.Button
        ImageDewarping_________________Label  matlab.ui.control.Label
        ResetMaskButton                matlab.ui.control.Button
        ImagePreprocessing________________Label  matlab.ui.control.Label
        WienerFilterButton             matlab.ui.control.StateButton
        kernalsizeEditFieldLabel       matlab.ui.control.Label
        kernalsizeEditField            matlab.ui.control.NumericEditField
        amountEditField_2Label         matlab.ui.control.Label
        amountEditField_2              matlab.ui.control.NumericEditField
        thresholdEditFieldLabel        matlab.ui.control.Label
        thresholdEditField             matlab.ui.control.NumericEditField
        LocalIntensityButton           matlab.ui.control.StateButton
        LoadMaskButton                 matlab.ui.control.Button
        LoadparamsfrompreviousButton   matlab.ui.control.Button
        General______________________Label  matlab.ui.control.Label
        ImageMasking___________________Label_2  matlab.ui.control.Label
        ImageDewarpingButton           matlab.ui.control.StateButton
        GridspacingmmEditFieldLabel    matlab.ui.control.Label
        GridspacingmmEditField         matlab.ui.control.NumericEditField
        MethodDropDownLabel            matlab.ui.control.Label
        MethodDropDown                 matlab.ui.control.DropDown
        ProcessTab                     matlab.ui.container.Tab
        GridSizeLabel                  matlab.ui.control.Label
        MedianFilterButton             matlab.ui.control.StateButton
        xEditFieldLabel                matlab.ui.control.Label
        xEditField                     matlab.ui.control.NumericEditField
        yEditFieldLabel                matlab.ui.control.Label
        yEditField                     matlab.ui.control.NumericEditField
        PxmmEditFieldLabel             matlab.ui.control.Label
        PxmmEditField                  matlab.ui.control.NumericEditField
        offramessecLabel               matlab.ui.control.Label
        offramessecEditField           matlab.ui.control.NumericEditField
        InterpolatorDropDownLabel      matlab.ui.control.Label
        InterpolatorDropDown           matlab.ui.control.DropDown
        KernelsizeEditFieldLabel       matlab.ui.control.Label
        KernelsizeEditField            matlab.ui.control.NumericEditField
        CalculateLagrangianTracksButton  matlab.ui.control.StateButton
        ScalingParameters_______________Label  matlab.ui.control.Label
        LagrangianParameters______________Label  matlab.ui.control.Label
        GriddedParameters_______________Label  matlab.ui.control.Label
        ScalefromreferenceButton       matlab.ui.control.Button
        ProcessingMode_________________Label  matlab.ui.control.Label
        PTVButton                      matlab.ui.control.StateButton
        SnapshotPIVButton              matlab.ui.control.StateButton
        Label_5                        matlab.ui.control.Label
        PostGridTab                    matlab.ui.container.Tab
        PODDEMFilterButton             matlab.ui.control.Button
        PODButton                      matlab.ui.control.Button
        MeanFieldsButton               matlab.ui.control.StateButton
        RANSQuantitiesButton           matlab.ui.control.StateButton
        SwirlVorticityButton           matlab.ui.control.StateButton
        Postvariable_mean              matlab.ui.control.DropDown
        DisplayButton_2                matlab.ui.control.Button
        ResetButton                    matlab.ui.control.Button
        KernelEditFieldLabel           matlab.ui.control.Label
        KernelEditField                matlab.ui.control.NumericEditField
        ShowDataLabel                  matlab.ui.control.Label
        SmoothingButton                matlab.ui.control.Button
        DMDButton                      matlab.ui.control.Button
        Label_4                        matlab.ui.control.Label
        QuantitiestoCalculate_______________Label  matlab.ui.control.Label
        ModalDecompositions_______________Label  matlab.ui.control.Label
        Label_3                        matlab.ui.control.Label
        Filters______________________Label  matlab.ui.control.Label
        PostLagTab                     matlab.ui.container.Tab
        MinpathlengthEditFieldLabel    matlab.ui.control.Label
        MinpathlengthEditField         matlab.ui.control.NumericEditField
        trackstodisplayEditFieldLabel  matlab.ui.control.Label
        trackstodisplayEditField       matlab.ui.control.NumericEditField
        DisplayButton                  matlab.ui.control.Button
        LagrangianPathsButton          matlab.ui.control.StateButton
        DiffusionCoefficientsButton    matlab.ui.control.StateButton
        QuantitiestoCalculate_______________Label_2  matlab.ui.control.Label
        PostvariableLag                matlab.ui.control.DropDown
        ShowDataLabelLag               matlab.ui.control.Label
        LoadvectorfieldsButton         matlab.ui.control.Button
        Postvariable                   matlab.ui.control.DropDown
        Slider_pod                     matlab.ui.control.Slider
        DisplayVariableLabel           matlab.ui.control.Label
        NumbertoprocessEditFieldLabel  matlab.ui.control.Label
        NumbertoprocessEditField       matlab.ui.control.NumericEditField
        StartFrameEditFieldLabel       matlab.ui.control.Label
        StartFrameEditField            matlab.ui.control.NumericEditField
        EndFrameEditField              matlab.ui.control.NumericEditField
        EndFrameEditFieldLabel         matlab.ui.control.Label
        Forward                        matlab.ui.control.Button
        Backward                       matlab.ui.control.Button
        PlayButton                     matlab.ui.control.Button
        PauseButton                    matlab.ui.control.Button
        ColourbarLabel                 matlab.ui.control.Label
        DropDown                       matlab.ui.control.DropDown
        minEditFieldLabel              matlab.ui.control.Label
        minEditField                   matlab.ui.control.NumericEditField
        maxEditFieldLabel              matlab.ui.control.Label
        maxEditField                   matlab.ui.control.NumericEditField
        PlaybackSpeedSliderLabel       matlab.ui.control.Label
        PlaybackSpeedSlider            matlab.ui.control.Slider
        LoadrawimagesButton            matlab.ui.control.Button
        Lagvariable                    matlab.ui.control.DropDown
        UIAxes                         matlab.ui.control.UIAxes
    end


    properties (Access = public)
        Property % Description
        vid
        VidName
        VidPath
        x
        y
        sx
        sy
        M
        flag
        width 
        height 
        post
        Uz
        Vz
        Um
        Vm
        vecpath
        RANS
        VAR
        MEAN
        Orig
        P_u
        P_v
        S
        A
        Names
        Path
        mask
        pause
        speed
        lag
        
        work_post
        work_pre
        pro_lims;
        lag_lims;
        
        run_flag_pro
        prev_pro
        run_flag_post
        vec_flag
        
        ul 
        vl 
        xl 
        yl
        
        images_raw
        tform
        warp

    end

 

    methods (Access = public)        
    end

    methods (Access = private)
    
        function [U,V] = print_results(app,points_old,u,v,frame,f_num,P1,P2);
            V=fliplr(rot90((reshape(griddata(double(points_old(:,1)),double(points_old(:,2)),...
                    double(v),app.y,app.x,app.InterpolatorDropDown.Value),app.sx,app.sy)),-1));
            U=fliplr(rot90((reshape(griddata(double(points_old(:,1)),double(points_old(:,2)),...
                    double(u),app.y,app.x,app.InterpolatorDropDown.Value),app.sx,app.sy)),-1));
            tmpmask = round(sum(app.mask,3))/3;
            tmpmask=imresize(tmpmask,[size(U,1),size(U,2)],'nearest');
            tmpmask(tmpmask==0)=nan;
            U=U+single(tmpmask);
            V=V+single(tmpmask);    
            
            if app.MedianFilterButton.Value==1;
                U=medfilt2(U,[app.KernelsizeEditField.Value,app.KernelsizeEditField.Value]);
                V=medfilt2(V,[app.KernelsizeEditField.Value,app.KernelsizeEditField.Value]);
            end
            if strcmp(app.DisplayResultsSwitch.Value,'On');
                if strcmp(app.DisplayVariableDropDown.Value,'Vector Length');
                    imagesc(app.UIAxes,sqrt(U.^2+V.^2)/2);
                    colorbar(app.UIAxes)
                    app.UIAxes.XLim=[1 size(U,2)];
                    app.UIAxes.YLim=[1 size(U,1)];
                end
                if strcmp(app.DisplayVariableDropDown.Value,'U')
                    imagesc(app.UIAxes,U);
                    colorbar(app.UIAxes)
                    app.UIAxes.XLim=[1 size(U,2)];
                    app.UIAxes.YLim=[1 size(U,1)];
                end
                if strcmp(app.DisplayVariableDropDown.Value,'V')
                    imagesc(app.UIAxes,V);
                    colorbar(app.UIAxes)
                    app.UIAxes.XLim=[1 size(U,2)];
                    app.UIAxes.YLim=[1 size(U,1)];
                end
                if strcmp(app.DisplayVariableDropDown.Value,'Quiver')
                    quiver(app.UIAxes,U,V,2,'k');
                    colorbar(app.UIAxes,'off')
                    app.UIAxes.XLim=[1 size(U,2)];
                    app.UIAxes.YLim=[1 size(U,1)];
                end
                if strcmp(app.DisplayVariableDropDown.Value,'Points')
                    frame = read(app.vid,[f_num]);
                    if app.ImageDewarpingButton.Value==1;
                        frame = imwarp(frame,app.tform);            
                    end
                    imagesc(app.UIAxes,frame);
                    vel = sqrt(u.^2+v.^2);
                    x = P1(:,1);
                    y = P1(:,2);
                    hold(app.UIAxes,'on')
                    scatter(app.UIAxes,x(vel~=0),y(vel~=0),4,vel(vel~=0));
                    app.UIAxes.XLim=[1 size(frame,2)];
                    app.UIAxes.YLim=[1 size(frame,1)];
%                     caxis(app.UIAxes,([1.2*min(vel(vel~=0)),0.8*max(vel(vel~=0))]));
                    colorbar(app.UIAxes)
                end
                if strcmp(app.DropDown.Value,'auto')
                    caxis(app.UIAxes,'auto');
                    tmp=app.UIAxes.CLim; 
                    app.minEditField.Value=tmp(1);
                    app.maxEditField.Value=tmp(2);
                end
                if strcmp(app.DropDown.Value,'manual')
                    caxis(app.UIAxes,[app.minEditField.Value,app.maxEditField.Value]);
                end
                hold(app.UIAxes,'off')
                colormap(app.UIAxes,jet);
            else
            end 
        end
        
    end


    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            imagesc(app.UIAxes,imread('loading.jpg'));
            axis(app.UIAxes,'tight');
%             app.UIAxes.YTick=[];
%             app.UIAxes.XTick=[];
%  
            set( findall(app.PTVResearchUIFigure, '-property', 'Enable'), 'Enable', 'off')
            app.LoadrawvideoButton.Enable='on';
            app.LoadrawimagesButton.Enable='on';
            app.LoadvectorfieldsButton.Enable='on';
            app.Slider_pod.Visible='off';

        end

        % Button pushed function: LoadrawvideoButton
        function LoadrawvideoButtonPushed(app, event)

            if isempty(app.images_raw)
            [file,path] = uigetfile([{'*.mov';'*.avi';'*.mp4'}]);
                if file == 0
            	   return;
                end
            else
                path = app.images_raw.path;
                file = 'files.avi';
            end
            app.vec_flag=0;
            app.images_raw=[];
            app.work_post=[];
            app.VidName = file; 
            app.VidPath = path;  
            app.tform =[];
            app.vid = VideoReader([path,file]);
            app.work_pre=1;
            app.UIAxes.YDir='normal';
            app.offramessecEditField.Value=app.vid.FrameRate;
            %begins new session
            nums=1;
            flag=0; 
            [~,out]=fileparts([app.vid.Name]);
            fn=sprintf('%s/%s_%0.4d',app.vid.Path,out,nums);
            if ~exist(fn);
                mkdir(fn);
            else 
                while flag==0; 
                    nums = nums + 1;
                    fn=sprintf('%s/%s_%0.4d',app.vid.Path,out,nums);
                    if ~exist(fn);
                        flag=1;
                    end
                end
                mkdir(fn);
            end
            app.Path = fn;
            if isempty(file); 
               return 
            end
            
            if isempty(app.mask);
                app.mask = uint8(ones(size(read(app.vid,1))));
            end
            if ~isempty(app.vid)
                app.DisplayVariableDropDown.Visible='on';
                app.Postvariable.Visible='off';
                app.UIAxes.Visible='on';
                imagesc(app.UIAxes,app.mask.*read(app.vid,[app.work_pre]))
                app.UIAxes.YDir='normal'
                colormap(app.UIAxes,gray);
                colorbar(app.UIAxes,'off')
                app.UIAxes.XLim=[1 app.vid.Width];
                app.UIAxes.YLim=[1 app.vid.Height];
                app.Slider.Limits=[1 app.vid.Duration*app.vid.FrameRate];
                app.width = app.vid.Width/8;
                app.height = app.vid.Height/8;
                [x y]=meshgrid(linspace(1,app.vid.Height,round(app.height)),linspace(1,app.vid.Width,round(app.width)));
                app.sx = size(x,1);app.sy=size(x,2);
                app.x = x(:); app.y = y(:);
                app.xEditField.Value=round(app.vid.Width/8);
                app.yEditField.Value=round(app.vid.Height/8);
                set(findall(app.PTVResearchUIFigure, '-property', 'Enable'), 'Enable', 'on')
                app.Lagvariable.Visible='off';
            end
            StartFrameEditFieldValueChanged(app, event)

        end

        % Value changed function: Slider
        function SliderValueChanged(app, event)
            value = app.Slider.Value;
            app.CurrentFrameEditField.Value=round(value);
            app.UIAxes.YDir='normal';
%             NumbertoprocessEditFieldValueChanged(app, event);
            
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -');
                app.Slider.Limits=app.pro_lims;
                U=app.Uz(:,:,round(value)-app.pro_lims(1)+1);
                V=app.Vz(:,:,round(value)-app.pro_lims(1)+1);
                
                u = U-app.Um;
                v = V-app.Vm;
                
%                 tmpmask = round(sum(app.mask,3))/3;
%                 tmpmask=imresize(tmpmask,[size(U,1),size(U,2)],'nearest');
%                 tmpmask(tmpmask==0)=nan;
%                 U=U+tmpmask;
%                 V=V+tmpmask;
%                 u=u+tmpmask;
%                 v=v+tmpmask;

                tmp = size(app.Uz,3);
                if strcmp(app.Postvariable.Value,'U')
                    imagesc(app.UIAxes,U);
                    app.UIAxes.YDir='normal'
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'V')
                    imagesc(app.UIAxes,V);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Vorticity');
                    imagesc(app.UIAxes,curl(U,V));
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Vector Length');
                    imagesc(app.UIAxes,sqrt(U.^2+V.^2)/2);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'uu');
                    imagesc(app.UIAxes,u.*u);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'uv');
                    imagesc(app.UIAxes,u.*v);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'vv');
                    imagesc(app.UIAxes,v.*v);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'TKE');
                    imagesc(app.UIAxes,(v.^2.+v.^2)/2);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Fluc Vorticity');
                    imagesc(app.UIAxes,curl(U,V));
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Swirl');
                    r = size(u,1);
                    c= size(u,2);
                    [du_dx,du_dy]=gradient(u);
                    [dv_dx,dv_dy]=gradient(v);                
                    Q=(du_dy.*dv_dx).*(-0.5); % Coherent structures identifier done using Q-criterion
                    lam=zeros(r,c);
                    for w=1:r;
                        for j=1:c;
                            D=[du_dx(w,j) du_dy(w,j);dv_dx(w,j) dv_dy(w,j)];
                            if isnan(det(D))~=1;
                               lamtemp=imag(eig(D));
                               pos=find(lamtemp>0);
                                if length(pos)~=0;
                                lam(w,j)=lamtemp(pos);
                                end
                            end
                        end
                    end
                    imagesc(app.UIAxes,lam);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Quiver')
                    quiver(app.UIAxes,U,V,5,'k');
                    colorbar(app.UIAxes,'off')
                end
                app.UIAxes.XLim=[1 size(U,2)];
                app.UIAxes.YLim=[1 size(U,1)];
                colormap(app.UIAxes,jet);
                app.work_post=value
            
            elseif strcmp(app.TabGroup.SelectedTab.Title,'Pre -') | strcmp(app.TabGroup.SelectedTab.Title,'Process')
                app.Slider.Limits=[1 round(app.vid.Duration.*app.vid.FrameRate)];
                frame = read(app.vid,round(value));
                if size(frame,3)>1
                    frame = rgb2gray(frame);
                end
                if app.BackgroundsubtractionButton.Value==1;
                    frame = single(frame)-app.M;
                end
                frame = uint8(app.mask(:,:,1)).*uint8(frame); 
                if app.SharpenButton.Value==1;
                    frame=imsharpen(frame,'Amount',app.amountEditField.Value);
                end
                if app.WienerFilterButton.Value==1;
                    frame=wiener2(frame,[app.kernalsizeEditField.Value,app.kernalsizeEditField.Value]);
                end
                if app.LocalIntensityButton.Value==1;
                    frame=localcontrast(frame,app.thresholdEditField.Value,app.amountEditField_2.Value);
                end
                if app.ImageDewarpingButton.Value==1;
                    frame = imwarp(frame,app.tform);            
                end

                imagesc(app.UIAxes,(frame));
                app.work_pre = round(value);
                colorbar(app.UIAxes,'off');
                colormap(app.UIAxes,gray);
%                 axis(app.UIAxes,'tight');
                app.UIAxes.XLim=[1 size(frame,2)];
                app.UIAxes.YLim=[1 size(frame,1)];
%                 
            elseif  strcmp(app.TabGroup.SelectedTab.Title,'Post Lag -')
                    cla(app.UIAxes);
                    app.UIAxes.YDir='normal';
                    app.Slider.Limits=app.lag_lims;
                    if strcmp(app.Lagvariable.Value,'Points vec length')
                        scatter(app.UIAxes,app.xl{round(value)},app.yl{round(value)},4,(abs(app.ul{round(value)})+abs(app.vl{round(value)}))/2);
                    end
                    if strcmp(app.Lagvariable.Value,'Points u-vel')
                        scatter(app.UIAxes,app.xl{round(value)},app.yl{round(value)},4,app.ul{round(value)});
                    end
                    if strcmp(app.Lagvariable.Value,'Points v-vel')
                        scatter(app.UIAxes,app.xl{round(value)},app.yl{round(value)},4,app.vl{round(value)});
                    end
                    colormap(app.UIAxes,jet);
                    app.UIAxes.XLim=[1 max(app.xl{round(value)})];
                    app.UIAxes.YLim=[1 max(app.yl{round(value)})];
            end
            if strcmp(app.DropDown.Value,'auto')
                tmp=app.UIAxes.CLim;
                caxis(app.UIAxes,'auto');
                app.minEditField.Value=tmp(1);
                app.maxEditField.Value=tmp(2);
            end
            if strcmp(app.DropDown.Value,'manual')
                caxis(app.UIAxes,[app.minEditField.Value,app.maxEditField.Value]);
            end
            
            app.UIAxes.YDir='normal';

            
        end

        % Value changed function: CurrentFrameEditField
        function CurrentFrameEditFieldValueChanged(app, event)
            value = app.CurrentFrameEditField.Value;
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -')  ;
                if app.CurrentFrameEditField.Value < app.pro_lims(2)
                    value = app.CurrentFrameEditField.Value;
                    app.Slider.Value=value;
                    app.work_post=value;
                    SliderValueChanged(app, event)
                else
                    app.CurrentFrameEditField.Value= app.pro_lims(2)
                    value = app.CurrentFrameEditField.Value;
                    app.Slider.Value=value;
                    app.work_post=value;
                    SliderValueChanged(app, event)
                end
            elseif strcmp(app.TabGroup.SelectedTab.Title,'Post Lag -')
                if app.CurrentFrameEditField.Value < app.lag_lims(2)
                    value = app.CurrentFrameEditField.Value;
                    app.Slider.Value=value;
                    app.work_post=value;
                    SliderValueChanged(app, event)
                else
                    app.CurrentFrameEditField.Value= app.lag_lims(2)
                    value = app.CurrentFrameEditField.Value;
                    app.Slider.Value=value;
                    app.work_post=value;
                    SliderValueChanged(app, event)
                end
            else 
                 if app.CurrentFrameEditField.Value < app.vid.Duration.*app.vid.FrameRate;
                    value = app.CurrentFrameEditField.Value;
                    app.Slider.Value=value;
                    app.work_pre=value;
                    SliderValueChanged(app, event)
                else
                    app.CurrentFrameEditField.Value=round(app.vid.Duration.*app.vid.FrameRate)
                    value = app.CurrentFrameEditField.Value;
                    app.Slider.Value=value;
                    app.work_pre=value;
                    SliderValueChanged(app, event)
                 end
            end         
        end

        % Button pushed function: RUNButton
        function RUNButtonPushed(app, event)
             
            if strcmp(app.TabGroup.SelectedTab.Title,'Pre -') | strcmp(app.TabGroup.SelectedTab.Title,'Process')           
                app.RANS=[];
                app.VAR=[];
                app.MEAN=[];
                app.run_flag_pro=1;
                app.prev_pro=app.TabGroup.SelectedTab.Title;
                app.work_post=[];
                f_num = app.StartFrameEditField.Value;
                app.flag=0;
            
                %pre_processing
                app.Path
                mkdir([app.Path,'/PARAMS'])
                fileID = fopen([app.Path,'/PARAMS/params.txt'],'w');
                fprintf(fileID,['app.StartFrameEditField.Value=',num2str(app.StartFrameEditField.Value),'\n']);
                fprintf(fileID,['app.NumbertoprocessEditField.Value=',num2str(app.NumbertoprocessEditField.Value),'\n']);
                fprintf(fileID,['app.EndFrameEditField.Value=',num2str(app.EndFrameEditField.Value),'\n']);
                if app.BackgroundsubtractionButton.Value==0;
                    fprintf(fileID,'app.BackgroundsubtractionButton.Value=0\n')
                else
                    fprintf(fileID,'app.BackgroundsubtractionButton.Value=1\n')
                    fprintf(fileID,['app.ofimagesformeanEditField.Value=',num2str(app.ofimagesformeanEditField.Value),'\n'])
                end
                if app.SharpenButton.Value==0;
                    fprintf(fileID,'app.SharpenButton.Value=0\n')
                else
                    fprintf(fileID,'app.SharpenButton.Value=1\n')
                    fprintf(fileID,['app.amountEditField.Value=',num2str(app.amountEditField.Value),'\n'])
                end
                if app.WienerFilterButton.Value==0;
                    fprintf(fileID,'app.WienerFilterButton.Value=0\n')
                else
                    fprintf(fileID,'app.WienerFilterButton.Value=1\n')
                    fprintf(fileID,['app.kernalsizeEditField.Value=',num2str(app.kernalsizeEditField.Value),'\n'])
                end
                if app.LocalIntensityButton.Value==0;
                    fprintf(fileID,'app.LocalIntensityButton.Value=0\n')
                else
                    fprintf(fileID,'app.WienerFilterButton.Value=1\n')
                    fprintf(fileID,['app.thresholdEditField.Value=',num2str(app.thresholdEditField.Value),'\n'])
                    fprintf(fileID,['app.amountEditField.Value=',num2str(app.amountEditField.Value),'\n'])
                end
                fprintf(fileID,['app.PxmmEditField.Value=',num2str(app.PxmmEditField.Value),'\n'])
                fprintf(fileID,['app.offramessecEditField.Value=',num2str(app.offramessecEditField.Value),'\n'])
                if app.CalculateLagrangianTracksButton.Value==0;
                    fprintf(fileID,'app.CalculateLagrangianTracksButton.Value=0\n')
                else
                    fprintf(fileID,'app.CalculateLagrangianTracksButton.Value=1\n')
                end
                fprintf(fileID,['app.xEditField.Value=',num2str(app.xEditField.Value),'\n'])
                fprintf(fileID,['app.yEditField.Value=',num2str(app.yEditField.Value),'\n'])
                fprintf(fileID,['app.InterpolatorDropDown.Value=','''',app.InterpolatorDropDown.Value,'''','\n'])
                if app.MedianFilterButton.Value==0;
                    fprintf(fileID,'app.MedianFilterButton.Value=0\n')
                else
                    fprintf(fileID,'app.MedianFilterButton.Value=1\n')
                    fprintf(fileID,['app.KernelEditField.Value=',num2str(app.KernelEditField.Value),'\n'])
                end
                if app.MedianFilterButton.Value==0;
                    fprintf(fileID,'app.MedianFilterButton.Value=0\n')
                end
                if app.PTVButton.Value==1
                    fprintf(fileID,'app.PTVButton.Value=1\n')
                    fprintf(fileID,'app.SnapshotPIVButton.Value=0\n')
                else
                    fprintf(fileID,'app.PTVButton.Value=0\n')
                    fprintf(fileID,'app.SnapshotPIVButton.Value=1\n')
                end
                if app.SnapshotPIVButton.Value==1
                    fprintf(fileID,'app.PTVButton.Value=0\n')
                    fprintf(fileID,'app.SnapshotPIVButton.Value=1\n')
                else
                    fprintf(fileID,'app.PTVButton.Value=1\n')
                    fprintf(fileID,'app.SnapshotPIVButton.Value=0\n')
                end
                if app.ImageDewarpingButton.Value==1; 
                    fprintf(fileID,'app.ImageDewarpingButton.Value=1\n')
                    fprintf(fileID,['app.MethodDropDown.Value=','''',app.MethodDropDown.Value,'''','\n'])
                    fprintf(fileID,['app.GridspacingmmEditField.Value=',num2str(app.GridspacingmmEditField.Value),'\n'])
                end
                
                fclose(fileID);
                
                if ~isempty(dir([app.Path,'/*.mat']))
                    [~,out]=fileparts([app.vid.Name]);
                    a = app.Path;
                    fn=sprintf('%s/%s_%0.4d',app.vid.Path,out,str2num(a(end-3:end))+1);
                    mkdir(fn);
                    mkdir([fn,'/PARAMS']);
                    copyfile([app.Path,'/PARAMS'],[fn,'/PARAMS']);
                    app.Path = fn;
                end
                
                if app.SnapshotPIVButton.Value==1;
                    U=[];
                    V=[];
                    app.Uz=[];
                    app.Vz=[]
                    count = 0;
                    start = f_num;
                    app.DisplayVariableDropDown.Value='U'
                    i = 0;
                    for loop = app.StartFrameEditField.Value:(app.EndFrameEditField.Value)/2;
                        i = i+1;
                        count = count + 1; 
                        if app.flag==1;
                            break
                        end
                        frame = app.mask.*read(app.vid,[f_num]);
                        if size(frame,3)>1
                            frame = app.mask(:,:,1).*rgb2gray(frame);
                        end
                        if app.BackgroundsubtractionButton.Value==1;
                            frame = single(app.mask(:,:,1)).*single(frame)-single(app.mask(:,:,1)).*app.M;
                        end         
                        if app.SharpenButton.Value==1;
                            frame=imsharpen(frame,'Amount',app.amountEditField.Value);
                        end
                        if app.WienerFilterButton.Value==1;
                            frame=wiener2(frame,[app.kernalsizeEditField.Value,app.kernalsizeEditField.Value]);
                        end
                        if app.LocalIntensityButton.Value==1;
                            frame=localcontrast(frame,app.thresholdEditField.Value,app.amountEditField_2.Value);
                        end
                        if app.ImageDewarpingButton.Value==1
                            frame = imwarp(frame,app.tform);            
                        end

                        points = detectMinEigenFeatures(frame);
                        points_old =  points.Location;
                        tracker = vision.PointTracker('MaxBidirectionalError',5,'MaxIterations',50);
                        initialize(tracker,points.Location,frame);
                        
                        f_num = f_num+1;
                        frame = app.mask.*read(app.vid,[f_num]);
                        if size(frame,3)>1
                            frame = app.mask(:,:,1).*rgb2gray(frame);
                        end
                        if app.BackgroundsubtractionButton.Value==1;
                            frame = single(app.mask(:,:,1)).*single(frame)-single(app.mask(:,:,1)).*app.M;
                        end         
                        if app.SharpenButton.Value==1;
                            frame=imsharpen(frame,'Amount',app.amountEditField.Value);
                        end
                        if app.WienerFilterButton.Value==1;
                            frame=wiener2(frame,[app.kernalsizeEditField.Value,app.kernalsizeEditField.Value]);
                        end
                        if app.LocalIntensityButton.Value==1;
                            frame=localcontrast(frame,app.thresholdEditField.Value,app.amountEditField_2.Value);
                        end
                        if app.ImageDewarpingButton.Value==1
                            frame = imwarp(frame,app.tform);            
                        end
                        [points,validity] = tracker(frame);                    
                        u = (points_old(:,1)-points(:,1)).*(app.offramessecEditField.Value/app.PxmmEditField.Value);
                        v = (points_old(:,2)-points(:,2)).*(app.offramessecEditField.Value/app.PxmmEditField.Value);
                        release(tracker)
                        f_num = f_num+1;
                        P1=[];
                        P2=[];
                        [U,V] = print_results(app,points_old,u,v,frame,f_num,P1,P2);
                        u = U;
                        v = V;
                        drawnow
                        save(sprintf('%s/output_%0.6d.mat',app.Path,count),'u','v');
                        app.Slider.Value=f_num-1;
                        app.CurrentFrameEditField.Value=f_num-1;
                        app.PROLabel.FontColor='k';
                        app.PROLabel.FontSize=28;
                        app.PROLabel.Text=sprintf('%0.0f %%',200*(i/(app.EndFrameEditField.Value-app.StartFrameEditField.Value)));
                    end
                    app.PROLabel.FontColor=[0.94 0.94 0.94];
                    app.run_flag_pro=0;
                end
                
                if app.SnapshotPIVButton.Value==0;
                    frame = app.mask.*read(app.vid,[f_num]);
                    
                    frame = padarray(frame(2:end-1,2:end-1,:),[1 1]);
                    
                    if size(frame,3)>1
                        frame = app.mask(:,:,1).*rgb2gray(frame);
                    end
                    if app.BackgroundsubtractionButton.Value==1;
                        frame = single(app.mask(:,:,1)).*single(frame)-single(app.mask(:,:,1)).*app.M;
                    end         
                    if app.SharpenButton.Value==1;
                        frame=imsharpen(frame,'Amount',app.amountEditField.Value);
                    end
                    if app.WienerFilterButton.Value==1;
                        frame=wiener2(frame,[app.kernalsizeEditField.Value,app.kernalsizeEditField.Value]);
                    end
                    if app.LocalIntensityButton.Value==1;
                        frame=localcontrast(frame,app.thresholdEditField.Value,app.amountEditField_2.Value);
                    end
                    if app.ImageDewarpingButton.Value==1
                        frame = imwarp(frame,app.tform,'OutputView',imref2d(size(frame)));            
                    end
                    points = detectMinEigenFeatures(frame);
                    points_old =  points.Location;
                    tracker = vision.PointTracker('MaxBidirectionalError',5,'MaxIterations',50);
                    initialize(tracker,points.Location,frame);
                    U=[];
                    V=[];
                    app.Uz=[];
                    app.Vz=[]
                    tmp_num = app.NumbertoprocessEditField.Value;
                    start = f_num
                    i=0;
                    for loop = app.StartFrameEditField.Value:app.EndFrameEditField.Value-1;
                        i=i+1;
                        app.NumbertoprocessEditField.Value=tmp_num;
                        if app.flag==1;
                            break
                        end
                        f_num = f_num+1;
                        frame = app.mask.*read(app.vid,[f_num]);
                        if size(frame,3)>1
                            frame = app.mask(:,:,1).*rgb2gray(frame);
                        end
                        if app.BackgroundsubtractionButton.Value==1;
                            frame = single(app.mask(:,:,1)).*single(frame)-single(app.mask(:,:,1)).*app.M;
                        end         
                        if app.SharpenButton.Value==1;
                            frame=imsharpen(frame,'Amount',app.amountEditField.Value);
                        end
                        if app.WienerFilterButton.Value==1;
                            frame=wiener2(frame,[app.kernalsizeEditField.Value,app.kernalsizeEditField.Value]);
                        end
                        if app.LocalIntensityButton.Value==1;
                            frame=localcontrast(frame,app.thresholdEditField.Value,app.amountEditField_2.Value);
                        end
                        if app.ImageDewarpingButton.Value==1
                            frame = imwarp(frame,app.tform,'OutputView',imref2d(size(frame)));            
                        end
    %                     frame = single(uint8(app.mask(:,:,1)).*uint8(frame));
                        
                        [points,validity] = tracker(frame);                    
                        u = (points_old(:,1)-points(:,1)).*(app.offramessecEditField.Value/app.PxmmEditField.Value);
                        v = (points_old(:,2)-points(:,2)).*(app.offramessecEditField.Value/app.PxmmEditField.Value);
                        release(tracker)
                        if app.CalculateLagrangianTracksButton.Value==1;
                            P1=[points_old(:,1),points_old(:,2)];
                            P2=[points(:,1),points(:,2)];
                        end
                        if ~exist('P1')
                            P1=[];
                            P2=[];
                        end
                        [U,V] = print_results(app,points_old,u,v,frame,f_num,P1,P2);                
                        points_new = detectMinEigenFeatures(frame);
                        points_old = points;
                        [val,l]=ismembertol(points_new.Location,points_old,0.01,'ByRows',true); 
                        points = single([points_new.Location(val==0,:);points_old]);
                        tracker = vision.PointTracker('MaxBidirectionalError',5,'MaxIterations',50);
                        initialize(tracker,abs(points),frame);
                        points_old=points;
                        drawnow
                        % saving line
                        u = U;
                        v = V;
                        if app.CalculateLagrangianTracksButton.Value==1;
                            save(sprintf('%s/output_%0.6d.mat',app.Path,f_num-1),'u','v','P1','P2');
                        else
                            save(sprintf('%s/output_%0.6d.mat',app.Path,f_num-1),'u','v');
                        end
                        app.Slider.Value=f_num-1;
                        app.CurrentFrameEditField.Value=f_num-1;
                        app.PROLabel.FontColor='k';
                        app.PROLabel.FontSize=28;
                        
                        app.PROLabel.Text=sprintf('%0.0f %%',100*(i/(app.EndFrameEditField.Value-app.StartFrameEditField.Value)));
                        
                    end  
                    app.PROLabel.FontColor=[0.94 0.94 0.94];
                    app.run_flag_pro=0;
                end
            end
            
            %post grid  
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -') 
                app.run_flag_post=1;
                start_val=max([app.StartFrameEditField.Value-app.Slider.Limits(1)+1]);
                end_val=min([app.EndFrameEditField.Value-app.Slider.Limits(1)+1,app.Slider.Limits(2)-app.Slider.Limits(1)+1]);
                names = dir([app.Path,'/','*.mat']);
                
                if app.MeanFieldsButton.Value==1; 
                    MEAN=[];
                    Um=mean(single(app.Uz(:,:,start_val:end_val)),3);
                    Vm=mean(single(app.Vz(:,:,start_val:end_val)),3);
                    MEAN.U_mean=Um;
                    MEAN.V_mean=Vm;
                    MEAN.VEC_mean = sqrt(Um.^2+Vm.^2);
                    MEAN.Vorticity=curl(Um,Vm);
                    fn=sprintf('%s/STATS',app.Path);
                    if ~exist(fn);
                        mkdir(fn);
                    end
                    save(sprintf('%s/MEAN_STATS.mat',fn),'MEAN');
                    app.MEAN=MEAN;
                end
                
                if app.RANSQuantitiesButton.Value==1;
                    U = single(app.Uz(:,:,start_val:end_val));
                    V = single(app.Vz(:,:,start_val:end_val));
                    Um=mean(single(app.Uz(:,:,start_val:end_val)),3);
                    Vm=mean(single(app.Vz(:,:,start_val:end_val)),3);
                    u = U - repmat(Um,[1 1 size(U,3)]);
                    v = V - repmat(Vm,[1 1 size(V,3)]);
                    app.flag=0;
                    app.PROLabel.FontColor='k';
                    app.PROLabel.FontSize=28;
                    app.PROLabel.Text=sprintf('%0.0f %%',0);
                    fn=sprintf('%s/RANS',app.Path);
                    if ~exist(fn);
                        mkdir(fn);
                    end         
                    uu=[];uv=[];vv=[];tke=[];
                    tmpuu=[];tmpuv=[];tmpvv=[];tmptke=[];                    
                    i = 0;
                    for l = start_val:end_val;
                        i = i+1; 
                        if app.flag==1;
                            break
                        end
                        tmpuu(:,:,i)=u(:,:,i).*u(:,:,i);
                        tmpuv(:,:,i)=u(:,:,i).*v(:,:,i);
                        tmpvv(:,:,i)=v(:,:,i).*v(:,:,i);
                        tmptke(:,:,i)=(u(:,:,i).^2+v(:,:,i).^2)/2;
                        app.PROLabel.Text=sprintf('%0.0f %%',100*(i/size(U,3)));
                        drawnow
                        uu = tmpuu(:,:,i);
                        uv = tmpuv(:,:,i);
                        vv = tmpvv(:,:,i);
                        tke = tmptke(:,:,i); 
                        save(sprintf('%s/%s',fn,names(l).name),'uu','uv','vv','tke'); 
                    end
                    UU = mean(tmpuu,3);
                    UV = mean(tmpuv,3);
                    VV = mean(tmpvv,3);
                    TKE = mean(tmptke,3);           
                    save(sprintf('%s/RANS_STATS.mat',fn),'UU','UV','VV','TKE');
                    app.PROLabel.FontColor=[0.94 0.94 0.94];   
                    app.RANS.UU=UU;
                    app.RANS.UV=UV;
                    app.RANS.VV=VV;
                    app.RANS.TKE=TKE;
                    uu=[];uv=[];vv=[];tke=[];
                    tmpuu=[];tmpuv=[];tmpvv=[];tmptke=[]; 
                end
                
                if app.SwirlVorticityButton.Value==1;
                    U = single(app.Uz(:,:,start_val:end_val));
                    V = single(app.Vz(:,:,start_val:end_val));
                    Um=mean(single(app.Uz(:,:,start_val:end_val)),3);
                    Vm=mean(single(app.Vz(:,:,start_val:end_val)),3);
                    u = U - repmat(Um,[1 1 size(U,3)]);
                    v = V - repmat(Vm,[1 1 size(V,3)]);
                    app.flag=0;
                    app.PROLabel.FontColor='k';
                    app.PROLabel.FontSize=28;
                    app.PROLabel.Text=sprintf('%0.0f %%',0);
                    fn=sprintf('%s/MISC',app.Path);
                    if ~exist(fn);
                        mkdir(fn);
                    end         
                    i = 0
                    for l = start_val:end_val;
                        i = i+1;
                        if app.flag==1;
                            break
                        end
                        tmpvort=[];
                        tmpvort(:,:,i)=curl(u(:,:,i),v(:,:,i));
                        r = size(u(:,:,i),1);
                        c= size(u(:,:,i),2);
                        [du_dx,du_dy]=gradient(u(:,:,i));
                        [dv_dx,dv_dy]=gradient(v(:,:,i));                
                        Q=(du_dy.*dv_dx).*(-0.5); % Coherent structures identifier done using Q-criterion
                        lam=zeros(r,c);
                        for w=1:r;
                            for j=1:c;
                                D=[du_dx(w,j) du_dy(w,j);dv_dx(w,j) dv_dy(w,j)];
                                if isnan(det(D))~=1;
                                   lamtemp=imag(eig(D));
                                   pos=find(lamtemp>0);
                                    if length(pos)~=0;
                                        lam(w,j)=lamtemp(pos);
                                    end
                                end
                            end
                        end
                        tmpswirl(:,:,i)=lam;
                        swirl = lam;
                        vort = tmpvort(:,:,i);
                        save(sprintf('%s/%s',fn,names(l).name),'swirl','vort'); 
                        app.PROLabel.Text=sprintf('%0.0f %%',100*(i/size(U,3)));
                        drawnow
                    end
                    MEAN_VAR.VORT = mean(tmpvort,3);
                    MEAN_VAR.SWIRL = mean(tmpswirl,3);
                    VORT = mean(tmpvort,3);
                    SWIRL = mean(tmpswirl,3);
                    app.VAR=MEAN_VAR;
                    save(sprintf('%s/VAR_STATS.mat',fn),'VORT','SWIRL');
                    app.PROLabel.FontColor=[0.94 0.94 0.94];  
                end
                TabGroupSelectionChanged(app, event)
                DisplayButton_2Pushed(app, event)
                app.run_flag_post=0
            end
            
            
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Lag -') 
                if app.LagrangianPathsButton.Value==1 | app.DiffusionCoefficientsButton.Value ==1;
                    start = app.StartFrameEditField.Value;  
                    start_write = 1; 
                    PT=0;
                    mkdir([app.Path,'/','LANG']);
                    delete([app.Path,'/','LANG','/','*.mat']);
                    count = 0; 
                    names = dir([app.Path,'/','*.mat']);
                    load([app.Path,'/',names(1).name]);
                    pathsx=[P1(:,1),P2(:,1)];
                    pathsy=[P1(:,2),P2(:,2)]; 
                    for i = start+1:min([start+1+app.NumbertoprocessEditField.Value,size(cat(1,names.name),1)])+app.MinpathlengthEditField.Value+1;
                        app.PROLabel.FontColor='k';
                        app.PROLabel.FontSize=28;
                        count = count + 1;
                        if i <= min([start+1+app.NumbertoprocessEditField.Value,size(cat(1,names.name),1)])
                            load([app.Path,'/',names(i).name]);
                            [v,l]=ismembertol(P1,[pathsx(:,end),pathsy(:,end)],.001,'ByRows',true); 
                            p1 = P1(l==0,:);
                            p2 = P2(l==0,:);
                            P1=P1(v,:);
                            P2=P2(v,:);
                            l=l(v);                
                            pathsx(l,end)=[P1(:,1)];
                            pathsy(l,end)=[P1(:,2)];
                            pathsx(end+1:end+size(p1,1),end)=p1(:,1);
                            pathsy(end+1:end+size(p1,1),end)=p1(:,2);                
                            pathsx(l,end+1)=[P2(:,1)];
                            pathsy(l,end+1)=[P2(:,2)];
                            pathsx(end-size(p1,1)+1:end,end)=p2(:,1);
                            pathsy(end-size(p1,1)+1:end,end)=p2(:,2);
                            pathsx(pathsx==0)=nan;
                            pathsy(pathsy==0)=nan;
                        end
                        if count >= app.MinpathlengthEditField.Value+1 
                            tracks = [pathsx(:,1),pathsy(:,1)];
                            save(sprintf('%s/%s',[app.Path,'/','LANG'],names(i-(app.MinpathlengthEditField.Value+1)).name),'tracks'); 
                            PT=PT+1;
                            pathsx(:,1)=[];
                            pathsy(:,1)=[];
                        end 
                        app.PROLabel.Text=sprintf('%0.0f %%',100*(i/(min([start+1+app.NumbertoprocessEditField.Value,size(cat(1,names.name),1)])+app.MinpathlengthEditField.Value+1)));                        
                        drawnow
                        app.PROLabel.FontColor=[0.94 0.94 0.94];  
                    end        
                    if app.LagrangianPathsButton.Value==1 
                        DisplayButtonPushed(app, event)
                    end
                end
            end
%             
            if app.DiffusionCoefficientsButton.Value==1;
                mkdir([app.Path,'/','LANG/','STATS'])
                app.Path
                names = dir([app.Path,'/','LANG','/','*.mat'])
                load([app.Path,'/','LANG','/',names(1).name]);
                sty=[];
                stx=[];
                xs=[];
                ys=[];
                for i = 2:min([app.NumbertoprocessEditField.Value,size(cat(1,names.name),1)]); 
                   load([app.Path,'/','LANG','/',names(i-1).name]);
                   t1 = tracks;
                   load([app.Path,'/','LANG','/',names(i).name]);
                   t2 = tracks; 
                   velx = abs(t2(1:size(t1,1),1)-t1(1:size(t1,1),1)).*app.offramessecEditField.Value;    
                   vely = abs(t2(1:size(t1,2),1)-t1(1:size(t1,2),1)).*app.offramessecEditField.Value;    
                   if i == 2; 
                      stx = velx;
                      sty = vely;
                      xs=0;
                      ys=0;
                   else
                      velx(1:length(stx)) = velx(1:length(stx)) + stx;
                      vely(1:length(sty)) = vely(1:length(sty)) + sty;
                      xs(i-2)=sum(velx(~isnan(velx)));
                      ys(i-2)=sum(vely(~isnan(vely)));
                      stx = velx; 
                      sty = vely;
                   end
                end
                xdiff = xs; 
                ydiff = ys;
                save([app.Path,'/','LANG/','STATS/DIFF.mat'],'xdiff','ydiff')
            end

            if exist([app.Path,'/','LANG']);
                app.PostvariableLag.Items=[{'Tracks Vec Length','Tracks u-vel','Tracks v-vel'}]; 
            end
            if exist([app.Path,'/','LANG/','STATS']);
                app.PostvariableLag.Items=[{'Diff-x','Diff-y','Tracks Vec Length','Tracks u-vel','Tracks v-vel'}]; 
            end
        end

        % Value changed function: BackgroundsubtractionButton
        function BackgroundsubtractionButtonValueChanged(app, event)
            value = app.BackgroundsubtractionButton.Value;
            app.RUNButton.Enable='off';
            if value == 1;
                app.PROLabel.FontColor='k';
                app.PROLabel.FontSize=28;
                num = app.NumbertoprocessEditField.Value;
                tmp = read(app.vid,1);    
                Ms = zeros(size(tmp,1),size(tmp,2),1); 
                c=0;
                for i = app.StartFrameEditField.Value:app.StartFrameEditField.Value+num-1;
                    c=c+1;
                    app.PROLabel.Text=sprintf('%0.0f %%',100*(c/num));
                    M=read(app.vid,i);    
                    if size(M,3)==3;
                        Ms = Ms+double(rgb2gray(M));
                    else
                        Ms = Ms+double(M);
                    end
                end
                M=Ms/c;
                app.M=M;
            end
            app.RUNButton.Enable='on';
            app.PROLabel.FontColor=[0.94 0.94 0.94];
            SliderValueChanged(app, event)
        end

        % Button pushed function: TestButton
        function TestButtonPushed(app, event)
            tmp = app.NumbertoprocessEditField.Value;
            app.NumbertoprocessEditField.Value = 1;
            RUNButtonPushed(app)
            app.NumbertoprocessEditField.Value=tmp;
        end

        % Button pushed function: STOPButton
        function STOPButtonPushed(app, event)
            app.flag = 1;
        end

        % Value changed function: xEditField
        function xEditFieldValueChanged(app, event)
            value = app.xEditField.Value 
            app.width = value;
            [x y]=meshgrid(linspace(1,app.vid.Height,round(app.height)),linspace(1,app.vid.Width,round(app.width)));
            app.sx = size(x,1);app.sy=size(x,2);
            app.x = x(:); app.y = y(:);
        end

        % Value changed function: yEditField
        function yEditFieldValueChanged(app, event)
            value = app.xEditField.Value 
            app.height = value;
            [x y]=meshgrid(linspace(1,app.vid.Height,round(app.height)),linspace(1,app.vid.Width,round(app.width)));
            app.sx = size(x,1);app.sy=size(x,2);
            app.x = x(:); app.y = y(:);
        end

        % Button pushed function: LoadvectorfieldsButton
        function LoadvectorfieldsButtonPushed(app, event)
            tmp = uigetdir(path);
            if tmp == 0;
    	       return;
            end
            app.Path=tmp;
            app.vec_flag = 1; 
            app.work_post=[];
            set(findall(app.PTVResearchUIFigure, '-property', 'Enable'), 'Enable', 'on')
            app.Lagvariable.Visible='off';
            app.TabGroup.SelectedTab=app.TabGroup.Children(3);
            TabGroupSelectionChanged(app, event)
        end

        % Value changed function: ofimagesformeanEditField
        function PreviewButtonPushed(app, event)
            if size(read(app.vid,round(app.Slider.Value)),3)==3
               data = double(app.mask(:,:,1).*rgb2gray(read(app.vid,round(app.Slider.Value))));
            else
               data = double(app.mask(:,:,1).*read(app.vid,round(app.Slider.Value)));
            end
            
            if app.BackgroundsubtractionButton.Value==1;
                data=data-app.M;
          
            end
            
            if app.SharpenButton.Value==1;
                data=imsharpen(data,'Amount',app.amountEditField.Value/100);
            end
            
            imagesc(app.UIAxes,single(app.mask(:,:,1)).*data);
        end

        % Callback function
        function ofimagesEditFieldValueChanged(app, event)
            BackgroundsubtractionButtonValueChanged(app)
            %             MeanFilterButtonValueChanged
        end

        % Button pushed function: PODDEMFilterButton
        function PODDEMFilterButtonPushed(app, event)
            U = single(app.Uz);
            V = single(app.Vz);
            mask = single(isnan(U));
            mask(mask==1)=nan;            
            U(isnan(U))=0;
            V(isnan(V))=0;
            [u s a]=svd([reshape(U,size(U,1)*size(U,2),size(U,3));...
                    reshape(V,size(V,1)*size(V,2),size(V,3))],'econ');
            a=double(a);
            for i = 1:length(a);
                a(:,i)=medfilt1(a(:,i),10);
            end
            recon = u*s*a';
            app.Uz = reshape(recon(1:size(U,1)*size(U,2),:),size(U,1),size(U,2),size(U,3));
            app.Vz = reshape(recon(size(U,1)*size(U,2)+1:end,:),size(U,1),size(U,2),size(U,3));
            app.Uz = app.Uz+mask;
            app.Vz = app.Vz+mask;
            SliderValueChanged(app);
            app.Postvariable_mean.Items={''}
            names = dir([app.Path,'/','*.mat']);
            mkdir([app.Path,'/','FILTERED'])
            for i = 1:size(cat(1,names.name),1) 
                load([app.Path,'/',names(i).name]);
                u=app.Uz(:,:,i);
                v=app.Vz(:,:,i);
                if exist('P1');
                    save([app.Path,'/','FILTERED/',names(i).name],'u','v','P1','P2');
                else
                    save([app.Path,'/','FILTERED/',names(i).name],'u','v');
                end
            end
        end

        % Value changed function: Postvariable_mean
        function Postvariable_meanValueChanged(app, event)
            DisplayButton_2Pushed(app)
        end

        % Button pushed function: DisplayButton_2
        function DisplayButton_2Pushed(app, event)
            if strcmp(app.Postvariable_mean.Value,'POD u') |...
                    strcmp(app.Postvariable_mean.Value,'POD v') |...
                    strcmp(app.Postvariable_mean.Value,'Coef vel') |...
                    strcmp(app.Postvariable_mean.Value,'Spectra vel');
                    app.Slider_pod.Visible='on'
                    app.Slider.Visible='off'
                    
                    num = round(app.Slider_pod.Value)
                    if strcmp(app.Postvariable_mean.Value,'POD u')
                        imagesc(app.UIAxes,app.P_u(:,:,num));
                        colorbar(app.UIAxes)
                    end
                    if strcmp(app.Postvariable_mean.Value,'POD v')
                        imagesc(app.UIAxes,app.P_v(:,:,num));
                        colorbar(app.UIAxes)
                    end
                    if strcmp(app.Postvariable_mean.Value,'Coef vel')
                        plot(app.UIAxes,app.A(:,num),'k','linewidth',1);
                        colorbar(app.UIAxes,'off')
                        axis(app.UIAxes,'tight')
                        
                    end
                    if strcmp(app.Postvariable_mean.Value,'Spectra vel')
                        plot(app.UIAxes,app.S,'k-o','linewidth',1);
                        colorbar(app.UIAxes,'off')
                        axis(app.UIAxes,'tight')
                    end

            else
                    app.Slider_pod.Visible='off'
                    app.Slider.Visible='on'
                
                if strcmp(app.Postvariable_mean.Value,'U')
                    imagesc(app.UIAxes,app.MEAN.U_mean);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'V')
                    imagesc(app.UIAxes,app.MEAN.V_mean);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'Vector Length')
                    imagesc(app.UIAxes,app.MEAN.VEC_mean);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'Vorticity')
                    imagesc(app.UIAxes,app.MEAN.Vorticity);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'Quiver')
                    quiver(app.UIAxes,app.MEAN.U_mean,app.MEAN.V_mean,10,'color','k');
    %                 colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'uu')
                    imagesc(app.UIAxes,app.RANS.UU);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'uv')
                    imagesc(app.UIAxes,app.RANS.UV);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable_mean.Value,'vv')
                    imagesc(app.UIAxes,app.RANS.VV);
                    colorbar(app.UIAxes)
                end            
                if strcmp(app.Postvariable_mean.Value,'tke')
                    imagesc(app.UIAxes,app.RANS.TKE);
                    colorbar(app.UIAxes)
                end                 
                if strcmp(app.Postvariable_mean.Value,'Fluc Vorticity')
                    imagesc(app.UIAxes,app.VAR.VORT);
                    colorbar(app.UIAxes)
                end                 
                if strcmp(app.Postvariable_mean.Value,'Swirling Strength')
                    imagesc(app.UIAxes,app.VAR.SWIRL);
                    colorbar(app.UIAxes)
                end      
            end
            caxis(app.UIAxes,'auto');
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            app.Uz = app.Orig.U;
            app.Vz = app.Orig.V;
            
            app.Postvariable_mean.Items={''}
            app.Slider_pod.Visible='off';
            app.Slider.Visible='on';
            SliderValueChanged(app);
        end

        % Callback function
        function MedianFilterButton_2ValueChanged(app, event)

        end

        % Button pushed function: SmoothingButton
        function SmoothingButtonPushed(app, event)
            U = single(app.Uz);
            V = single(app.Vz);
            app.PROLabel.FontColor='k';
            app.PROLabel.FontSize=28;
            names = dir([app.Path,'/','*.mat']);
            mkdir([app.Path,'/','FILTERED'])
            for i = 1:size(cat(1,names.name),1);
                U(:,:,i)=medfilt2(U(:,:,i),[app.KernelEditField.Value,app.KernelEditField.Value]);
                V(:,:,i)=medfilt2(V(:,:,i),[app.KernelEditField.Value,app.KernelEditField.Value]);
                app.PROLabel.Text=sprintf('%0.0f %%',100*(i/size(U,3)));
                load([app.Path,'/',names(i).name]);
                u=U(:,:,i);
                v=V(:,:,i);
                if exist('P1');
                    save([app.Path,'/','FILTERED/',names(i).name],'u','v','P1','P2');
                else
                    save([app.Path,'/','FILTERED/',names(i).name],'u','v');
                end
            end
            app.Uz=U;
            app.Vz=V;
            app.PROLabel.FontColor=[0.94 0.94 0.94];
            SliderValueChanged(app);
            app.Postvariable_mean.Items={''}
        end

        % Value changed function: Postvariable
        function PostvariableValueChanged(app, event)
            value = app.Postvariable.Value;
            SliderValueChanged(app);
        end

        % Button pushed function: PODButton
        function PODButtonPushed(app, event)
            answer = questdlg('You are about to run a POD, this is computationally heavy and needs evenly spaces data!',...
                'POD','Continue','Cancel','')
            if strcmp(answer,'Continue');
                U = single(app.Uz);
                V = single(app.Vz);
                U(isnan(U))=0;
                V(isnan(V))=0;
                U=U-repmat(app.Um,[1 1 size(U,3)]);
                V=V-repmat(app.Vm,[1 1 size(U,3)]);
                U(isnan(U))=0;
                V(isnan(V))=0;
                [P S A]=svd([reshape(U,size(U,1)*size(U,2),size(U,3));...
                        reshape(V,size(V,1)*size(V,2),size(V,3))],'econ');
            end
            mx = size(U,3);
            app.P_u=reshape(P(1:size(U,1)*size(U,2),1:min([mx,100])),size(U,1),size(U,2),min([mx,100]));
            app.P_v=reshape(P(size(U,1)*size(U,2)+1:end,1:min([mx,100])),size(U,1),size(U,2),min([mx,100]));
            app.S=diag(S(1:min([mx,100]),1:min([mx,100])));
            app.A=A(:,1:min([mx,100]));
            imagesc(app.UIAxes,app.P_u(:,:,1))
            fn=sprintf('%s/POD',app.Path);  
            mkdir(fn)
            Pu = app.P_u;   
            Pv = app.P_v;        
            A = app.A(:,1:min([mx,100]));  
            S=diag(S(1:min([mx,100]),1:min([mx,100])));
            save(sprintf('%s/POD.mat',fn),'Pu','Pv','A','S');
            
            items={''};
            if ~isempty(app.MEAN);
                add = [{'U','V','Vorticity','Vector Length','Quiver'}]
                items = cat(2,items,add);
            end
            if ~isempty(app.RANS);
                add = [{'uu','uv','vv','tke'}]
                items = cat(2,items,add);
            end
            if ~isempty(app.VAR);
                add = [{'Fluc Vorticity','Swirling Strength'}]
                items = cat(2,items,add);
            end
            if ~isempty(app.P_u);
                add = [{'POD u','POD v','Coef vel','Spectra vel'}]
                items = cat(2,items,add);
            end                
            items(1)=[];
            app.Postvariable_mean.Items=items;
            app.Postvariable_mean.Value='POD u'
            app.Slider_pod.Visible='on';
            app.Slider_pod.Limits=[1 size(app.A,2)]
            app.Slider_pod.Value=1;
            app.Slider_pod.Visible='on';
            app.Slider.Visible='off';
%             app.CurrentFrameEditFieldLabel='POD Modes';
        end

        % Value changing function: Slider_pod
        function Slider_podValueChanging(app, event)
            changingValue = event.Value;
        end

        % Value changed function: Slider_pod
        function Slider_podValueChanged(app, event)
            num = round(app.Slider_pod.Value);
%             app.CurrentFrameEditField.Value=round(num);

            if strcmp(app.Postvariable_mean.Value,'POD u')
                imagesc(app.UIAxes,app.P_u(:,:,num));
                colorbar(app.UIAxes)
            end
            if strcmp(app.Postvariable_mean.Value,'POD v')
                imagesc(app.UIAxes,app.P_v(:,:,num));
                colorbar(app.UIAxes)
            end
            if strcmp(app.Postvariable_mean.Value,'Coef vel')
                plot(app.UIAxes,app.A(:,num),'k','linewidth',1);
                colorbar(app.UIAxes,'off')
                axis(app.UIAxes,'tight')
                
            end
            if strcmp(app.Postvariable_mean.Value,'Spectra vel')
                plot(app.UIAxes,app.S,'k-o','linewidth',1);
                colorbar(app.UIAxes,'off')
                axis(app.UIAxes,'tight')
            end
            app.CurrentFrameEditField.Value=round(num);
        end

        % Value changed function: DropDown
        function DropDownValueChanged(app, event)
            SliderValueChanged(app, event)
        end

        % Value changed function: minEditField
        function minEditFieldValueChanged(app, event)
            SliderValueChanged(app, event)
        end

        % Value changed function: maxEditField
        function maxEditFieldValueChanged(app, event)
            SliderValueChanged(app, event)
        end

        % Value changed function: StartFrameEditField
        function StartFrameEditFieldValueChanged(app, event)
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -')
                tmp = size(app.Uz,3) + app.Slider.Limits(1);
                if app.StartFrameEditField.Value<app.Slider.Limits(1);
                    app.StartFrameEditField.Value = app.Slider.Limits(1);
                end
                app.EndFrameEditField.Value=min([app.Slider.Limits(2),app.Slider.Limits(1)+app.NumbertoprocessEditField.Value-1])
            else      
                value = app.StartFrameEditField.Value;
                app.CurrentFrameEditField.Value=value;
%                 app.Slider.Value=value;   
                app.EndFrameEditField.Value=min([round(value+app.NumbertoprocessEditField.Value)-1,round(app.vid.Duration.*app.vid.FrameRate)]);
                imagesc(app.UIAxes,app.mask.*read(app.vid,round(app.Slider.Value)))
            end
            
        end

        % Value changed function: NumbertoprocessEditField
        function NumbertoprocessEditFieldValueChanged(app, event)
%             value = app.NumbertoprocessEditField.Value;
            if ~isempty(app.vid)
                if app.NumbertoprocessEditField.Value+app.StartFrameEditField.Value-1 > app.vid.Duration.*app.vid.FrameRate;
                   app.EndFrameEditField.Value=app.vid.Duration.*app.vid.FrameRate;
                else 
                    app.EndFrameEditField.Value=app.NumbertoprocessEditField.Value+app.StartFrameEditField.Value-1;
                end
            else
                app.EndFrameEditField.Value=round(min([app.Slider.Limits(2),app.Slider.Limits(1)+app.NumbertoprocessEditField.Value]));

%                 app.NumbertoprocessEditField.Value = 1;
            end

        end

        % Selection change function: TabGroup
        function TabGroupSelectionChanged(app, event)
            if app.vec_flag == 1
                if strcmp(app.TabGroup.SelectedTab.Title,'Process');
                    app.TabGroup.SelectedTab=app.TabGroup.Children(3);
                end
                if strcmp(app.TabGroup.SelectedTab.Title,'Pre -');
                    app.TabGroup.SelectedTab=app.TabGroup.Children(3);
                end
            end
            
            if isempty(dir([app.Path,'/*.mat']))
                if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -');
                    app.TabGroup.SelectedTab=app.TabGroup.Children(2);
                end
                if strcmp(app.TabGroup.SelectedTab.Title,'Post Lag -');
                    app.TabGroup.SelectedTab=app.TabGroup.Children(2);
                end
            else
                if strcmp(app.TabGroup.SelectedTab.Title,'Post Lag -')
                    tmp = dir([app.Path,'/','*.mat']); 
                    load([app.Path,'/',tmp(1).name])
                    if ~exist('P1');
                        app.TabGroup.SelectedTab=app.TabGroup.Children(3);
                    end
                end
            end
            
            if app.run_flag_pro==1;
                if strcmp(app.prev_pro,'Process');
                    app.TabGroup.SelectedTab=app.TabGroup.Children(2);
                end
                if strcmp(app.prev_pro,'Pre -');
                    app.TabGroup.SelectedTab=app.TabGroup.Children(1);
                end
            end
            if app.run_flag_post==1;
                app.TabGroup.SelectedTab=app.TabGroup.Children(3);
            end
            
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -')
                app.Postvariable.Visible='on';
                app.DisplayVariableDropDown.Visible='off';
                app.Lagvariable.Visible='off';
            end
                       
            if strcmp(app.TabGroup.SelectedTab.Title,'Pre -') | strcmp(app.TabGroup.SelectedTab.Title,'Process')
                app.Postvariable.Visible='off';
                app.DisplayVariableDropDown.Visible='on';
                app.Lagvariable.Visible='off';
            end
            
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Lag -')
                app.Slider.Enable='on';  
                app.Postvariable.Visible='off';
                app.Lagvariable.Visible='on';
                tmp = dir([app.Path,'/','*.mat']); 
                start_num = str2num(tmp(1).name(8:end-4));
                end_num = str2num(tmp(end).name(8:end-4));
                tmp = dir([app.Path,'/','*.mat']); 
                app.PROLabel.FontColor='k';
                count = 0
                for i = start_num:end_num;
                    count = count + 1
                    load([app.Path,'/',tmp(count).name]);
                    app.ul{i}=P2(1:size(P1,1),1)-P1(:,1);
                    app.vl{i}=P2(1:size(P1,1),2)-P1(:,2);
                    app.xl{i}=P1(:,1);
                    app.yl{i}=P1(:,2);
                    app.PROLabel.Text=sprintf('%0.0f %%',100*(i/size(cat(1,tmp.name),1)));
                end
                app.lag_lims=[start_num,end_num];
                app.Slider.Limits=[start_num,end_num];
                app.PROLabel.FontColor=[0.94 0.94 0.94];
            end
            
            
            if strcmp(app.TabGroup.SelectedTab.Title,'Post Grid -')
                app.MeanFieldsButton.Value=0;
                app.SwirlVorticityButton.Value=0;
                app.RANSQuantitiesButton.Value=0;
                app.DisplayVariableDropDown.Visible='off';
                app.Postvariable_mean.Items={''};
                
                app.Slider.Enable='on';    
                app.DisplayVariableDropDown.Enable='on';
                app.post=1; % this needs to be removed
                app.PROLabel.FontColor='k';
                app.PROLabel.FontSize=28;
                tmp = dir([app.Path,'/','*.mat']); 
                for i = 1:size(cat(1,tmp.name),1);
                    load([app.Path,'/',tmp(i).name]);
                    app.Uz(:,:,i)=single(u);
                    app.Vz(:,:,i)=single(v);
                    app.PROLabel.Text=sprintf('%0.0f %%',100*(i/size(cat(1,tmp.name),1)));
                end
                app.Orig.U = app.Uz;
                app.Orig.V = app.Vz;
                app.PROLabel.FontColor=[0.94 0.94 0.94];
                app.RUNButton.Enable='on';
                app.STOPButton.Enable='on';
                app.UIAxes.XLim=[1 size(app.Uz,2)];
                app.UIAxes.YLim=[1 size(app.Uz,1)];
                start_num = str2num(tmp(1).name(8:end-4));
                end_num = str2num(tmp(end).name(8:end-4));
                app.pro_lims=[start_num,end_num];

                if isempty(app.work_post)                    
                    app.work_post=start_num;
                    app.Slider.Value=app.work_post;
                end

                app.Slider.Limits=app.pro_lims
                app.CurrentFrameEditField.Value=app.work_post;
                app.StartFrameEditField.Value=app.work_post;
                app.EndFrameEditField.Value=end_num;
                app.NumbertoprocessEditField.Value=end_num - start_num -1;
                
                
                U=app.Uz(:,:,ceil(app.work_post)-app.pro_lims(1)+1);
                V=app.Vz(:,:,ceil(app.work_post)-app.pro_lims(1)+1); 
                
                app.Um = mean(app.Uz,3);
                app.Vm = mean(app.Vz,3);
                u = U-app.Um;
                v = V-app.Vm;
                if strcmp(app.Postvariable.Value,'U')
                    imagesc(app.UIAxes,U);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'V')
                    imagesc(app.UIAxes,V);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Vorticity');
                    imagesc(app.UIAxes,curl(U,V));
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Vector Length');
                    imagesc(app.UIAxes,sqrt(U.^2+V.^2)/2);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'uu');
                    imagesc(app.UIAxes,u.*u);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'uv');
                    imagesc(app.UIAxes,u.*v);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'vv');
                    imagesc(app.UIAxes,v.*v);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'TKE');
                    imagesc(app.UIAxes,(v.^2.+v.^2)/2);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Fluc Vorticity');
                    imagesc(app.UIAxes,curl(u,v));
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Swirl');
                    r = size(u,1);
                    c= size(u,2);
                    [du_dx,du_dy]=gradient(u);
                    [dv_dx,dv_dy]=gradient(v);                
                    Q=(du_dy.*dv_dx).*(-0.5); % Coherent structures identifier done using Q-criterion
                    lam=zeros(r,c);
                    for w=1:r;
                        for j=1:c;
                            D=[du_dx(w,j) du_dy(w,j);dv_dx(w,j) dv_dy(w,j)];
                            if isnan(det(D))~=1;
                               lamtemp=imag(eig(D));
                               pos=find(lamtemp>0);
                                if length(pos)~=0;
                                lam(w,j)=lamtemp(pos);
                                end
                            end
                        end
                    end
                    imagesc(app.UIAxes,lam);
                    colorbar(app.UIAxes)
                end
                if strcmp(app.Postvariable.Value,'Quiver')
                    quiver(app.UIAxes,U,V,2,'k');
                    colorbar(app.UIAxes,'off')
                end
                colormap(app.UIAxes,'jet')    
                caxis(app.UIAxes,'auto');
            else   
                app.Slider.Limits=[1,app.vid.Duration.*app.vid.FrameRate];
                app.Slider.Value=app.work_pre;
            end            
            SliderValueChanged(app, event)
            
            items={''};
            if ~isempty(app.MEAN);
                add = [{'U','V','Vorticity','Vector Length','Quiver'}]
                items = cat(2,items,add);
            end
            if ~isempty(app.RANS);
                add = [{'uu','uv','vv','tke'}]
                items = cat(2,items,add);
            end
            if ~isempty(app.VAR);
                add = [{'Fluc Vorticity','Swirling Strength'}]
                items = cat(2,items,add);
            end
            if ~isempty(app.P_u);
                add = [{'POD u','POD v','Coef vel','Spectra vel'}]
                items = cat(2,items,add);
            end       
            items(1)=[];
            app.Postvariable_mean.Items=items;
        end

        % Value changed function: CalculateLagrangianTracksButton
        function CalculateLagrangianTracksButtonValueChanged(app, event)
            if app.CalculateLagrangianTracksButton.Value==1
                app.DisplayVariableDropDown.Items=[{'Points','Vector Length','U','V','Quiver'}]
                app.DisplayVariableDropDown.Value='Points'
            else
                app.DisplayVariableDropDown.Items=[{'Vector Length','U','V','Quiver'}]
            end
        end

        % Button pushed function: MaskButton
        function MaskButtonPushed(app, event)
            waitfor(msgbox('Select the desired points once chosen, press o to mask outside, i to mask inside or c to cancel. You can keep adding to the mask as many times as you need.'))
            figure(1) 
            imagesc(uint8(app.mask).*read(app.vid,round(app.CurrentFrameEditField.Value)))
            set(gca,'YDir','normal')
            count=0;
            hold on
            X=[];Y=[];
            flag = 0; 
            while flag == 0;
                if strcmp(get(gcf,'CurrentKey'),'i'); flag = 1; end
                if strcmp(get(gcf,'CurrentKey'),'o'); flag = 1; end
                if strcmp(get(gcf,'CurrentKey'),'c'); flag = 1; end
                count = count +1; 
                [xs ys]=ginput(1);
                X(count)=xs; 
                Y(count)=ys;
                plot(X(count),Y(count),'or')
                drawnow
                if strcmp(get(gcf,'CurrentKey'),'i'); flag = 1; end
                if strcmp(get(gcf,'CurrentKey'),'o'); flag = 1; end
                if strcmp(get(gcf,'CurrentKey'),'c'); flag = 1; end
            end
            X(end+1)=X(1);Y(end+1)=Y(1);
            [xy,yv]=meshgrid(1:size(read(app.vid,round(app.CurrentFrameEditField.Value)),2),1:size(read(app.vid,round(app.CurrentFrameEditField.Value)),1));
            
            if strcmp(get(gcf,'CurrentKey'),'o')
                mask = inpolygon(xy,yv,X(:),Y(:));
                mask = repmat(uint8(mask),[1 1 size(read(app.vid,round(app.CurrentFrameEditField.Value)),3)]);
                app.mask = app.mask.*mask;
                mkdir([app.Path,'/','PARAMS']);
                save([app.Path,'/','PARAMS/','mask.mat'],'mask');
            end
            if strcmp(get(gcf,'CurrentKey'),'i')
                mask = 1-inpolygon(xy,yv,X(:),Y(:));
                mask = repmat(uint8(mask),[1 1 size(read(app.vid,round(app.CurrentFrameEditField.Value)),3)]);
                app.mask = app.mask.*mask;
                mkdir([app.Path,'/','PARAMS']);
                save([app.Path,'/','PARAMS/','mask.mat'],'mask');
            end
            if strcmp(get(gcf,'CurrentKey'),'c')
                mask = app.mask;
                close(figure(1));
            end
            close(figure(1))
            imagesc(app.UIAxes,...
                app.mask.*read(app.vid,round(app.CurrentFrameEditField.Value)))
            app.UIAxes.YDir='normal';
           
        end

        % Button pushed function: PlayButton
        function PlayButtonPushed(app, event)
             while app.Slider.Value < app.Slider.Limits(2)-1;   
                if app.pause == 1;
                    break
                end
                app.Slider.Value=app.Slider.Value + 1; 
                SliderValueChanged(app, event)
                speeds = fliplr([0.01 0.1 0.2 0.25 0.3 0.5 0.75 1 1.5 2]);
                pause(speeds(ceil(app.speed)))
            end
            app.pause = 0;
        end

        % Button pushed function: PauseButton
        function PauseButtonPushed(app, event)
            app.pause = 1;
        end

        % Value changed function: PlaybackSpeedSlider
        function PlaybackSpeedSliderValueChanged(app, event)
            value = app.PlaybackSpeedSlider.Value;
            app.speed=value;
         
        end

        % Button pushed function: Forward
        function ForwardButtonPushed(app, event)
            if app.Slider.Value < app.Slider.Limits(2)-1
                app.Slider.Value=app.Slider.Value + 1; 
                SliderValueChanged(app, event)
            end
        end

        % Button pushed function: Backward
        function BackwardButtonPushed(app, event)
            if app.Slider.Value > app.Slider.Limits(1)+1;
                app.Slider.Value=app.Slider.Value - 1;
                SliderValueChanged(app, event)
            end
        end

        % Callback function
        function CalculateLagrangianPathsButtonPushed(app, event)
            
        
        end

        % Button pushed function: DisplayButton
        function DisplayButtonPushed(app, event)
            cla(app.UIAxes);

            if strcmp(app.PostvariableLag.Value,'Diff-x') | strcmp(app.PostvariableLag.Value,'Diff-y');
                cla(app.UIAxes);
                app.UIAxes.YDir='normal';
                load([app.Path,'/','LANG/','STATS/DIFF.mat'])
                if strcmp(app.PostvariableLag.Value,'Diff-x')
                    plot(app.UIAxes,1:length(xdiff),xdiff)
                end
                if strcmp(app.PostvariableLag.Value,'Diff-y')
                    plot(app.UIAxes,1:length(ydiff),ydiff)
                end
                axis(app.UIAxes,'tight');
            end

            if strcmp(app.PostvariableLag.Value,'Tracks Vec Length') | strcmp(app.PostvariableLag.Value,'Tracks u-vel') | strcmp(app.PostvariableLag.Value,'Tracks v-vel')
                names = dir([app.Path,'/','LANG','/','*.mat'])
                load([app.Path,'/','LANG','/',names(1).name]);
                tmps = round(linspace(1,size(tracks,1),size(tracks,1).*(app.trackstodisplayEditField.Value/100)));
                size(tmps)
                colormap(app.UIAxes,'jet')
                for i = 2:size(cat(1,names.name),1);
                   load([app.Path,'/','LANG','/',names(i-1).name]);
                   t1 = tracks;
                   load([app.Path,'/','LANG','/',names(i).name]);
                   t2 = tracks; 
                   vel = (abs(t2(1:size(t1,1),1)-t1(1:size(t1,1),1)) + abs(t2(1:size(t1,1),2)-t1(1:size(t1,1),2)))./2;     
                   velu = t2(1:size(t1,1),1)-t1(1:size(t1,1),1);   
                   velv = t2(1:size(t1,1),2)-t1(1:size(t1,1),2);   
                   if strcmp(app.PostvariableLag.Value,'Tracks Vec Length');
                       scatter(app.UIAxes,t1(tmps,1),t1(tmps,2),2,vel(tmps))
                        hold(app.UIAxes,'on')
%                        caxis(app.UIAxes,([1.2*min(vel(~isnan(vel))),0.8*max(vel(~isnan(vel)))]));
                   end
                   if strcmp(app.PostvariableLag.Value,'Tracks u-vel');
                       scatter(app.UIAxes,t1(tmps,1),t1(tmps,2),2,velu(tmps))
                       hold(app.UIAxes,'on')
                   end
                   if strcmp(app.PostvariableLag.Value,'Tracks v-vel');
                       scatter(app.UIAxes,t1(tmps,1),t1(tmps,2),2,velv(tmps))
                       hold(app.UIAxes,'on')
                   end
                end
                colorbar(app.UIAxes);
                if strcmp(app.DropDown.Value,'auto')
                    tmp=double(app.UIAxes.CLim);
                    caxis(app.UIAxes,'auto');
                    app.minEditField.Value=tmp(1);
                    app.maxEditField.Value=tmp(2);
                end
                if strcmp(app.DropDown.Value,'manual')
                    caxis(app.UIAxes,[app.minEditField.Value,app.maxEditField.Value]);
                end
            end
          
        end

        % Value changed function: SharpenButton
        function SharpenButtonValueChanged(app, event)
            SliderValueChanged(app, event)            
        end

        % Value changed function: WienerFilterButton
        function WienerFilterButtonValueChanged(app, event)
            SliderValueChanged(app, event)
        end

        % Callback function
        function LocalIntensityButtonPushed(app, event)
            SliderValueChanged(app, event)
        end

        % Value changed function: LocalIntensityButton
        function LocalIntensityButtonValueChanged(app, event)
            value = app.LocalIntensityButton.Value;
            SliderValueChanged(app, event)
        end

        % Callback function
        function LoadMaskButtonPushed(app, event)
            [file,path] = uigetfile([{'*.mat'}])
            load([path,file])
            if size(app.mask)==mask;
                app.mask = mask; 
            else
                msgbox('Size of mask does not match image');
            end
            
        end

        % Button pushed function: LoadparamsfrompreviousButton
        function LoadparamsfrompreviousButtonPushed(app, event)
            [file,path] = uigetfile([{'*.txt'}]);
            FID = fopen([path,file]);
            p = textscan(FID,'%s')
            for i = 1:size(p{1},1);
                p{1}(i);
                eval(p{1}{i});
            end
            fclose(FID);
            try
                load([path,'tform.mat']);
                app.tform=tform;
            end
 
        end

        % Button pushed function: ScalefromreferenceButton
        function ScalefromreferenceButtonPushed(app, event)

            answer = questdlg('Load from external image / video?','Select two points','Yes','No','No');
            waitfor(answer)
            if strcmp(answer,'No')
                figure(1) 
                correct = 0
                while correct == 0
                    cla
                    imagesc(uint8(app.mask).*read(app.vid,round(app.CurrentFrameEditField.Value)))
                    set(gca,'YDir','normal')
                    hold on
                    X=[];Y=[];
                    for i  = 1:2 
                        [x y]=ginput(1);
                        X(i)=x;
                        Y(i)=y;
                        plot(X,Y,'-ro','linewidth',2);
                    end
                    scale = sqrt(diff(X).^2 + diff(Y).^2); 
                    a=inputdlg('Distance in mm');
                    text(mean(X),mean(Y),[a{1},'mm'],'fontsize',25,'color','r');
                    answer = questdlg('Are these OK?',...
                        'Scale','Yes','No','')
                    waitfor(answer);
                    if strcmp(answer,'Yes');
                        correct = 1;
                        close(figure(1))
                    end
                    app.PxmmEditField.Value = scale/str2num(a{1});
                end
                close(figure(1));
            end
            if strcmp(answer,'Yes');
                [file,path]=uigetfile('*.*');
                try
                    v = VideoReader([path,file]);
                    a=(read(v,[1]));
                catch
                    a=imread([path,file]);
                end
                imagesc(a);
                hold on
                set(gca,'YDir','normal')
                X=[];Y=[];
                for i  = 1:2 
                    [x y]=ginput(1);
                    X(i)=x;
                    Y(i)=y;
                    plot(X,Y,'-ro','linewidth',2);
                end
                scale = sqrt(diff(X).^2 + diff(Y).^2); 
                a=inputdlg('Distance in mm');
                text(mean(X),mean(Y),[a{1},'mm'],'fontsize',25,'color','r');
                answer = questdlg('Are these OK?',...
                    'Scale','Yes','No','')
                waitfor(answer);
                if strcmp(answer,'Yes');
                    correct = 1;
                    close(figure(1))
                end
                app.PxmmEditField.Value = scale/str2num(a{1});
                close(figure(1));
            end
            
                            
        end

        % Button pushed function: ResetMaskButton
        function ResetMaskButtonPushed(app, event)
            app.mask = uint8(ones(size(app.mask)));
            SliderValueChanged(app, event);
        end

        % Value changed function: Lagvariable
        function LagvariableValueChanged(app, event)
            value = app.Lagvariable.Value;
            SliderValueChanged(app, event)
        end

        % Button pushed function: LoadrawimagesButton
        function LoadrawimagesButtonPushed(app, event)
            [file,path] = uigetfile('MultiSelect','on','*.*');
%             if file==0;
%     	       return;
%             end
            vid = VideoWriter([path,'files.avi']);
            open(vid)
            app.PROLabel.FontColor='k';
            app.PROLabel.FontSize=28;
            for i = 1:size(file,2)
                tmp = imread([path,file{i}]);
                writeVideo(vid,uint8(tmp))
                app.PROLabel.Text=sprintf('%0.0f %%',100*(i/size(file,2)));

            end
            close(vid)
            app.PROLabel.FontColor=[0.94 0.94 0.94];
            app.images_raw.path = path;
            app.images_raw.file = file;
            LoadrawvideoButtonPushed(app, event)
        end

        % Button pushed function: LoadMaskButton
        function LoadMaskButtonPushed2(app, event)
            [file,path] = uigetfile([{'*.mat'}])
            load([path,file])
            if size(app.mask)==size(mask);
                app.mask = mask; 
            else
                msgbox('Size of mask does not match image');
            end
            SliderValueChanged(app, event)
        end

        % Value changed function: SnapshotPIVButton
        function SnapshotPIVButtonValueChanged(app, event)
            if app.SnapshotPIVButton.Value==1
                app.LagrangianPathsButton.Value=0;
                app.PTVButton.Value=0; 
                app.CalculateLagrangianTracksButton.Enable='off'
            end
            if app.SnapshotPIVButton.Value==0 & app.PTVButton.Value==0
                
                app.SnapshotPIVButton.Value=1
            end
        end

        % Value changed function: PTVButton
        function PTVButtonValueChanged(app, event)
            if app.PTVButton.Value==1
                app.SnapshotPIVButton.Value=0;
                app.CalculateLagrangianTracksButton.Enable='on'
            end
            if app.SnapshotPIVButton.Value==0 & app.PTVButton.Value==0
                app.PTVButton.Value=1
            end
        end

        % Value changed function: ImageDewarpingButton
        function ImageDewarpingButtonValueChanged(app, event)
            if app.ImageDewarpingButton.Value==1; 
                answer3='';
                if ~isempty(app.tform);
                    answer3 = questdlg('Use Previous Points','','Yes','No','Yes');
                    waitfor(answer3)
                end
                if isempty(app.tform) | strcmp(answer3,'No') ;
                    % use previous points line here 
                    [file,path]=uigetfile('*.*');
                    if file == 0;
                        app.ImageDewarpingButton.Value=0;
                        return;
                    end
                    try
                        v = VideoReader([path,file]);
                        img=(read(v,[1]));
                    catch
                        img=imread([path,file]);
                    end
                    answer2 = 'No';
                    waitfor(msgbox('Select all point in each row from left to right, starting in the bottom left hand corner. Press ''e'' when complete'));
                    while strcmp(answer2, 'No')
                        close(figure(1));
                        figure(1)
                        X=[];Y=[];
                        imagesc(img)
                        set(gca,'YDir','normal')
                        hold on
                        flag = 0;
                        i = 0
                        while flag == 0;
                            if strcmp(get(gcf,'CurrentKey'),'e'); break;end
                            i = i+1; 
                            [x y]=ginput(1);
                            X(i)=x;
                            Y(i)=y;
                            plot(X,Y,'ro')
                        end
                        answer2 = questdlg('Are these points correct?','','No','Yes','Yes');
                        waitfor(answer2)
                    end
                    X(end)=[];Y(end)=[];
                    a=inputdlg('How many columns in the grid?');
                    waitfor(a);
                    c = str2num(a{1});
                    r = length(X)/c;
                    xs = linspace(X(1),X(end),c);
                    [x y]=meshgrid(xs,linspace(Y(1),Y(1)+r*mean(diff(xs)),r)); 
                    app.warp.x=x; 
                    app.warp.y=y; 
                    app.warp.Y=Y;
                    app.warp.X=X;
                    close(figure(1))
                    MethodDropDownValueChanged(app, event)
                    SliderValueChanged(app, event)
                    app.PxmmEditField.Value = (mean(diff(xs))./app.GridspacingmmEditField.Value);
                end
            end
            SliderValueChanged(app, event);
        end

        % Value changed function: MethodDropDown
        function MethodDropDownValueChanged(app, event)
            if strcmp(app.MethodDropDown.Value,'Poly2');
                try
                    tform = fitgeotrans([app.warp.X(:),app.warp.Y(:)],[app.warp.y(:),app.warp.x(:)],'poly',2)
                catch msgbox('Not enough grid points');
                    app.ImageDewarpingButton.Value=0;
                    return;
                end
            end
            if strcmp(app.MethodDropDown.Value,'Poly3');
                try
                    tform = fitgeotrans([app.warp.X(:),app.warp.Y(:)],[app.warp.y(:),app.warp.x(:)],'poly',3);
                catch msgbox('Not enough grid points');
                    app.ImageDewarpingButton.Value=0;
                    return;
                end
            end
            if strcmp(app.MethodDropDown.Value,'PWL');
                try
                    tform = fitgeotrans([app.warp.X(:),app.warp.Y(:)],[app.warp.y(:),app.warp.x(:)],'pwl');
                    catch msgbox('Not enough grid points');
                    app.ImageDewarpingButton.Value=0;
                    return;
                end
            end
            app.tform = tform;
            mkdir([app.Path,'/','PARAMS']);
            save([app.Path,'/','PARAMS/','tform.mat'],'tform');
            SliderValueChanged(app, event)
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PTVResearchUIFigure
            app.PTVResearchUIFigure = uifigure;
            app.PTVResearchUIFigure.Position = [100 100 915 668];
            app.PTVResearchUIFigure.Name = 'PTVResearch';

            % Create LoadrawvideoButton
            app.LoadrawvideoButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.LoadrawvideoButton.ButtonPushedFcn = createCallbackFcn(app, @LoadrawvideoButtonPushed, true);
            app.LoadrawvideoButton.Position = [528 628 94 22];
            app.LoadrawvideoButton.Text = 'Load raw video';

            % Create RUNButton
            app.RUNButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.RUNButton.ButtonPushedFcn = createCallbackFcn(app, @RUNButtonPushed, true);
            app.RUNButton.BackgroundColor = [0 1 0];
            app.RUNButton.Enable = 'off';
            app.RUNButton.Position = [767 30 100 53];
            app.RUNButton.Text = 'RUN';

            % Create TestButton
            app.TestButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.TestButton.ButtonPushedFcn = createCallbackFcn(app, @TestButtonPushed, true);
            app.TestButton.Enable = 'off';
            app.TestButton.Position = [646 58 100 25];
            app.TestButton.Text = 'Test';

            % Create Slider
            app.Slider = uislider(app.PTVResearchUIFigure);
            app.Slider.MajorTickLabels = {'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Slider.Enable = 'off';
            app.Slider.Position = [57 155 401 3];

            % Create CurrentFrameEditFieldLabel
            app.CurrentFrameEditFieldLabel = uilabel(app.PTVResearchUIFigure);
            app.CurrentFrameEditFieldLabel.HorizontalAlignment = 'right';
            app.CurrentFrameEditFieldLabel.Enable = 'off';
            app.CurrentFrameEditFieldLabel.Position = [373 82 83 22];
            app.CurrentFrameEditFieldLabel.Text = 'Current Frame';

            % Create CurrentFrameEditField
            app.CurrentFrameEditField = uieditfield(app.PTVResearchUIFigure, 'numeric');
            app.CurrentFrameEditField.Limits = [1 Inf];
            app.CurrentFrameEditField.ValueDisplayFormat = '%.0f';
            app.CurrentFrameEditField.ValueChangedFcn = createCallbackFcn(app, @CurrentFrameEditFieldValueChanged, true);
            app.CurrentFrameEditField.Enable = 'off';
            app.CurrentFrameEditField.Position = [371 107 87 22];
            app.CurrentFrameEditField.Value = 1;

            % Create DisplayResultsSwitchLabel
            app.DisplayResultsSwitchLabel = uilabel(app.PTVResearchUIFigure);
            app.DisplayResultsSwitchLabel.HorizontalAlignment = 'center';
            app.DisplayResultsSwitchLabel.Enable = 'off';
            app.DisplayResultsSwitchLabel.Position = [33 32 88 22];
            app.DisplayResultsSwitchLabel.Text = 'Display Results';

            % Create DisplayResultsSwitch
            app.DisplayResultsSwitch = uiswitch(app.PTVResearchUIFigure, 'rocker');
            app.DisplayResultsSwitch.Orientation = 'horizontal';
            app.DisplayResultsSwitch.Enable = 'off';
            app.DisplayResultsSwitch.Position = [57 61 51 23];
            app.DisplayResultsSwitch.Value = 'On';

            % Create DisplayVariableDropDown
            app.DisplayVariableDropDown = uidropdown(app.PTVResearchUIFigure);
            app.DisplayVariableDropDown.Items = {'Vector Length', 'U', 'V', 'Quiver'};
            app.DisplayVariableDropDown.Editable = 'on';
            app.DisplayVariableDropDown.Enable = 'off';
            app.DisplayVariableDropDown.BackgroundColor = [1 1 1];
            app.DisplayVariableDropDown.Position = [154 53 138 22];
            app.DisplayVariableDropDown.Value = 'U';

            % Create STOPButton
            app.STOPButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.STOPButton.ButtonPushedFcn = createCallbackFcn(app, @STOPButtonPushed, true);
            app.STOPButton.BackgroundColor = [1 0 0];
            app.STOPButton.Enable = 'off';
            app.STOPButton.Position = [645 30 100 22];
            app.STOPButton.Text = 'STOP';

            % Create PROLabel
            app.PROLabel = uilabel(app.PTVResearchUIFigure);
            app.PROLabel.HorizontalAlignment = 'right';
            app.PROLabel.FontSize = 28;
            app.PROLabel.FontColor = [0.9412 0.9412 0.9412];
            app.PROLabel.Position = [531 150 100 37];
            app.PROLabel.Text = 'PRO';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.PTVResearchUIFigure);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Position = [528 198 343 409];

            % Create PreTab
            app.PreTab = uitab(app.TabGroup);
            app.PreTab.Title = 'Pre -';

            % Create BackgroundsubtractionButton
            app.BackgroundsubtractionButton = uibutton(app.PreTab, 'state');
            app.BackgroundsubtractionButton.ValueChangedFcn = createCallbackFcn(app, @BackgroundsubtractionButtonValueChanged, true);
            app.BackgroundsubtractionButton.Enable = 'off';
            app.BackgroundsubtractionButton.Text = 'Background subtraction';
            app.BackgroundsubtractionButton.Position = [13 113 145 22];

            % Create ofimagesformeanEditFieldLabel
            app.ofimagesformeanEditFieldLabel = uilabel(app.PreTab);
            app.ofimagesformeanEditFieldLabel.HorizontalAlignment = 'right';
            app.ofimagesformeanEditFieldLabel.Enable = 'off';
            app.ofimagesformeanEditFieldLabel.Position = [161 113 124 22];
            app.ofimagesformeanEditFieldLabel.Text = '% of images for mean';

            % Create ofimagesformeanEditField
            app.ofimagesformeanEditField = uieditfield(app.PreTab, 'numeric');
            app.ofimagesformeanEditField.ValueDisplayFormat = '%.0f';
            app.ofimagesformeanEditField.ValueChangedFcn = createCallbackFcn(app, @PreviewButtonPushed, true);
            app.ofimagesformeanEditField.Enable = 'off';
            app.ofimagesformeanEditField.Position = [291 113 38 22];
            app.ofimagesformeanEditField.Value = 10;

            % Create SharpenButton
            app.SharpenButton = uibutton(app.PreTab, 'state');
            app.SharpenButton.ValueChangedFcn = createCallbackFcn(app, @SharpenButtonValueChanged, true);
            app.SharpenButton.Enable = 'off';
            app.SharpenButton.Text = 'Sharpen';
            app.SharpenButton.Position = [13 79 83 22];

            % Create amountEditFieldLabel
            app.amountEditFieldLabel = uilabel(app.PreTab);
            app.amountEditFieldLabel.HorizontalAlignment = 'right';
            app.amountEditFieldLabel.Enable = 'off';
            app.amountEditFieldLabel.Position = [220 80 62 22];
            app.amountEditFieldLabel.Text = '% amount';

            % Create amountEditField
            app.amountEditField = uieditfield(app.PreTab, 'numeric');
            app.amountEditField.ValueDisplayFormat = '%.0f';
            app.amountEditField.Enable = 'off';
            app.amountEditField.Position = [288 80 38 22];
            app.amountEditField.Value = 2;

            % Create MaskButton
            app.MaskButton = uibutton(app.PreTab, 'push');
            app.MaskButton.ButtonPushedFcn = createCallbackFcn(app, @MaskButtonPushed, true);
            app.MaskButton.BackgroundColor = [0 1 0];
            app.MaskButton.Position = [11 173 84 22];
            app.MaskButton.Text = 'Mask';

            % Create ImageDewarping_________________Label
            app.ImageDewarping_________________Label = uilabel(app.PreTab);
            app.ImageDewarping_________________Label.Position = [9 294 331 22];
            app.ImageDewarping_________________Label.Text = '__________________  Image Dewarping   _________________';

            % Create ResetMaskButton
            app.ResetMaskButton = uibutton(app.PreTab, 'push');
            app.ResetMaskButton.ButtonPushedFcn = createCallbackFcn(app, @ResetMaskButtonPushed, true);
            app.ResetMaskButton.BackgroundColor = [0.851 0.3294 0.102];
            app.ResetMaskButton.Position = [230 173 87 22];
            app.ResetMaskButton.Text = 'Reset Mask';

            % Create ImagePreprocessing________________Label
            app.ImagePreprocessing________________Label = uilabel(app.PreTab);
            app.ImagePreprocessing________________Label.Position = [9 143 334 22];
            app.ImagePreprocessing________________Label.Text = '________________ Image Pre-processing   ________________';

            % Create WienerFilterButton
            app.WienerFilterButton = uibutton(app.PreTab, 'state');
            app.WienerFilterButton.ValueChangedFcn = createCallbackFcn(app, @WienerFilterButtonValueChanged, true);
            app.WienerFilterButton.Text = 'Wiener Filter';
            app.WienerFilterButton.Position = [13 45 98 22];

            % Create kernalsizeEditFieldLabel
            app.kernalsizeEditFieldLabel = uilabel(app.PreTab);
            app.kernalsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.kernalsizeEditFieldLabel.Position = [215 47 62 22];
            app.kernalsizeEditFieldLabel.Text = 'kernal size';

            % Create kernalsizeEditField
            app.kernalsizeEditField = uieditfield(app.PreTab, 'numeric');
            app.kernalsizeEditField.Limits = [1 20];
            app.kernalsizeEditField.ValueDisplayFormat = '%.0f';
            app.kernalsizeEditField.Position = [288 47 38 22];
            app.kernalsizeEditField.Value = 2;

            % Create amountEditField_2Label
            app.amountEditField_2Label = uilabel(app.PreTab);
            app.amountEditField_2Label.HorizontalAlignment = 'right';
            app.amountEditField_2Label.Position = [221 10 62 22];
            app.amountEditField_2Label.Text = '% amount';

            % Create amountEditField_2
            app.amountEditField_2 = uieditfield(app.PreTab, 'numeric');
            app.amountEditField_2.Limits = [0 1];
            app.amountEditField_2.ValueDisplayFormat = '%.2f';
            app.amountEditField_2.Position = [289 10 38 22];
            app.amountEditField_2.Value = 0.1;

            % Create thresholdEditFieldLabel
            app.thresholdEditFieldLabel = uilabel(app.PreTab);
            app.thresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.thresholdEditFieldLabel.Position = [121 10 56 22];
            app.thresholdEditFieldLabel.Text = 'threshold';

            % Create thresholdEditField
            app.thresholdEditField = uieditfield(app.PreTab, 'numeric');
            app.thresholdEditField.Limits = [0 1];
            app.thresholdEditField.ValueDisplayFormat = '%.2f';
            app.thresholdEditField.Position = [183 10 38 22];
            app.thresholdEditField.Value = 0.1;

            % Create LocalIntensityButton
            app.LocalIntensityButton = uibutton(app.PreTab, 'state');
            app.LocalIntensityButton.ValueChangedFcn = createCallbackFcn(app, @LocalIntensityButtonValueChanged, true);
            app.LocalIntensityButton.Text = 'Local Intensity';
            app.LocalIntensityButton.Position = [13 10 98 22];

            % Create LoadMaskButton
            app.LoadMaskButton = uibutton(app.PreTab, 'push');
            app.LoadMaskButton.ButtonPushedFcn = createCallbackFcn(app, @LoadMaskButtonPushed2, true);
            app.LoadMaskButton.BackgroundColor = [0 1 0];
            app.LoadMaskButton.Position = [122 173 84 22];
            app.LoadMaskButton.Text = 'Load Mask';

            % Create LoadparamsfrompreviousButton
            app.LoadparamsfrompreviousButton = uibutton(app.PreTab, 'push');
            app.LoadparamsfrompreviousButton.ButtonPushedFcn = createCallbackFcn(app, @LoadparamsfrompreviousButtonPushed, true);
            app.LoadparamsfrompreviousButton.Position = [12 324 163 22];
            app.LoadparamsfrompreviousButton.Text = 'Load params from previous';

            % Create General______________________Label
            app.General______________________Label = uilabel(app.PreTab);
            app.General______________________Label.Position = [7 361 327 22];
            app.General______________________Label.Text = '_____________________  General   ______________________';

            % Create ImageMasking___________________Label_2
            app.ImageMasking___________________Label_2 = uilabel(app.PreTab);
            app.ImageMasking___________________Label_2.Position = [13 206 331 22];
            app.ImageMasking___________________Label_2.Text = '__________________  Image Masking   ___________________';

            % Create ImageDewarpingButton
            app.ImageDewarpingButton = uibutton(app.PreTab, 'state');
            app.ImageDewarpingButton.ValueChangedFcn = createCallbackFcn(app, @ImageDewarpingButtonValueChanged, true);
            app.ImageDewarpingButton.Text = 'Image Dewarping';
            app.ImageDewarpingButton.Position = [15 266 110 22];

            % Create GridspacingmmEditFieldLabel
            app.GridspacingmmEditFieldLabel = uilabel(app.PreTab);
            app.GridspacingmmEditFieldLabel.HorizontalAlignment = 'right';
            app.GridspacingmmEditFieldLabel.Position = [157 265 102 22];
            app.GridspacingmmEditFieldLabel.Text = 'Grid spacing /mm';

            % Create GridspacingmmEditField
            app.GridspacingmmEditField = uieditfield(app.PreTab, 'numeric');
            app.GridspacingmmEditField.Limits = [0 Inf];
            app.GridspacingmmEditField.Position = [264 265 33 22];
            app.GridspacingmmEditField.Value = 10;

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.PreTab);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.Position = [20 233 47 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create MethodDropDown
            app.MethodDropDown = uidropdown(app.PreTab);
            app.MethodDropDown.Items = {'PWL', 'Poly2', 'Poly3'};
            app.MethodDropDown.ValueChangedFcn = createCallbackFcn(app, @MethodDropDownValueChanged, true);
            app.MethodDropDown.Position = [82 233 100 22];
            app.MethodDropDown.Value = 'Poly2';

            % Create ProcessTab
            app.ProcessTab = uitab(app.TabGroup);
            app.ProcessTab.Title = 'Process';

            % Create GridSizeLabel
            app.GridSizeLabel = uilabel(app.ProcessTab);
            app.GridSizeLabel.Position = [22 75 54 22];
            app.GridSizeLabel.Text = 'Grid Size';

            % Create MedianFilterButton
            app.MedianFilterButton = uibutton(app.ProcessTab, 'state');
            app.MedianFilterButton.Text = 'Median Filter';
            app.MedianFilterButton.Position = [11 15 95 22];

            % Create xEditFieldLabel
            app.xEditFieldLabel = uilabel(app.ProcessTab);
            app.xEditFieldLabel.HorizontalAlignment = 'right';
            app.xEditFieldLabel.Position = [105 74 25 22];
            app.xEditFieldLabel.Text = 'x';

            % Create xEditField
            app.xEditField = uieditfield(app.ProcessTab, 'numeric');
            app.xEditField.ValueChangedFcn = createCallbackFcn(app, @xEditFieldValueChanged, true);
            app.xEditField.Position = [134 74 35 22];

            % Create yEditFieldLabel
            app.yEditFieldLabel = uilabel(app.ProcessTab);
            app.yEditFieldLabel.HorizontalAlignment = 'right';
            app.yEditFieldLabel.Position = [171 75 25 22];
            app.yEditFieldLabel.Text = 'y';

            % Create yEditField
            app.yEditField = uieditfield(app.ProcessTab, 'numeric');
            app.yEditField.ValueChangedFcn = createCallbackFcn(app, @yEditFieldValueChanged, true);
            app.yEditField.Position = [199 75 31 22];

            % Create PxmmEditFieldLabel
            app.PxmmEditFieldLabel = uilabel(app.ProcessTab);
            app.PxmmEditFieldLabel.HorizontalAlignment = 'right';
            app.PxmmEditFieldLabel.Position = [6 322 51 22];
            app.PxmmEditFieldLabel.Text = 'Px / mm';

            % Create PxmmEditField
            app.PxmmEditField = uieditfield(app.ProcessTab, 'numeric');
            app.PxmmEditField.Position = [68 322 52 22];
            app.PxmmEditField.Value = 1;

            % Create offramessecLabel
            app.offramessecLabel = uilabel(app.ProcessTab);
            app.offramessecLabel.Position = [152 322 89 22];
            app.offramessecLabel.Text = '# of frames/sec';

            % Create offramessecEditField
            app.offramessecEditField = uieditfield(app.ProcessTab, 'numeric');
            app.offramessecEditField.ValueDisplayFormat = '%.0f';
            app.offramessecEditField.Position = [244 322 50 22];

            % Create InterpolatorDropDownLabel
            app.InterpolatorDropDownLabel = uilabel(app.ProcessTab);
            app.InterpolatorDropDownLabel.HorizontalAlignment = 'right';
            app.InterpolatorDropDownLabel.Position = [15 44 67 22];
            app.InterpolatorDropDownLabel.Text = 'Interpolator';

            % Create InterpolatorDropDown
            app.InterpolatorDropDown = uidropdown(app.ProcessTab);
            app.InterpolatorDropDown.Items = {'linear', 'nearest', 'cubic', 'natural'};
            app.InterpolatorDropDown.Position = [129 44 101 22];
            app.InterpolatorDropDown.Value = 'cubic';

            % Create KernelsizeEditFieldLabel
            app.KernelsizeEditFieldLabel = uilabel(app.ProcessTab);
            app.KernelsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.KernelsizeEditFieldLabel.Position = [121 15 62 22];
            app.KernelsizeEditFieldLabel.Text = 'Kernel size';

            % Create KernelsizeEditField
            app.KernelsizeEditField = uieditfield(app.ProcessTab, 'numeric');
            app.KernelsizeEditField.Position = [192 15 37 22];
            app.KernelsizeEditField.Value = 2;

            % Create CalculateLagrangianTracksButton
            app.CalculateLagrangianTracksButton = uibutton(app.ProcessTab, 'state');
            app.CalculateLagrangianTracksButton.ValueChangedFcn = createCallbackFcn(app, @CalculateLagrangianTracksButtonValueChanged, true);
            app.CalculateLagrangianTracksButton.Text = 'Calculate Lagrangian Tracks';
            app.CalculateLagrangianTracksButton.Position = [16 143 167 22];

            % Create ScalingParameters_______________Label
            app.ScalingParameters_______________Label = uilabel(app.ProcessTab);
            app.ScalingParameters_______________Label.Position = [9 353 313 22];
            app.ScalingParameters_______________Label.Text = '________________ Scaling Parameters   _______________';

            % Create LagrangianParameters______________Label
            app.LagrangianParameters______________Label = uilabel(app.ProcessTab);
            app.LagrangianParameters______________Label.Position = [15 182 319 22];
            app.LagrangianParameters______________Label.Text = '______________  Lagrangian Parameters   ______________';

            % Create GriddedParameters_______________Label
            app.GriddedParameters_______________Label = uilabel(app.ProcessTab);
            app.GriddedParameters_______________Label.Position = [13 106 321 22];
            app.GriddedParameters_______________Label.Text = '________________  Gridded Parameters   _______________';

            % Create ScalefromreferenceButton
            app.ScalefromreferenceButton = uibutton(app.ProcessTab, 'push');
            app.ScalefromreferenceButton.ButtonPushedFcn = createCallbackFcn(app, @ScalefromreferenceButtonPushed, true);
            app.ScalefromreferenceButton.Position = [15 287 130 22];
            app.ScalefromreferenceButton.Text = 'Scale from reference ';

            % Create ProcessingMode_________________Label
            app.ProcessingMode_________________Label = uilabel(app.ProcessTab);
            app.ProcessingMode_________________Label.Position = [13 254 319 24];
            app.ProcessingMode_________________Label.Text = '_______________   Processing Mode   _________________';

            % Create PTVButton
            app.PTVButton = uibutton(app.ProcessTab, 'state');
            app.PTVButton.ValueChangedFcn = createCallbackFcn(app, @PTVButtonValueChanged, true);
            app.PTVButton.Text = 'PTV';
            app.PTVButton.Position = [15 220 100 22];
            app.PTVButton.Value = true;

            % Create SnapshotPIVButton
            app.SnapshotPIVButton = uibutton(app.ProcessTab, 'state');
            app.SnapshotPIVButton.ValueChangedFcn = createCallbackFcn(app, @SnapshotPIVButtonValueChanged, true);
            app.SnapshotPIVButton.Text = 'Snapshot PIV';
            app.SnapshotPIVButton.Position = [144 220 100 22];

            % Create Label_5
            app.Label_5 = uilabel(app.ProcessTab);
            app.Label_5.Position = [73 287 25 22];
            app.Label_5.Text = '';

            % Create PostGridTab
            app.PostGridTab = uitab(app.TabGroup);
            app.PostGridTab.Title = 'Post Grid -';

            % Create PODDEMFilterButton
            app.PODDEMFilterButton = uibutton(app.PostGridTab, 'push');
            app.PODDEMFilterButton.ButtonPushedFcn = createCallbackFcn(app, @PODDEMFilterButtonPushed, true);
            app.PODDEMFilterButton.BackgroundColor = [0 1 1];
            app.PODDEMFilterButton.Position = [11 315 100 22];
            app.PODDEMFilterButton.Text = 'PODDEM Filter';

            % Create PODButton
            app.PODButton = uibutton(app.PostGridTab, 'push');
            app.PODButton.ButtonPushedFcn = createCallbackFcn(app, @PODButtonPushed, true);
            app.PODButton.BackgroundColor = [0.4706 0.6706 0.1882];
            app.PODButton.Position = [14 138 50 22];
            app.PODButton.Text = 'POD';

            % Create MeanFieldsButton
            app.MeanFieldsButton = uibutton(app.PostGridTab, 'state');
            app.MeanFieldsButton.Text = 'Mean Fields';
            app.MeanFieldsButton.BackgroundColor = [0.302 0.749 0.9294];
            app.MeanFieldsButton.Position = [13 210 93 22];

            % Create RANSQuantitiesButton
            app.RANSQuantitiesButton = uibutton(app.PostGridTab, 'state');
            app.RANSQuantitiesButton.Text = 'RANS Quantities';
            app.RANSQuantitiesButton.BackgroundColor = [0.302 0.749 0.9294];
            app.RANSQuantitiesButton.Position = [118 210 101 22];

            % Create SwirlVorticityButton
            app.SwirlVorticityButton = uibutton(app.PostGridTab, 'state');
            app.SwirlVorticityButton.Text = 'Swirl & Vorticity';
            app.SwirlVorticityButton.BackgroundColor = [0.302 0.749 0.9294];
            app.SwirlVorticityButton.Position = [230 210 97 22];

            % Create Postvariable_mean
            app.Postvariable_mean = uidropdown(app.PostGridTab);
            app.Postvariable_mean.Items = {'U', 'V', 'Vector Length', 'Voricity', 'Quiver', 'uu', 'uv', 'vv', 'TKE', 'Fluc Vorticity', 'Swirl'};
            app.Postvariable_mean.Editable = 'on';
            app.Postvariable_mean.ValueChangedFcn = createCallbackFcn(app, @Postvariable_meanValueChanged, true);
            app.Postvariable_mean.BackgroundColor = [1 1 1];
            app.Postvariable_mean.Position = [121 17 99 22];
            app.Postvariable_mean.Value = 'U';

            % Create DisplayButton_2
            app.DisplayButton_2 = uibutton(app.PostGridTab, 'push');
            app.DisplayButton_2.ButtonPushedFcn = createCallbackFcn(app, @DisplayButton_2Pushed, true);
            app.DisplayButton_2.Position = [228 17 100 22];
            app.DisplayButton_2.Text = 'Display';

            % Create ResetButton
            app.ResetButton = uibutton(app.PostGridTab, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.BackgroundColor = [1 0 0];
            app.ResetButton.Position = [264 315 63 30];
            app.ResetButton.Text = 'Reset';

            % Create KernelEditFieldLabel
            app.KernelEditFieldLabel = uilabel(app.PostGridTab);
            app.KernelEditFieldLabel.HorizontalAlignment = 'right';
            app.KernelEditFieldLabel.Position = [124 280 40 22];
            app.KernelEditFieldLabel.Text = 'Kernel';

            % Create KernelEditField
            app.KernelEditField = uieditfield(app.PostGridTab, 'numeric');
            app.KernelEditField.Position = [174 280 38 22];
            app.KernelEditField.Value = 2;

            % Create ShowDataLabel
            app.ShowDataLabel = uilabel(app.PostGridTab);
            app.ShowDataLabel.Position = [52 17 64 22];
            app.ShowDataLabel.Text = 'Show Data';

            % Create SmoothingButton
            app.SmoothingButton = uibutton(app.PostGridTab, 'push');
            app.SmoothingButton.ButtonPushedFcn = createCallbackFcn(app, @SmoothingButtonPushed, true);
            app.SmoothingButton.BackgroundColor = [0 1 1];
            app.SmoothingButton.Position = [11 280 100 22];
            app.SmoothingButton.Text = 'Smoothing';

            % Create DMDButton
            app.DMDButton = uibutton(app.PostGridTab, 'push');
            app.DMDButton.BackgroundColor = [0.4706 0.6706 0.1882];
            app.DMDButton.Enable = 'off';
            app.DMDButton.Position = [75 138 43 22];
            app.DMDButton.Text = 'DMD';

            % Create Label_4
            app.Label_4 = uilabel(app.PostGridTab);
            app.Label_4.Position = [9 89 25 22];
            app.Label_4.Text = '';

            % Create QuantitiestoCalculate_______________Label
            app.QuantitiestoCalculate_______________Label = uilabel(app.PostGridTab);
            app.QuantitiestoCalculate_______________Label.Position = [11 243 319 22];
            app.QuantitiestoCalculate_______________Label.Text = '_____________   Quantities to Calculate   _______________';

            % Create ModalDecompositions_______________Label
            app.ModalDecompositions_______________Label = uilabel(app.PostGridTab);
            app.ModalDecompositions_______________Label.Position = [10 174 318 22];
            app.ModalDecompositions_______________Label.Text = '_____________   Modal Decompositions  _______________';

            % Create Label_3
            app.Label_3 = uilabel(app.PostGridTab);
            app.Label_3.Position = [11 114 321 22];
            app.Label_3.Text = '____________________________________________________';

            % Create Filters______________________Label
            app.Filters______________________Label = uilabel(app.PostGridTab);
            app.Filters______________________Label.Position = [9 353 320 22];
            app.Filters______________________Label.Text = '_____________________   Filters   ______________________';

            % Create PostLagTab
            app.PostLagTab = uitab(app.TabGroup);
            app.PostLagTab.Title = 'Post Lag -';

            % Create MinpathlengthEditFieldLabel
            app.MinpathlengthEditFieldLabel = uilabel(app.PostLagTab);
            app.MinpathlengthEditFieldLabel.HorizontalAlignment = 'right';
            app.MinpathlengthEditFieldLabel.Position = [179 318 89 22];
            app.MinpathlengthEditFieldLabel.Text = 'Min path length';

            % Create MinpathlengthEditField
            app.MinpathlengthEditField = uieditfield(app.PostLagTab, 'numeric');
            app.MinpathlengthEditField.Limits = [1 Inf];
            app.MinpathlengthEditField.ValueDisplayFormat = '%.0f';
            app.MinpathlengthEditField.Position = [283 318 40 22];
            app.MinpathlengthEditField.Value = 2;

            % Create trackstodisplayEditFieldLabel
            app.trackstodisplayEditFieldLabel = uilabel(app.PostLagTab);
            app.trackstodisplayEditFieldLabel.HorizontalAlignment = 'right';
            app.trackstodisplayEditFieldLabel.Position = [165 52 109 22];
            app.trackstodisplayEditFieldLabel.Text = '% tracks to display';

            % Create trackstodisplayEditField
            app.trackstodisplayEditField = uieditfield(app.PostLagTab, 'numeric');
            app.trackstodisplayEditField.Limits = [1 100];
            app.trackstodisplayEditField.ValueDisplayFormat = '%.0f';
            app.trackstodisplayEditField.Position = [286 52 40 22];
            app.trackstodisplayEditField.Value = 2;

            % Create DisplayButton
            app.DisplayButton = uibutton(app.PostLagTab, 'push');
            app.DisplayButton.ButtonPushedFcn = createCallbackFcn(app, @DisplayButtonPushed, true);
            app.DisplayButton.Position = [228 17 100 22];
            app.DisplayButton.Text = 'Display';

            % Create LagrangianPathsButton
            app.LagrangianPathsButton = uibutton(app.PostLagTab, 'state');
            app.LagrangianPathsButton.Text = 'Lagrangian Paths';
            app.LagrangianPathsButton.Position = [12 317 109 22];

            % Create DiffusionCoefficientsButton
            app.DiffusionCoefficientsButton = uibutton(app.PostLagTab, 'state');
            app.DiffusionCoefficientsButton.Text = 'Diffusion Coefficients';
            app.DiffusionCoefficientsButton.Position = [12 280 129 22];

            % Create QuantitiestoCalculate_______________Label_2
            app.QuantitiestoCalculate_______________Label_2 = uilabel(app.PostLagTab);
            app.QuantitiestoCalculate_______________Label_2.Position = [11 353 319 22];
            app.QuantitiestoCalculate_______________Label_2.Text = '_____________   Quantities to Calculate   _______________';

            % Create PostvariableLag
            app.PostvariableLag = uidropdown(app.PostLagTab);
            app.PostvariableLag.Items = {};
            app.PostvariableLag.Editable = 'on';
            app.PostvariableLag.BackgroundColor = [1 1 1];
            app.PostvariableLag.Position = [121 17 99 22];
            app.PostvariableLag.Value = {};

            % Create ShowDataLabelLag
            app.ShowDataLabelLag = uilabel(app.PostLagTab);
            app.ShowDataLabelLag.Position = [52 17 64 22];
            app.ShowDataLabelLag.Text = 'Show Data';

            % Create LoadvectorfieldsButton
            app.LoadvectorfieldsButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.LoadvectorfieldsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadvectorfieldsButtonPushed, true);
            app.LoadvectorfieldsButton.Position = [764 628 112 22];
            app.LoadvectorfieldsButton.Text = 'Load vector fields';

            % Create Postvariable
            app.Postvariable = uidropdown(app.PTVResearchUIFigure);
            app.Postvariable.Items = {'U', 'V', 'Vector Length', 'Voricity', 'Quiver', 'uu', 'uv', 'vv', 'TKE', 'Fluc Vorticity', 'Swirl'};
            app.Postvariable.Editable = 'on';
            app.Postvariable.ValueChangedFcn = createCallbackFcn(app, @PostvariableValueChanged, true);
            app.Postvariable.BackgroundColor = [1 1 1];
            app.Postvariable.Position = [154 53 138 22];
            app.Postvariable.Value = 'U';

            % Create Slider_pod
            app.Slider_pod = uislider(app.PTVResearchUIFigure);
            app.Slider_pod.MajorTickLabels = {'', '', '', '', '', ''};
            app.Slider_pod.ValueChangedFcn = createCallbackFcn(app, @Slider_podValueChanged, true);
            app.Slider_pod.ValueChangingFcn = createCallbackFcn(app, @Slider_podValueChanging, true);
            app.Slider_pod.Position = [57 155 403 3];

            % Create DisplayVariableLabel
            app.DisplayVariableLabel = uilabel(app.PTVResearchUIFigure);
            app.DisplayVariableLabel.Position = [175 29 91 22];
            app.DisplayVariableLabel.Text = 'Display Variable';

            % Create NumbertoprocessEditFieldLabel
            app.NumbertoprocessEditFieldLabel = uilabel(app.PTVResearchUIFigure);
            app.NumbertoprocessEditFieldLabel.HorizontalAlignment = 'right';
            app.NumbertoprocessEditFieldLabel.Enable = 'off';
            app.NumbertoprocessEditFieldLabel.Position = [648 129 109 22];
            app.NumbertoprocessEditFieldLabel.Text = 'Number to process';

            % Create NumbertoprocessEditField
            app.NumbertoprocessEditField = uieditfield(app.PTVResearchUIFigure, 'numeric');
            app.NumbertoprocessEditField.ValueDisplayFormat = '%.0f';
            app.NumbertoprocessEditField.ValueChangedFcn = createCallbackFcn(app, @NumbertoprocessEditFieldValueChanged, true);
            app.NumbertoprocessEditField.Enable = 'off';
            app.NumbertoprocessEditField.Position = [772 129 95 22];
            app.NumbertoprocessEditField.Value = 100;

            % Create StartFrameEditFieldLabel
            app.StartFrameEditFieldLabel = uilabel(app.PTVResearchUIFigure);
            app.StartFrameEditFieldLabel.HorizontalAlignment = 'right';
            app.StartFrameEditFieldLabel.Enable = 'off';
            app.StartFrameEditFieldLabel.Position = [687 156 69 22];
            app.StartFrameEditFieldLabel.Text = 'Start Frame';

            % Create StartFrameEditField
            app.StartFrameEditField = uieditfield(app.PTVResearchUIFigure, 'numeric');
            app.StartFrameEditField.Limits = [1 Inf];
            app.StartFrameEditField.ValueDisplayFormat = '%.0f';
            app.StartFrameEditField.ValueChangedFcn = createCallbackFcn(app, @StartFrameEditFieldValueChanged, true);
            app.StartFrameEditField.Enable = 'off';
            app.StartFrameEditField.Position = [772 156 93 22];
            app.StartFrameEditField.Value = 1;

            % Create EndFrameEditField
            app.EndFrameEditField = uieditfield(app.PTVResearchUIFigure, 'numeric');
            app.EndFrameEditField.ValueDisplayFormat = '%.0f';
            app.EndFrameEditField.Enable = 'off';
            app.EndFrameEditField.Position = [771 103 95 22];
            app.EndFrameEditField.Value = 1;

            % Create EndFrameEditFieldLabel
            app.EndFrameEditFieldLabel = uilabel(app.PTVResearchUIFigure);
            app.EndFrameEditFieldLabel.HorizontalAlignment = 'right';
            app.EndFrameEditFieldLabel.Enable = 'off';
            app.EndFrameEditFieldLabel.Position = [695 103 64 22];
            app.EndFrameEditFieldLabel.Text = 'End Frame';

            % Create Forward
            app.Forward = uibutton(app.PTVResearchUIFigure, 'push');
            app.Forward.ButtonPushedFcn = createCallbackFcn(app, @ForwardButtonPushed, true);
            app.Forward.Position = [92 106 28 27];
            app.Forward.Text = '>';

            % Create Backward
            app.Backward = uibutton(app.PTVResearchUIFigure, 'push');
            app.Backward.ButtonPushedFcn = createCallbackFcn(app, @BackwardButtonPushed, true);
            app.Backward.Position = [62 106 28 27];
            app.Backward.Text = '<';

            % Create PlayButton
            app.PlayButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.Position = [133 105 49 27];
            app.PlayButton.Text = 'Play';

            % Create PauseButton
            app.PauseButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.PauseButton.ButtonPushedFcn = createCallbackFcn(app, @PauseButtonPushed, true);
            app.PauseButton.Position = [196 104 46 27];
            app.PauseButton.Text = 'Pause';

            % Create ColourbarLabel
            app.ColourbarLabel = uilabel(app.PTVResearchUIFigure);
            app.ColourbarLabel.Position = [21 639 59 22];
            app.ColourbarLabel.Text = 'Colourbar';

            % Create DropDown
            app.DropDown = uidropdown(app.PTVResearchUIFigure);
            app.DropDown.Items = {'auto', 'manual'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            app.DropDown.Position = [19 619 78 22];
            app.DropDown.Value = 'auto';

            % Create minEditFieldLabel
            app.minEditFieldLabel = uilabel(app.PTVResearchUIFigure);
            app.minEditFieldLabel.HorizontalAlignment = 'right';
            app.minEditFieldLabel.Position = [111 619 25 22];
            app.minEditFieldLabel.Text = 'min';

            % Create minEditField
            app.minEditField = uieditfield(app.PTVResearchUIFigure, 'numeric');
            app.minEditField.ValueDisplayFormat = '%.3f';
            app.minEditField.ValueChangedFcn = createCallbackFcn(app, @minEditFieldValueChanged, true);
            app.minEditField.Position = [139 619 58 22];

            % Create maxEditFieldLabel
            app.maxEditFieldLabel = uilabel(app.PTVResearchUIFigure);
            app.maxEditFieldLabel.HorizontalAlignment = 'right';
            app.maxEditFieldLabel.Position = [200 619 28 22];
            app.maxEditFieldLabel.Text = 'max';

            % Create maxEditField
            app.maxEditField = uieditfield(app.PTVResearchUIFigure, 'numeric');
            app.maxEditField.ValueDisplayFormat = '%.3f';
            app.maxEditField.ValueChangedFcn = createCallbackFcn(app, @maxEditFieldValueChanged, true);
            app.maxEditField.Position = [232 619 59 22];

            % Create PlaybackSpeedSliderLabel
            app.PlaybackSpeedSliderLabel = uilabel(app.PTVResearchUIFigure);
            app.PlaybackSpeedSliderLabel.HorizontalAlignment = 'right';
            app.PlaybackSpeedSliderLabel.Position = [251 93 93 22];
            app.PlaybackSpeedSliderLabel.Text = 'Playback Speed';

            % Create PlaybackSpeedSlider
            app.PlaybackSpeedSlider = uislider(app.PTVResearchUIFigure);
            app.PlaybackSpeedSlider.Limits = [1 10];
            app.PlaybackSpeedSlider.MajorTicks = [];
            app.PlaybackSpeedSlider.MajorTickLabels = {};
            app.PlaybackSpeedSlider.ValueChangedFcn = createCallbackFcn(app, @PlaybackSpeedSliderValueChanged, true);
            app.PlaybackSpeedSlider.MinorTicks = [];
            app.PlaybackSpeedSlider.Position = [259 118 85 3];
            app.PlaybackSpeedSlider.Value = 5;

            % Create LoadrawimagesButton
            app.LoadrawimagesButton = uibutton(app.PTVResearchUIFigure, 'push');
            app.LoadrawimagesButton.ButtonPushedFcn = createCallbackFcn(app, @LoadrawimagesButtonPushed, true);
            app.LoadrawimagesButton.Position = [642 628 108 22];
            app.LoadrawimagesButton.Text = 'Load raw images';

            % Create Lagvariable
            app.Lagvariable = uidropdown(app.PTVResearchUIFigure);
            app.Lagvariable.Items = {'Points vec length', 'Points u-vel', 'Points v-vel'};
            app.Lagvariable.Editable = 'on';
            app.Lagvariable.ValueChangedFcn = createCallbackFcn(app, @LagvariableValueChanged, true);
            app.Lagvariable.BackgroundColor = [1 1 1];
            app.Lagvariable.Position = [154 52 138 22];
            app.Lagvariable.Value = 'Points vec length';

            % Create UIAxes
            app.UIAxes = uiaxes(app.PTVResearchUIFigure);
            app.UIAxes.PlotBoxAspectRatio = [1 0.823045267489712 0.823045267489712];
            app.UIAxes.Box = 'on';
            app.UIAxes.Position = [4 185 516 427];
        end
    end

    methods (Access = public)

        % Construct app
        function app = PTVResearch

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PTVResearchUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PTVResearchUIFigure)
        end
    end
end