# gui_taus_matlab
GUI using app designer in matlab for relaxation analysis

classdef analysisAD < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure                   % UI Figure
        TabGroup               matlab.ui.container.TabGroup       % Prepari...
        prep                   matlab.ui.container.Tab            % Preparing
        UIAxes1                matlab.ui.control.UIAxes           % Title
        init                   matlab.ui.container.Panel          % Initial...
        LabelSlider            matlab.ui.control.Label            % Begin
        cutBegin               matlab.ui.control.Slider           % [1 100]
        Label                  matlab.ui.control.Label            % End
        cutEnd                 matlab.ui.control.Slider           % [0 100]
        setP                   matlab.ui.control.Button           % Set
        findk                  matlab.ui.control.Button           % Find K's
        findpulse              matlab.ui.control.Button           % Find pu...
        qvstime                matlab.ui.control.Button           % SeeQ
        pulsenumtext           matlab.ui.control.TextArea        
        qends                  matlab.ui.control.Button           % Qends
        LabelNumericEditField4 matlab.ui.control.Label            % Poly de...
        polydegree             matlab.ui.control.NumericEditField % [1 10]
        applypoly              matlab.ui.control.Button           % Apply poly
        OpenFile               matlab.ui.control.Button           % open file
        LabelNumericEditField6 matlab.ui.control.Label            % SFluid ...
        sfluid                 matlab.ui.control.NumericEditField % [0 1]
        LabelTextArea          matlab.ui.control.Label            % K values
        textfork               matlab.ui.control.TextArea         % K_T = 0
        LabelEditField         matlab.ui.control.Label            % Folder ...
        path                   matlab.ui.control.EditField        % D:\ther...
        tempSF                 matlab.ui.control.Button           % Tempera...
        Tab2                   matlab.ui.container.Tab            % Filtering
        UIAxes2                matlab.ui.control.UIAxes           % Title
        Panel                  matlab.ui.container.Panel          % Pulse s...
        LabelNumericEditField  matlab.ui.control.Label            % Pulse n...
        pulsenum               matlab.ui.control.NumericEditField % [0 Inf]
        LabelDropDown          matlab.ui.control.Label            % Filters
        filt                   matlab.ui.control.DropDown         % Median,...
        LabelNumericEditField2 matlab.ui.control.Label            % 1st fie...
        filt_n1                matlab.ui.control.NumericEditField % [0 Inf]
        apply1                 matlab.ui.control.Button           % Apply f...
        incline                matlab.ui.control.Button           % Incline
        LabelSlider2           matlab.ui.control.Label            % Begin
        cutBegin2              matlab.ui.control.Slider           % [1 100]
        Label2                 matlab.ui.control.Label            % End
        cutEnd2                matlab.ui.control.Slider           % [1 100]
        apply_cut              matlab.ui.control.Button           % Apply cut
        Tab3                   matlab.ui.container.Tab            % Fitting
        UIAxes3                matlab.ui.control.UIAxes           % Title
        Panel2                 matlab.ui.container.Panel          % Filter ...
        LabelDropDown2         matlab.ui.control.Label            % Fitting...
        fitfuns                matlab.ui.control.DropDown         % Fit1, F...
        Label3                 matlab.ui.control.Label            % Fit 1: ...
        Label4                 matlab.ui.control.Label            % Fit 2: ...
        Label5                 matlab.ui.control.Label            % Fit 3: ...
        exportF                matlab.ui.control.Button           % Export ...
        UIAxes4                matlab.ui.control.UIAxes           % Title
        Panel3                 matlab.ui.container.Panel          % Fitting...
        printfit               matlab.ui.control.TextArea         % 1
        Tab4                   matlab.ui.container.Tab            % Residuals
        UIAxes5                matlab.ui.control.UIAxes           % Residuals
        UIAxes6                matlab.ui.control.UIAxes           % Title
        errtext                matlab.ui.control.TextArea         % res_n
        LabelNumericEditField3 matlab.ui.control.Label            % Bins
        nbin                   matlab.ui.control.NumericEditField % [1 Inf]
        clear_err              matlab.ui.control.Button           % Clear
        Closeapp               matlab.ui.control.Button           % Close
        LabelNumericEditField5 matlab.ui.control.Label            % Edit Field
        k_f2                   matlab.ui.control.NumericEditField % [-Inf Inf]
        Tab5                   matlab.ui.container.Tab            % Re-model
        UIAxes7                matlab.ui.control.UIAxes           % Re-model
    end


    properties (Access = public)
        ff; % import data function
        ff1; % data after cutting
        beg1; % cut begin in ini
        end1; % cut form end in ini
        beg2; % cut after filtering
        end2;
        c; % indecies of begin and end pulse
        pulse_one; % Description
        n_pul; % number of pulse
        af; % filtered y axis
        f_inc; % data after incline
        k_q; %incline of Q data
        b_q;
        tim; % time after filtering
        tem; % temperature of pulse
        fit1; % fit 1 parametrs
        fit2; % fit 2 parameters
        fit3;
        fit4;
        poly_off; 
    end

    


    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.path.Value=pwd;
            [file,ruta]=uigetfile({'*.dat'},'Select a dat file for second Fork');
            cd(ruta);
            ff2=importdata(file);
            [b1,~]=size(ff2.data);
            if b1>1
                app.ff=ff2.data;
            else
                app.ff=str2double(ff2.textdata(2:end,2:end));
            end
            app.ff(:,1)=app.ff(:,1)-app.ff(1,1);
            app.ff(:,end)=[];
            %% initial parameters
            cd(app.path.Value)
            app.end1=length(app.ff);
            app.beg1=1;
            app.cutBegin.Limits = [1 app.end1];
            app.cutEnd.Limits = [1 app.end1];
            app.cutEnd.Value = app.end1;
            app.ff1=app.ff;
            Pkb=[0 0];
            app.fit2=[1 0 1000 0];
            app.b_q=0;
            p_n=1;
            %% plot temperature vs time
                        
            app.UIAxes1.cla;
            plot(app.UIAxes1,app.ff(:,1),app.ff(:,12));
            title(app.UIAxes1, 'Print after import');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Temperature [mK]');
            
            mintime=1;
            maxtime=length(app.ff);

        end

        % cutBegin value changed function
        function cutBeginValueChanged(app)
            app.beg1 = round(app.cutBegin.Value);
            app.UIAxes1.cla;
            plot(app.UIAxes1,app.ff(app.beg1:app.end1,1),app.ff(app.beg1:app.end1,12));
            title(app.UIAxes1, 'Modifying');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Temperature [mK]');
        end

        % cutEnd value changed function
        function cutEndValueChanged(app)
            app.end1 = round(app.cutEnd.Value);
           app.UIAxes1.cla;
            plot(app.UIAxes1,app.ff(app.beg1:app.end1,1),app.ff(app.beg1:app.end1,12));
            title(app.UIAxes1, 'Modifying');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Temperature [mK]');
        end

        % setP button pushed function
        function setPButtonPushed(app)
            app.ff1=app.ff;
            app.ff1(app.end1:end,:)=[];
            app.ff1(1:app.beg1,:)=[];
           app.UIAxes1.cla;
            plot(app.UIAxes1,app.ff1(:,1),app.ff1(:,12));
            title(app.UIAxes1, 'Modifyed');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Temperature [mK]');
            ylim(app.UIAxes1,'auto');
            
        end

        % qvstime button pushed function
        function qvstimeButtonPushed(app)
            app.UIAxes1.cla;
            plot(app.UIAxes1,app.ff1(:,1),app.ff1(:,5));
            title(app.UIAxes1, 'Modifyed');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Temperature [mK]');
        end

        % findk button pushed function
        function findkButtonPushed(app)
            Pkb = polyfit(app.ff1(:,1),app.ff1(:,12),1);
            x1=mean(app.ff1(1:100,1));
            x2=mean(app.ff1(end-100:end,1));
            y1=mean(app.ff1(1:100,5));
            y2=mean(app.ff1(end-100:end,5));
            app.k_q=(y2-y1)/(x2-x1);
            st1=num2str(Pkb(1));
            st=['K_T is ',st1];
            st2=num2str(app.k_q);
            stt=['K_Q is ',st2];
            app.textfork.Value = {st; stt};
        end

        % findpulse button pushed function
        function findpulseButtonPushed(app)
            app.c=pulse_finder(app.ff1,10);
            [c1,~]=size(app.c);
            app.pulsenum.Limits = [1 c1];
            st1=num2str(c1);
            app.pulsenumtext.Value={st1};
            app.UIAxes1.cla;
            hold(app.UIAxes1,'on');
            plot(app.UIAxes1,app.ff1(:,1),app.ff1(:,5));
            plot(app.UIAxes1,app.ff1(app.c(:,1),1),app.ff1(app.c(:,1),5),'or');
        end

        % pulsenum value changed function
        function pulsenumValueChanged(app)
            app.n_pul = round(app.pulsenum.Value);
            st1=num2str(app.n_pul);
            app.tem=num2str(app.ff1(app.c(app.n_pul,1),12));
            st=['Unfiltered ',st1,'s pulse; T = ',app.tem,' mK'];           
            app.UIAxes2.cla;
            plot(app.UIAxes2,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1),app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5));
            
            title(app.UIAxes2,st);
            xlabel(app.UIAxes2,'Time [sec]');
            ylabel(app.UIAxes2,'Q');
        end

        % filt value changed function
        function filtValueChanged(app, event)
            val = app.filt.Value;
            switch val
                case 'Median'
                    %app.af=abs(medfilt1(app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5),app.pulse_one));
                    app.af=medfilt1(app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5),app.pulse_one);
                    app.af(1)=app.af(2);
                    app.UIAxes2.cla;
                    hold(app.UIAxes2,'on');
                    plot(app.UIAxes2,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1),app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5));
                    plot(app.UIAxes2,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1),app.af,'-r','LineWidth',1);
                    title(app.UIAxes2,'Median');
                    xlabel(app.UIAxes2,'Time [sec]');
                    ylabel(app.UIAxes2,'Q');
                    app.tim=app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1);
                case 'None'
                    st1=num2str(app.n_pul);
                    app.tem=num2str(app.ff1(app.c(app.n_pul,1),12));
                    st=['Unfiltered ',st1,'s pulse; T = ',app.tem,' mK'];   
                    app.UIAxes2.cla;
                    plot(app.UIAxes2,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1),app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5));
                    title(app.UIAxes2,st);
                    xlabel(app.UIAxes2,'Time [sec]');
                    ylabel(app.UIAxes2,'Q');
                    app.tim=app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1);
                case 'Moving'
                    b=(1/app.pulse_one)*ones(1,app.pulse_one);
                    a=1;
                    app.af=filter(b,a,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5));
                    app.af(1:app.pulse_one-1)=app.af(app.pulse_one);
                    app.UIAxes2.cla;
                    hold(app.UIAxes2,'on');
                    plot(app.UIAxes2,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1),app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),5));
                    plot(app.UIAxes2,app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1),app.af,'-r','LineWidth',1);
                    title(app.UIAxes2,'Moving');
                    xlabel(app.UIAxes2,'Time [sec]');
                    ylabel(app.UIAxes2,'Q');
                    app.tim=app.ff1(app.c(app.n_pul,1):app.c(app.n_pul,2),1);
            end
            
        end

        % filt_n1 value changed function
        function filt_n1ValueChanged(app)
            app.pulse_one = round(app.filt_n1.Value);
            
        end

        % apply1 button pushed function
        function apply1ButtonPushed(app)
            app.filtValueChanged;
        end

        % incline button pushed function
        function inclineButtonPushed(app)
            for ii=1:length(app.af)
                app.f_inc(ii,1)=app.ff1(app.c(app.n_pul,1)+ii-1,1);
                %app.f_inc(ii,2)= app.af(ii)-app.k_q*app.f_inc(ii,1);
                app.f_inc(ii,2)= app.af(ii);
            end
            app.cutBegin2.Limits = [1 length(app.f_inc)];
            app.cutEnd2.Limits = [1 length(app.f_inc)];
            app.cutEnd2.Value = length(app.f_inc);
            app.UIAxes2.cla;
            plot(app.UIAxes2,app.f_inc(:,1),app.f_inc(:,2),'-r','LineWidth',1);
                    title(app.UIAxes2,'After incline');
                    xlabel(app.UIAxes2,'Time [sec]');
                    ylabel(app.UIAxes2,'Q');
        end

        % cutBegin2 value changed function
        function cutBegin2ValueChanged(app)
            app.beg2 = round(app.cutBegin2.Value);
            app.UIAxes2.cla;
            plot(app.UIAxes2,app.f_inc(app.beg2:app.end2,1),app.f_inc(app.beg2:app.end2,2));
            title(app.UIAxes2, 'Modifying');
            xlabel(app.UIAxes2, 'Time [sec]');
            ylabel(app.UIAxes2, 'Q');
        end

        % cutEnd2 value changed function
        function cutEnd2ValueChanged(app)
            app.end2 = round(app.cutEnd2.Value);
            app.UIAxes2.cla;
            plot(app.UIAxes2,app.f_inc(app.beg2:app.end2,1),app.f_inc(app.beg2:app.end2,2));
            title(app.UIAxes2, 'Modifying');
            xlabel(app.UIAxes2, 'Time [sec]');
            ylabel(app.UIAxes2, 'Q');
        end

        % Closeapp button pushed function
        function CloseappButtonPushed(app)
            delete(app);
        end

        % apply_cut button pushed function
        function apply_cutButtonPushed(app)
            app.tim(app.end2:end,:)=[];
            app.tim(1:app.beg2,:)=[];
            app.af(app.end2:end,:)=[];
            app.af(1:app.beg2,:)=[];
            app.f_inc(app.end2:end,:)=[];
            app.f_inc(1:app.beg2,:)=[];
            app.f_inc(:,1)=app.f_inc(:,1)-app.f_inc(1,1);
            %app.f_inc(:,2)=abs(app.f_inc(:,2));
            app.f_inc(:,2)=app.f_inc(:,2);
            m1=min(app.f_inc(:,2));
                    app.b_q=m1;
            app.tim(:)=app.tim(:)-app.tim(1);
            y2(:)=app.k_q.*app.tim(:)+mean(app.af(end-100:end));
            y2(:)=abs(y2(:)-(y2(end)-app.af(end)));
            app.UIAxes2.cla;
            hold(app.UIAxes2,'on');
            plot(app.UIAxes2,app.f_inc(:,1),app.f_inc(:,2));
            plot(app.UIAxes2,app.tim,app.af);
            plot(app.UIAxes2,app.tim,y2);
            title(app.UIAxes2, 'Modified');
            xlabel(app.UIAxes2, 'Time [sec]');
            ylabel(app.UIAxes2, 'Q');
        end

        % fitfuns value changed function
        function fitfunsValueChanged(app, event)
            val = app.fitfuns.Value;
            switch val
                case 'Fit1'
                    % fit
                    m1=min(app.f_inc(:,2));
                    %m1=0;
                    app.f_inc(:,2)=app.f_inc(:,2)-m1;
                    x0=[1 0 1000 0.1];
                    options = optimoptions(@fminunc,'MaxFunEvals',5000);
                    lb = [-400 -100 0 0];
                    ub = [400 100 400000 400];
                    [app.fit1,resnorm,residuals] = lsqcurvefit(@expfitfunUI,x0,app.f_inc(:,1)',app.f_inc(:,2)',lb,ub,options);
                    % text fittning parametrs
                    s1=num2str(app.fit1(1));
                    s2=num2str(app.fit1(2));
                    s3=num2str(app.fit1(3));
                    s4=num2str(app.fit1(4));
                    s5=num2str(app.tem);
                    s6=num2str(resnorm);
                    ss1=['A = ',s1];
                    ss2=['t_0 = ',s2];
                    ss3=['tau = ',s3];
                    ss4=['b = ',s4];
                    ss5=['Temperature = ',s5];
                    ss6=['Res Norm = ',s6];
                    app.printfit.Value = {ss1;ss2;ss3;ss4;ss5;ss6};
                    y=expfitfunUI(app.fit1,app.f_inc(:,1)');
                    app.UIAxes3.cla;
                    hold(app.UIAxes3,'on');
                    plot(app.UIAxes3,app.f_inc(:,1),app.f_inc(:,2));
                    plot(app.UIAxes3,app.f_inc(:,1),y);
                    title(app.UIAxes3, 'Fit1');
                    xlabel(app.UIAxes3, 'Time [sec]');
                    ylabel(app.UIAxes3, 'Q');
                    app.UIAxes4.cla;
                    hold(app.UIAxes4,'on');
                    plot(app.UIAxes4,app.f_inc(:,1),log(app.f_inc(:,2)));
                    plot(app.UIAxes4,app.f_inc(:,1),log(y));
                    title(app.UIAxes4, ' Log of Fit1');
                    xlabel(app.UIAxes4, 'Time [sec]');
                    ylabel(app.UIAxes4, 'ln(Q)');
                    app.UIAxes5.cla;
                    plot(app.UIAxes5,1:length(residuals),residuals);
                    [a,bc,rod,h,p]=errans(y,app.f_inc(:,2)',0,residuals,app.nbin.Value);
                   
                    %hold(app.UIAxes6,'on');
                    
                    plot(app.UIAxes6,bc,rod,'-xk');
                    pd = fitdist(residuals','Normal');
                    fgaus = fit(bc,rod,'gauss1');
                    y_gaus(:)=fgaus.a1.*exp(-((bc(:)-fgaus.b1)/fgaus.c1).^2);
                    plot(app.UIAxes6,bc,y_gaus,'--k');
                    sa6=num2str(fgaus.a1);
                    sa7=num2str(fgaus.b1);
                    sa8=num2str(fgaus.c1);
                    saa6=['A_gaus = ',sa6];
                    saa7=['mu_gaus = ',sa7];
                    saa8=['signma_gaus = ',sa8];
                    sa1=num2str(a);
                    sa2=num2str(h);
                    sa3=num2str(p);
                    sa4=num2str(pd.mu);
                    sa5=num2str(pd.sigma);
                    saa1=['chi2 = ',sa1];
                    saa2=['h = ',sa2];
                    saa3=['p = ',sa3];
                    saa4=['mu = ',sa4];
                    saa5=['sigma = ',sa5];
                    app.errtext.Value = {ss6;saa1;saa2;saa3;saa4;saa5;saa6;saa7;saa8};
                    [yfit,xfit,nb,ne]=remodel(y,m1,app.tim,app.sfluid.Value,app.beg2,app.ff,app.poly_off,app.c(app.n_pul,1),app.c(app.n_pul,2),app.beg1);
                    app.UIAxes7.cla;
                    hold(app.UIAxes7,'on');
                    plot(app.UIAxes7,app.ff(nb:ne,1),app.ff(nb:ne,5),'-b');
                    plot(app.UIAxes7,xfit,yfit,'r');
                case 'None'
                    ly(:)=log(app.f_inc(:,2));
                    app.UIAxes3.cla;
                    plot(app.UIAxes3,app.f_inc(:,1),app.f_inc(:,2));                    
                    title(app.UIAxes3, 'None');
                    xlabel(app.UIAxes3, 'Time [sec]');
                    ylabel(app.UIAxes3, 'Q');
                    % log
                    app.UIAxes4.cla;
                    plot(app.UIAxes4,app.f_inc(:,1),ly(:));                    
                    title(app.UIAxes4, 'Log');
                    xlabel(app.UIAxes4, 'Time [sec]');
                    ylabel(app.UIAxes4, 'ln(Q)');
                     
                case 'Fit2'
                    m1=min(app.f_inc(:,2));
                    %m1=0;
                    app.f_inc(:,2)=app.f_inc(:,2)-m1;
                    x0=[app.b_q app.k_q app.fit1(3) 0];
                    if app.polydegree.Value == 1
                        lb = [-10 -app.k_f2.Value 0 -40000];
                        ub = [10 app.k_f2.Value 400000 40000];
                    else
                        lb = [-10 -app.k_f2.Value 0 -40000];
                        ub = [10 app.k_f2.Value 400000 40000];
                    end
                    options = optimoptions(@fminunc,'MaxFunEvals',5000);
%                     lb = [-4000 0 0 -40000];
%                     ub = [4000 1 400000 40000];
                    [app.fit2,resnorm,residuals] = lsqcurvefit(@exp2fitfunUI,x0,app.tim',app.af',lb,ub,options);
                    y1(:)=app.af(:)-app.fit2(2).*app.tim(:)-app.fit2(1); % substract
                    y2(:)=app.fit2(2).*app.tim(:)+app.fit2(1);
                    s1=num2str(app.fit2(1));
                    s2=num2str(app.fit2(2));
                    s3=num2str(app.fit2(3));
                   s4=num2str(app.fit2(4));
                    s5=num2str(app.tem);
                    s6=num2str(resnorm);
                    ss1=['A = ',s1];
                    ss2=['B = ',s2];
                    ss3=['tau = ',s3];
                    ss4=['t0 = ',s4];
                    ss5=['Temperature = ',s5];
                    ss6=['Res Norm = ',s6];
                    app.printfit.Value = {ss1;ss2;ss3;ss4;ss5;ss6};
                    y=exp2fitfunUI(app.fit2,app.tim');
                    app.UIAxes3.cla;
                    hold(app.UIAxes3,'on');
                    plot(app.UIAxes3,app.tim,app.af);
                    plot(app.UIAxes3,app.tim,y);
                    plot(app.UIAxes3,app.tim,y2);
                    title(app.UIAxes3, 'Fit2');
                    xlabel(app.UIAxes3, 'Time [sec]');
                    ylabel(app.UIAxes3, 'Q');
                    
                    app.UIAxes4.cla;
                    hold(app.UIAxes4,'on');
                    plot(app.UIAxes4,app.f_inc(:,1),app.f_inc(:,2));
                    plot(app.UIAxes4,app.tim,y1);
                    title(app.UIAxes4, ' Log of Fit1');
                    xlabel(app.UIAxes4, 'Time [sec]');
                    ylabel(app.UIAxes4, 'Q');
                    app.UIAxes5.cla;
                    plot(app.UIAxes5,1:length(residuals),residuals);
                    
                    [a,bc,rod,h,p]=errans(y,app.af',5,residuals,app.nbin.Value);
                    
                    %hold(app.UIAxes6,'on');
                    plot(app.UIAxes6,bc,rod,'-xr');
                    pd = fitdist(residuals','Normal');
                    fgaus = fit(bc,rod,'gauss1');
                    y_gaus(:)=fgaus.a1.*exp(-((bc(:)-fgaus.b1)/fgaus.c1).^2);
                    plot(app.UIAxes6,bc,y_gaus,'--r');
                    sa6=num2str(fgaus.a1);
                    sa7=num2str(fgaus.b1);
                    sa8=num2str(fgaus.c1);
                    saa6=['A_gaus = ',sa6];
                    saa7=['mu_gaus = ',sa7];
                    saa8=['signma_gaus = ',sa8];
                    sa1=num2str(a);
                    sa2=num2str(h);
                    sa3=num2str(p);
                    sa4=num2str(pd.mu);
                    sa5=num2str(pd.sigma);
                    saa1=['chi2 = ',sa1];
                    saa2=['h = ',sa2];
                    saa3=['p = ',sa3];
                    saa4=['mu = ',sa4];
                    saa5=['sigma = ',sa5];
                    app.errtext.Value = {ss6;saa1;saa2;saa3;saa4;saa5;saa6;saa7;saa8};
                    [yfit,xfit,nb,ne]=remodel(y,m1,app.tim,app.sfluid.Value,app.beg2,app.ff,app.poly_off,app.c(app.n_pul,1),app.c(app.n_pul,2),app.beg1);
                    app.UIAxes7.cla;
                    hold(app.UIAxes7,'on');
                    plot(app.UIAxes7,app.ff(nb:ne,1),app.ff(nb:ne,5),'-b');
                    plot(app.UIAxes7,xfit,yfit,'r');
                case 'Fit3'
                    m1=min(app.f_inc(:,2));
                    app.f_inc(:,2)=app.f_inc(:,2)-m1;
                    x0=[app.fit2(1) app.fit2(2) app.fit2(3) app.fit2(4) 0];
                    options = optimoptions(@fminunc,'MaxFunEvals',5000);
                    lb=[0 -1000 0 -40000 0];
                    ub=[40000 1000 40000 40000 1000];
                    [app.fit3,resnorm,residuals] = lsqcurvefit(@exp3fitfunUI,x0,app.tim',app.af',lb,ub,options);
                    y1(:)=app.af(:)-app.fit3(2).*app.tim(:)-app.fit3(1); % substract
                    y2(:)=app.fit3(2).*app.tim(:)+app.fit3(1);
                    taue=app.fit3(3)+app.fit3(5)*app.tim(end);
                    s1=num2str(app.fit3(1));
                    s2=num2str(app.fit3(2));
                    s3=num2str(app.fit3(3));
                   s4=num2str(app.fit3(4));
                    s5=num2str(app.tem);
                    s6=num2str(resnorm);
                    s7=num2str(app.fit3(5));
                    s8=num2str(taue);
                    ss1=['A = ',s1];
                    ss2=['B = ',s2];
                    ss3=['tau = ',s3];
                    ss4=['t0 = ',s4];
                    ss7=['C = ',s7];
                    ss8=['tau_e = ',s8];
                    ss5=['Temperature = ',s5];
                    ss6=['Res Norm = ',s6];
                    app.printfit.Value = {ss1;ss2;ss3;ss4;ss7;ss8;ss5;ss6};
                    y=exp3fitfunUI(app.fit3,app.tim');
                    app.UIAxes3.cla;
                    hold(app.UIAxes3,'on');
                    plot(app.UIAxes3,app.tim,app.af);
                    plot(app.UIAxes3,app.tim,y);
                    plot(app.UIAxes3,app.tim,y2);
                    title(app.UIAxes3, 'Fit3');
                    xlabel(app.UIAxes3, 'Time [sec]');
                    ylabel(app.UIAxes3, 'Q');
                    
                    app.UIAxes4.cla;
                    hold(app.UIAxes4,'on');
                    plot(app.UIAxes4,app.f_inc(:,1),app.f_inc(:,2));
                    plot(app.UIAxes4,app.tim,y1);
                    title(app.UIAxes4, ' Log of Fit1');
                    xlabel(app.UIAxes4, 'Time [sec]');
                    ylabel(app.UIAxes4, 'Q');
                    app.UIAxes5.cla;
                    plot(app.UIAxes5,1:length(residuals),residuals);
                    
                    [a,bc,rod,h,p]=errans(y,app.af',5,residuals,app.nbin.Value);
                    
                    %hold(app.UIAxes6,'on');
                    plot(app.UIAxes6,bc,rod,'-xg','LineWidth',1);
                    pd = fitdist(residuals','Normal');
                    fgaus = fit(bc,rod,'gauss1');
                    y_gaus(:)=fgaus.a1.*exp(-((bc(:)-fgaus.b1)/fgaus.c1).^2);
                    plot(app.UIAxes6,bc,y_gaus,'--g');
                    sa6=num2str(fgaus.a1);
                    sa7=num2str(fgaus.b1);
                    sa8=num2str(fgaus.c1);
                    saa6=['A_gaus = ',sa6];
                    saa7=['mu_gaus = ',sa7];
                    saa8=['signma_gaus = ',sa8];
                    sa1=num2str(a);
                    sa2=num2str(h);
                    sa3=num2str(p);
                    sa4=num2str(pd.mu);
                    sa5=num2str(pd.sigma);
                    saa1=['chi2 = ',sa1];
                    saa2=['h = ',sa2];
                    saa3=['p = ',sa3];
                    saa4=['mu = ',sa4];
                    saa5=['sigma = ',sa5];
                    app.errtext.Value = {ss6;saa1;saa2;saa3;saa4;saa5;saa6;saa7;saa8};
                    [yfit,xfit,nb,ne]=remodel(y,m1,app.tim,app.sfluid.Value,app.beg2,app.ff,app.poly_off,app.c(app.n_pul,1),app.c(app.n_pul,2),app.beg1);
                    app.UIAxes7.cla;
                    hold(app.UIAxes7,'on');
                    plot(app.UIAxes7,app.ff(nb:ne,1),app.ff(nb:ne,5),'-b');
                    plot(app.UIAxes7,xfit,yfit,'r');
                    
                case 'Fit4'
                    m1=min(app.f_inc(:,2));
                    app.f_inc(:,2)=app.f_inc(:,2)-m1;
                    x0=[app.fit2(1) app.fit2(2) app.fit2(3) app.fit2(4) 0];
                    options = optimoptions(@fminunc,'MaxFunEvals',5000);
                    lb = [0 0 0 -40000 0];
                    ub = [4000 1 400000 40000 1000];
                    [app.fit4,resnorm,residuals] = lsqcurvefit(@exp4fitfunUI,x0,app.tim',app.af',lb,ub,options);
                    y1(:)=app.af(:)-app.fit4(2).*app.tim(:)-app.fit4(1); % substract
                    y2(:)=app.fit4(2).*app.tim(:)+app.fit4(1);
                    %taue=app.fit4(3)+app.fit4(5)*app.tim(end);
                    s1=num2str(app.fit4(1));
                    s2=num2str(app.fit4(2));
                    s3=num2str(app.fit4(3));
                   s4=num2str(app.fit4(4));
                    s5=num2str(app.tem);
                    s6=num2str(resnorm);
                    s7=num2str(app.fit4(5));
                    %s8=num2str(taue);
                    ss1=['A = ',s1];
                    ss2=['B = ',s2];
                    ss3=['tau = ',s3];
                    ss4=['t0 = ',s4];
                    ss7=['C = ',s7];
                    %ss8=['tau_e = ',s8];
                    ss5=['Temperature = ',s5];
                    ss6=['Res Norm = ',s6];
                    app.printfit.Value = {ss1;ss2;ss3;ss4;ss7;ss5;ss6};
                    y=exp4fitfunUI(app.fit4,app.tim');
                    app.UIAxes3.cla;
                    hold(app.UIAxes3,'on');
                    plot(app.UIAxes3,app.tim,app.af);
                    plot(app.UIAxes3,app.tim,y);
                    plot(app.UIAxes3,app.tim,y2);
                    title(app.UIAxes3, 'Fit4');
                    xlabel(app.UIAxes3, 'Time [sec]');
                    ylabel(app.UIAxes3, 'Q');
                    
                    app.UIAxes4.cla;
                    hold(app.UIAxes4,'on');
                    plot(app.UIAxes4,app.f_inc(:,1),app.f_inc(:,2));
                    plot(app.UIAxes4,app.tim,y1);
                    title(app.UIAxes4, ' Log of Fit1');
                    xlabel(app.UIAxes4, 'Time [sec]');
                    ylabel(app.UIAxes4, 'Q');
                    app.UIAxes5.cla;
                    plot(app.UIAxes5,1:length(residuals),residuals);
            end
        end

        % qends button pushed function
        function qendsButtonPushed(app)
            [c1,~]=size(app.c);
            y=zeros(c1,2);
            num_off=90;
            for ii=1:c1
                y(ii,1)=mean(app.ff1(app.c(ii,2)-num_off-50:app.c(ii,2)-num_off,1));
                y(ii,2)=mean(app.ff1(app.c(ii,2)-num_off-50:app.c(ii,2)-num_off,5));
            end
            app.poly_off=polyfit(y(:,1),y(:,2),app.polydegree.Value);
            %y1 = polyval(app.poly_off,y(:,1));
            y1 = polyval(app.poly_off,app.ff1(:,1));
            s1=num2str(app.poly_off(1));
            app.k_q=app.poly_off(1);
            st=['K_q = ',s1];
            app.textfork.Value = {st};
            app.UIAxes1.cla;
            hold(app.UIAxes1,'on');
            plot(app.UIAxes1,y(:,1),y(:,2),'-o');
            plot(app.UIAxes1,app.ff1(:,1),y1,'-');
            title(app.UIAxes1, ' Time dependence of Q at the end of the pulses');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Q');
           % ylim(app.UIAxes1,[-1 10]);
        end

        % clear_err button pushed function
        function clear_errButtonPushed(app)
            app.UIAxes6.cla;
        end

        % applypoly button pushed function
        function applypolyButtonPushed(app)
            y1 = polyval(app.poly_off,app.ff1(:,1));
            if  app.sfluid.Value == 0
                 app.ff1(:,5)=app.ff1(:,5)-y1(:);
            else
                app.ff1(:,5)=y1(:)-app.ff1(:,5);
            end
                 plot(app.UIAxes1,app.ff1(:,1),app.ff1(:,5),'--');
                 ylim(app.UIAxes1,[-1 3]);
            if app.polydegree.Value ~= 1
                 
                 app.k_q=0;
                 app.b_q=0;
                 
            else
                
                 app.k_q=y1(1);
                 app.b_q=y1(2);
                 
            end
        end

        % OpenFile button pushed function
        function OpenFileButtonPushed(app)
            runStartupFcn(app, @startupFcn);
        end

        % tempSF button pushed function
        function tempSFButtonPushed(app)
            yf=testL(app.ff1);
            
            app.UIAxes1.cla;
            hold(app.UIAxes1,'on');
            plot(app.UIAxes1,app.ff1(:,1),app.ff1(:,12),'.b');
            plot(app.UIAxes1,app.ff1(:,1),yf,'-r');
            title(app.UIAxes1, 'New temperature');
            xlabel(app.UIAxes1, 'Time [sec]');
            ylabel(app.UIAxes1, 'Temperature [mK]');
            ylim(app.UIAxes1,'auto');
            app.ff1(:,12)=yf(:);
        end

        % exportF button pushed function
        function exportFButtonPushed(app)
            y=exp2fitfunUI(app.fit2,app.tim');
            
            %        plot(app.UIAxes3,app.tim,app.af);
            %        plot(app.UIAxes3,app.tim,y);
            fileID = fopen('fitAndSub.dat','w');
            for ii=1:length(y)
                fprintf(fileID,'%f \t %f \t %f \n',app.tim(ii),app.af(ii),y(ii));
            end
            
            
            fclose(fileID);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 640 482];
            app.UIFigure.Name = 'UI Figure';
            setAutoResize(app, app.UIFigure, true)

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Units = 'pixels';
            app.TabGroup.Position = [0 -1 640 483];

            % Create prep
            app.prep = uitab(app.TabGroup);
            app.prep.Units = 'pixels';
            app.prep.Title = 'Preparing';

            % Create UIAxes1
            app.UIAxes1 = uiaxes(app.prep);
            title(app.UIAxes1, 'Title');
            xlabel(app.UIAxes1, 'X');
            ylabel(app.UIAxes1, 'Y');
            app.UIAxes1.Box = 'on';
            app.UIAxes1.XGrid = 'on';
            app.UIAxes1.YGrid = 'on';
            app.UIAxes1.Position = [277 155 347 292];

            % Create init
            app.init = uipanel(app.prep);
            app.init.BorderType = 'line';
            app.init.TitlePosition = 'centertop';
            app.init.Title = 'Initial data manipulation';
            app.init.FontName = 'Helvetica';
            app.init.FontUnits = 'pixels';
            app.init.FontSize = 16;
            app.init.Units = 'pixels';
            app.init.Position = [0 30 260 431];

            % Create LabelSlider
            app.LabelSlider = uilabel(app.init);
            app.LabelSlider.HorizontalAlignment = 'right';
            app.LabelSlider.FontSize = 14;
            app.LabelSlider.Position = [98 372 36 18];
            app.LabelSlider.Text = 'Begin';

            % Create cutBegin
            app.cutBegin = uislider(app.init);
            app.cutBegin.Limits = [1 100];
            app.cutBegin.ValueChangedFcn = createCallbackFcn(app, @cutBeginValueChanged);
            app.cutBegin.Position = [37 362 150 3];
            app.cutBegin.Value = 1;

            % Create Label
            app.Label = uilabel(app.init);
            app.Label.HorizontalAlignment = 'right';
            app.Label.FontSize = 14;
            app.Label.Position = [109 295 25 18];
            app.Label.Text = 'End';

            % Create cutEnd
            app.cutEnd = uislider(app.init);
            app.cutEnd.ValueChangedFcn = createCallbackFcn(app, @cutEndValueChanged);
            app.cutEnd.Position = [37 285 150 3];

            % Create setP
            app.setP = uibutton(app.init, 'push');
            app.setP.ButtonPushedFcn = createCallbackFcn(app, @setPButtonPushed);
            app.setP.Position = [31 212 100 22];
            app.setP.Text = 'Set';

            % Create findk
            app.findk = uibutton(app.init, 'push');
            app.findk.ButtonPushedFcn = createCallbackFcn(app, @findkButtonPushed);
            app.findk.Position = [31 157 100 22];
            app.findk.Text = 'Find K''s';

            % Create findpulse
            app.findpulse = uibutton(app.init, 'push');
            app.findpulse.ButtonPushedFcn = createCallbackFcn(app, @findpulseButtonPushed);
            app.findpulse.Position = [31 108 100 22];
            app.findpulse.Text = 'Find pulses';

            % Create qvstime
            app.qvstime = uibutton(app.init, 'push');
            app.qvstime.ButtonPushedFcn = createCallbackFcn(app, @qvstimeButtonPushed);
            app.qvstime.Position = [148 212 100 22];
            app.qvstime.Text = 'SeeQ';

            % Create pulsenumtext
            app.pulsenumtext = uitextarea(app.init);
            app.pulsenumtext.FontSize = 14;
            app.pulsenumtext.Position = [148 109 100 26];

            % Create qends
            app.qends = uibutton(app.init, 'push');
            app.qends.ButtonPushedFcn = createCallbackFcn(app, @qendsButtonPushed);
            app.qends.Position = [22 56 100 22];
            app.qends.Text = 'Qends';

            % Create LabelNumericEditField4
            app.LabelNumericEditField4 = uilabel(app.init);
            app.LabelNumericEditField4.HorizontalAlignment = 'right';
            app.LabelNumericEditField4.Position = [158 83 66 15];
            app.LabelNumericEditField4.Text = 'Poly degree';

            % Create polydegree
            app.polydegree = uieditfield(app.init, 'numeric');
            app.polydegree.Limits = [1 10];
            app.polydegree.Position = [148 56 100 22];
            app.polydegree.Value = 1;

            % Create applypoly
            app.applypoly = uibutton(app.init, 'push');
            app.applypoly.ButtonPushedFcn = createCallbackFcn(app, @applypolyButtonPushed);
            app.applypoly.Position = [22 15 100 22];
            app.applypoly.Text = 'Apply poly';

            % Create OpenFile
            app.OpenFile = uibutton(app.init, 'push');
            app.OpenFile.ButtonPushedFcn = createCallbackFcn(app, @OpenFileButtonPushed);
            app.OpenFile.Position = [148 15 100 22];
            app.OpenFile.Text = 'open file';

            % Create LabelNumericEditField6
            app.LabelNumericEditField6 = uilabel(app.init);
            app.LabelNumericEditField6.HorizontalAlignment = 'right';
            app.LabelNumericEditField6.Position = [158 178 66 15];
            app.LabelNumericEditField6.Text = 'SFluid true?';

            % Create sfluid
            app.sfluid = uieditfield(app.init, 'numeric');
            app.sfluid.Limits = [0 1];
            app.sfluid.ValueDisplayFormat = '%.0f';
            app.sfluid.Position = [157 156 100 22];

            % Create LabelTextArea
            app.LabelTextArea = uilabel(app.prep);
            app.LabelTextArea.HorizontalAlignment = 'right';
            app.LabelTextArea.FontSize = 14;
            app.LabelTextArea.Position = [403 128 54 18];
            app.LabelTextArea.Text = 'K values';

            % Create textfork
            app.textfork = uitextarea(app.prep);
            app.textfork.FontSize = 14;
            app.textfork.Position = [398 81 134 43];
            app.textfork.Value = {'K_T = 0'; 'K_Q = 0'};

            % Create LabelEditField
            app.LabelEditField = uilabel(app.prep);
            app.LabelEditField.HorizontalAlignment = 'right';
            app.LabelEditField.Position = [426 50 62 15];
            app.LabelEditField.Text = 'Folder path';

            % Create path
            app.path = uieditfield(app.prep, 'text');
            app.path.Position = [412 22 192 20];
            app.path.Value = 'D:\therm_transport\';

            % Create tempSF
            app.tempSF = uibutton(app.prep, 'push');
            app.tempSF.ButtonPushedFcn = createCallbackFcn(app, @tempSFButtonPushed);
            app.tempSF.Position = [264.5 102 109 22];
            app.tempSF.Text = 'Temperature at SO';

            % Create Tab2
            app.Tab2 = uitab(app.TabGroup);
            app.Tab2.Units = 'pixels';
            app.Tab2.Title = 'Filtering';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.Tab2);
            title(app.UIAxes2, 'Title');
            xlabel(app.UIAxes2, 'X');
            ylabel(app.UIAxes2, 'Y');
            app.UIAxes2.Box = 'on';
            app.UIAxes2.Position = [199 119 440 342];

            % Create Panel
            app.Panel = uipanel(app.Tab2);
            app.Panel.BorderType = 'line';
            app.Panel.TitlePosition = 'centertop';
            app.Panel.Title = 'Pulse selection';
            app.Panel.FontName = 'Helvetica';
            app.Panel.FontUnits = 'pixels';
            app.Panel.FontSize = 12;
            app.Panel.Units = 'pixels';
            app.Panel.Position = [0 73 200 388];

            % Create LabelNumericEditField
            app.LabelNumericEditField = uilabel(app.Panel);
            app.LabelNumericEditField.HorizontalAlignment = 'center';
            app.LabelNumericEditField.VerticalAlignment = 'center';
            app.LabelNumericEditField.Position = [3 318 76 15];
            app.LabelNumericEditField.Text = 'Pulse number';

            % Create pulsenum
            app.pulsenum = uieditfield(app.Panel, 'numeric');
            app.pulsenum.ValueChangedFcn = createCallbackFcn(app, @pulsenumValueChanged);
            app.pulsenum.Limits = [0 Inf];
            app.pulsenum.ValueDisplayFormat = '%.0f';
            app.pulsenum.Position = [87 314 100 22];

            % Create LabelDropDown
            app.LabelDropDown = uilabel(app.Panel);
            app.LabelDropDown.HorizontalAlignment = 'right';
            app.LabelDropDown.Position = [39 284 33 15];
            app.LabelDropDown.Text = 'Filters';

            % Create filt
            app.filt = uidropdown(app.Panel);
            app.filt.Items = {'Median', 'Moving', 'None'};
            app.filt.ValueChangedFcn = createCallbackFcn(app, @filtValueChanged, true);
            app.filt.Position = [87 280 100 22];
            app.filt.Value = 'Median';

            % Create LabelNumericEditField2
            app.LabelNumericEditField2 = uilabel(app.Panel);
            app.LabelNumericEditField2.HorizontalAlignment = 'right';
            app.LabelNumericEditField2.Position = [13 248 57 15];
            app.LabelNumericEditField2.Text = '1st field filt';

            % Create filt_n1
            app.filt_n1 = uieditfield(app.Panel, 'numeric');
            app.filt_n1.ValueChangedFcn = createCallbackFcn(app, @filt_n1ValueChanged);
            app.filt_n1.Limits = [0 Inf];
            app.filt_n1.Position = [85 244 100 22];
            app.filt_n1.Value = 1;

            % Create apply1
            app.apply1 = uibutton(app.Panel, 'push');
            app.apply1.ButtonPushedFcn = createCallbackFcn(app, @apply1ButtonPushed);
            app.apply1.Position = [99 205 86 22];
            app.apply1.Text = 'Apply filter';

            % Create incline
            app.incline = uibutton(app.Panel, 'push');
            app.incline.ButtonPushedFcn = createCallbackFcn(app, @inclineButtonPushed);
            app.incline.Position = [3 205 81 22];
            app.incline.Text = 'Incline';

            % Create LabelSlider2
            app.LabelSlider2 = uilabel(app.Panel);
            app.LabelSlider2.HorizontalAlignment = 'right';
            app.LabelSlider2.Position = [77 175 32 15];
            app.LabelSlider2.Text = 'Begin';

            % Create cutBegin2
            app.cutBegin2 = uislider(app.Panel);
            app.cutBegin2.Limits = [1 100];
            app.cutBegin2.ValueChangedFcn = createCallbackFcn(app, @cutBegin2ValueChanged);
            app.cutBegin2.Position = [18 167 150 3];
            app.cutBegin2.Value = 1;

            % Create Label2
            app.Label2 = uilabel(app.Panel);
            app.Label2.HorizontalAlignment = 'right';
            app.Label2.Position = [87 105 22 15];
            app.Label2.Text = 'End';

            % Create cutEnd2
            app.cutEnd2 = uislider(app.Panel);
            app.cutEnd2.Limits = [1 100];
            app.cutEnd2.ValueChangedFcn = createCallbackFcn(app, @cutEnd2ValueChanged);
            app.cutEnd2.Position = [18 97 150 3];
            app.cutEnd2.Value = 1;

            % Create apply_cut
            app.apply_cut = uibutton(app.Panel, 'push');
            app.apply_cut.ButtonPushedFcn = createCallbackFcn(app, @apply_cutButtonPushed);
            app.apply_cut.Position = [85 24 100 22];
            app.apply_cut.Text = 'Apply cut';

            % Create Tab3
            app.Tab3 = uitab(app.TabGroup);
            app.Tab3.Units = 'pixels';
            app.Tab3.Title = 'Fitting';

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.Tab3);
            title(app.UIAxes3, 'Title');
            xlabel(app.UIAxes3, 'X');
            ylabel(app.UIAxes3, 'Y');
            app.UIAxes3.Box = 'on';
            app.UIAxes3.XGrid = 'on';
            app.UIAxes3.YGrid = 'on';
            app.UIAxes3.Position = [252 230 386 231];

            % Create Panel2
            app.Panel2 = uipanel(app.Tab3);
            app.Panel2.BorderType = 'line';
            app.Panel2.TitlePosition = 'centertop';
            app.Panel2.Title = 'Filter tools';
            app.Panel2.FontName = 'Helvetica';
            app.Panel2.FontUnits = 'pixels';
            app.Panel2.FontSize = 12;
            app.Panel2.Units = 'pixels';
            app.Panel2.Position = [0 240 236 221];

            % Create LabelDropDown2
            app.LabelDropDown2 = uilabel(app.Panel2);
            app.LabelDropDown2.HorizontalAlignment = 'right';
            app.LabelDropDown2.Position = [78 170 70 15];
            app.LabelDropDown2.Text = 'Fitting model';

            % Create fitfuns
            app.fitfuns = uidropdown(app.Panel2);
            app.fitfuns.Items = {'Fit1', 'Fit2', 'Fit3', 'Fit4', 'None'};
            app.fitfuns.ValueChangedFcn = createCallbackFcn(app, @fitfunsValueChanged, true);
            app.fitfuns.Position = [67 135 100 22];
            app.fitfuns.Value = 'Fit1';

            % Create Label3
            app.Label3 = uilabel(app.Panel2);
            app.Label3.Position = [37 107 135 19];
            app.Label3.Text = 'Fit 1: A*exp(-t+t0/tau) + B';

            % Create Label4
            app.Label4 = uilabel(app.Panel2);
            app.Label4.Position = [37 87 150 15];
            app.Label4.Text = 'Fit 2: A + b*t + exp(-t+t0/tau)';

            % Create Label5
            app.Label5 = uilabel(app.Panel2);
            app.Label5.Position = [25 65 182 15];
            app.Label5.Text = 'Fit 3: A + b*t + exp(-t+t0/(tau+C*t))';

            % Create exportF
            app.exportF = uibutton(app.Panel2, 'push');
            app.exportF.ButtonPushedFcn = createCallbackFcn(app, @exportFButtonPushed);
            app.exportF.Position = [63 22 100 22];
            app.exportF.Text = 'Export graphs';

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.Tab3);
            title(app.UIAxes4, 'Title');
            xlabel(app.UIAxes4, 'X');
            ylabel(app.UIAxes4, 'Y');
            app.UIAxes4.Box = 'on';
            app.UIAxes4.XGrid = 'on';
            app.UIAxes4.YGrid = 'on';
            app.UIAxes4.Position = [252 6 386 225];

            % Create Panel3
            app.Panel3 = uipanel(app.Tab3);
            app.Panel3.BorderType = 'line';
            app.Panel3.TitlePosition = 'centertop';
            app.Panel3.Title = 'Fitting parameters';
            app.Panel3.FontName = 'Helvetica';
            app.Panel3.FontUnits = 'pixels';
            app.Panel3.FontSize = 12;
            app.Panel3.Units = 'pixels';
            app.Panel3.Position = [0 10 236 221];

            % Create printfit
            app.printfit = uitextarea(app.Panel3);
            app.printfit.Position = [42 70 150 119];
            app.printfit.Value = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'};

            % Create Tab4
            app.Tab4 = uitab(app.TabGroup);
            app.Tab4.Units = 'pixels';
            app.Tab4.Title = 'Residuals';

            % Create UIAxes5
            app.UIAxes5 = uiaxes(app.Tab4);
            title(app.UIAxes5, 'Residuals');
            xlabel(app.UIAxes5, 'Time');
            ylabel(app.UIAxes5, 'Deviation');
            app.UIAxes5.Box = 'on';
            app.UIAxes5.XGrid = 'on';
            app.UIAxes5.YGrid = 'on';
            app.UIAxes5.Position = [190 236 449 225];

            % Create UIAxes6
            app.UIAxes6 = uiaxes(app.Tab4);
            title(app.UIAxes6, 'Title');
            xlabel(app.UIAxes6, 'X');
            ylabel(app.UIAxes6, 'Y');
            app.UIAxes6.Box = 'on';
            app.UIAxes6.NextPlot = 'add';
            app.UIAxes6.XGrid = 'on';
            app.UIAxes6.YGrid = 'on';
            app.UIAxes6.Position = [190 3 448 234];

            % Create errtext
            app.errtext = uitextarea(app.Tab4);
            app.errtext.Position = [15 236 150 209];
            app.errtext.Value = {'res_n'; 'chi'};

            % Create LabelNumericEditField3
            app.LabelNumericEditField3 = uilabel(app.Tab4);
            app.LabelNumericEditField3.HorizontalAlignment = 'right';
            app.LabelNumericEditField3.Position = [26 128 24 15];
            app.LabelNumericEditField3.Text = 'Bins';

            % Create nbin
            app.nbin = uieditfield(app.Tab4, 'numeric');
            app.nbin.Limits = [1 Inf];
            app.nbin.Position = [65 124 100 22];
            app.nbin.Value = 21;

            % Create clear_err
            app.clear_err = uibutton(app.Tab4, 'push');
            app.clear_err.ButtonPushedFcn = createCallbackFcn(app, @clear_errButtonPushed);
            app.clear_err.Position = [57 78 100 22];
            app.clear_err.Text = 'Clear';

            % Create Closeapp
            app.Closeapp = uibutton(app.Tab4, 'push');
            app.Closeapp.ButtonPushedFcn = createCallbackFcn(app, @CloseappButtonPushed);
            app.Closeapp.Position = [57 24 100 22];
            app.Closeapp.Text = 'Close';

            % Create LabelNumericEditField5
            app.LabelNumericEditField5 = uilabel(app.Tab4);
            app.LabelNumericEditField5.HorizontalAlignment = 'right';
            app.LabelNumericEditField5.Position = [74 195 52 15];
            app.LabelNumericEditField5.Text = 'Edit Field';

            % Create k_f2
            app.k_f2 = uieditfield(app.Tab4, 'numeric');
            app.k_f2.Position = [57 160 100 22];
            app.k_f2.Value = 0.0002;

            % Create Tab5
            app.Tab5 = uitab(app.TabGroup);
            app.Tab5.Units = 'pixels';
            app.Tab5.Title = 'Re-model';

            % Create UIAxes7
            app.UIAxes7 = uiaxes(app.Tab5);
            title(app.UIAxes7, 'Re-model');
            xlabel(app.UIAxes7, 'Original time');
            ylabel(app.UIAxes7, 'Original Q');
            app.UIAxes7.Box = 'on';
            app.UIAxes7.XGrid = 'on';
            app.UIAxes7.YGrid = 'on';
            app.UIAxes7.Position = [0 2 639 459];
        end
    end

    methods (Access = public)

        % Construct app
        function app = analysisAD()

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
