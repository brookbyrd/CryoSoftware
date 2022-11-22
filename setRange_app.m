classdef setRange_app < matlab.apps.AppBase
    
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        ChannelPanel_760       matlab.ui.container.Panel
        ButtonGroup_4        matlab.ui.container.ButtonGroup
        ICGButton            matlab.ui.control.RadioButton
        IR800Button          matlab.ui.control.RadioButton
        max530       matlab.ui.control.NumericEditField
        max530Label  matlab.ui.control.Label
        ChannelPanel_635       matlab.ui.container.Panel
        ButtonGroup_3        matlab.ui.container.ButtonGroup
        AF680Button          matlab.ui.control.RadioButton
        OF550Button          matlab.ui.control.RadioButton   
        RhdButton          matlab.ui.control.RadioButton
        OF650Button          matlab.ui.control.RadioButton
        IR680Button          matlab.ui.control.RadioButton
        AF647Button          matlab.ui.control.RadioButton
        PPIXButton          matlab.ui.control.RadioButton
        min760       matlab.ui.control.NumericEditField
        min760Label  matlab.ui.control.Label
        max635       matlab.ui.control.NumericEditField
        max635Label  matlab.ui.control.Label
        min635       matlab.ui.control.NumericEditField
        min635Label  matlab.ui.control.Label
        ChannelPanel_530       matlab.ui.container.Panel
        ButtonGroup_2        matlab.ui.container.ButtonGroup
        tdTomatoButton       matlab.ui.control.RadioButton
        max470       matlab.ui.control.NumericEditField
        max760 matlab.ui.control.NumericEditField
        ChannelPanel         matlab.ui.container.Panel
        ButtonGroup          matlab.ui.container.ButtonGroup
        AF488Button          matlab.ui.control.RadioButton
        BodipyButton         matlab.ui.control.RadioButton
        FITCButton           matlab.ui.control.RadioButton
        min530       matlab.ui.control.NumericEditField
        min530Label  matlab.ui.control.Label
        max470Label    matlab.ui.control.Label
        max760Label    matlab.ui.control.Label
        min470       matlab.ui.control.NumericEditField
        minLabel             matlab.ui.control.Label
        agent470             string
        agent530             string
        agent635             string
        agent760             string
        acceptButton         matlab.ui.control.Button
        running              logical
        range                double
        minWL                double
        maxWL                double
        channels             string
        agents               string
        
    end
    
    
    
    
    % Callbacks that handle component events
    methods (Access = private)
        
        % Callback function: ButtonGroup, ButtonGroup_3,
        % max635, min635
        
        % Function to stop running the app once the user accepts the
        % selection
        function AcceptRange(app, event)
            disp('Ranges accepted');
            
            % set the range to the given ranges
            for c = 1:length(app.channels)
                switch(app.channels(c))
                    case "470"
                        app.range(c,1) = app.min470.Value;
                        app.range(c,2) = app.max470.Value;
                        app.agents(c) = app.agent470;
                    case "530"
                        app.range(c,1) = app.min530.Value;
                        app.range(c,2) = app.max530.Value;
                        app.agents(c) = app.agent530;
                    case "635"
                        app.range(c,1) = app.min635.Value;
                        app.range(c,2) = app.max635.Value;
                        app.agents(c) = app.agent635;
                    case "760"
                        app.range(c,1) = app.min760.Value;
                        app.range(c,2) = app.max760.Value;
                        app.agents(c) = app.agent760;
                    case "WL"
                        app.range(c,1) = app.minWL;
                        app.range(c,2) = app.maxWL;
                        app.agents(c) = "";
                    otherwise
                        disp('Channel not recognized');
                end
            end
            
            app.running = 0;
        end
        
        % Function of a switch case to call on different call backs
        % depending on which button gets selected.
        function ButtonGroupSelectionChanged(app, event)
            disp(event.Source.SelectedObject.Text);
            
            switch event.Source.SelectedObject.Text
                case 'FITC'
                    FITC_Callback(app, event);
                case 'Bodipy'
                    Bodipy_Callback(app, event);
                case 'AF488'
                    AF488_Callback(app, event);
                case 'tdTomato'
                    tdTomato_Callback(app, event);
                case 'OF550'
                    OF550_Callback(app, event);
               case 'Rhd'
                    Rhd_Callback(app, event);
                case 'AF680'
                    AF680_Callback(app, event);
                case 'AF647'
                    AF647_Callback(app, event);
                case 'IR680'
                    IR680_Callback(app, event);
               case 'OF650'
                    OF650_Callback(app, event);
                case 'PPIX'
                    PPIX_Callback(app, event);
                case 'ICG'
                    ICG_Callback(app, event);
                case 'IR800'
                    IR800_Callback(app, event);
                    
            end
            disp('button push');
        end
        
        function FITC_Callback(app, event)
            app.min470.Value = 510;
            app.max470.Value = 550;
            app.agent470 = 'FITC';
            
        end
        
        function Bodipy_Callback(app, event)
            app.min470.Value = 520;
            app.max470.Value = 550;
            app.agent470 = 'Bodipy';
        end
        
        function AF488_Callback(app, event)
            app.min470.Value = 520;
            app.max470.Value = 550;
            app.agent470 = 'AF488';
        end
        
        
        function tdTomato_Callback(app, event)
            app.min530.Value = 570;
            app.max530.Value = 650;
            app.agent530 = 'tdTomato';
        end
        
                
        function Rhd_Callback(app, event)
            app.min530.Value = 575;
            app.max530.Value = 610;
            app.agent530 = 'Rhd';
        end
        
        function OF550_Callback(app, event)
            app.min530.Value = 570;
            app.max530.Value = 640;
            app.agent530 = 'OF550';
        end
        
       function OF650_Callback(app, event)
            app.min635.Value = 660;
            app.max635.Value = 720;
            app.agent635 = 'OF650';
       end
        
        function AF680_Callback(app, event)
            app.min635.Value = 680;
            app.max635.Value = 720;
            app.agent635 = 'AF680';
        end
        
        
        function IR680_Callback(app, event)
            app.min635.Value = 680;
            app.max635.Value = 720;
            app.agent635 = 'IR680';
        end
        function AF647_Callback(app, event)
            app.min635.Value = 660;
            app.max635.Value = 680;
            app.agent635 = 'AF647';
        end
        
        function PPIX_Callback(app, event)
            app.min635.Value = 680;
            app.max635.Value = 720;
            app.agent635 = 'PPIX';
        end
        
        function IR800_Callback(app, event)
            app.agent760 = 'IR800';
        end
        
        function ICG_Callback(app, event)
            app.agent760 = 'ICG';
        end
        
        
        
    end
    
    % Component initialization
    methods (Access = private)
        
        
        % Create UIFigure and components
        function createComponents(app)
            
            % Boolean for running the app
            app.running = 1;
            
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 367 600];
            app.UIFigure.Name = 'MATLAB App';
            
            
            % Create Accept button
            app.acceptButton = uibutton(app.UIFigure);
            app.acceptButton.Position = [100 55 150 35];
            app.acceptButton.Text= 'Accept Agent Selection';
            app.acceptButton.ButtonPushedFcn = createCallbackFcn(app, @AcceptRange, true);
            
            
            % Create ChannelPanel
            app.ChannelPanel = uipanel(app.UIFigure);
            app.ChannelPanel.Title = '470 Channel';
            app.ChannelPanel.Position = [19 430 322 97];
            app.ChannelPanel.Enable = 'off';
            
            % Create minLabel
            app.minLabel = uilabel(app.ChannelPanel);
            app.minLabel.HorizontalAlignment = 'right';
            app.minLabel.Position = [97 31 34 22];
            app.minLabel.Text = 'min:';
            
            % Create min470
            app.min470 = uieditfield(app.ChannelPanel, 'numeric');
            app.min470.Position = [141 31 52 22];
            app.min470.Value = 510;
            
            % Create max470Label
            app.max470Label = uilabel(app.ChannelPanel);
            app.max470Label.HorizontalAlignment = 'right';
            app.max470Label.Position = [215 31 34 22];
            app.max470Label.Text = 'max:';
            
            % Create max470
            app.max470 = uieditfield(app.ChannelPanel, 'numeric');
            app.max470.Position = [259 31 52 22];
            app.max470.Value = 550;
            
            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.ChannelPanel);
            app.ButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup.Position = [2 0 100 74];
            
            % Create FITCButton
            app.FITCButton = uiradiobutton(app.ButtonGroup);
            app.FITCButton.Text = 'FITC';
            app.FITCButton.Position = [11 47 58 22];
            app.FITCButton.Value = true;
            
            % Create BodipyButton
            app.BodipyButton = uiradiobutton(app.ButtonGroup);
            app.BodipyButton.Text = 'Bodipy';
            app.BodipyButton.Position = [11 25 65 22];
            
            % Create AF488Button
            app.AF488Button = uiradiobutton(app.ButtonGroup);
            app.AF488Button.Text = 'AF488';
            app.AF488Button.Position = [11 3 65 22];
            
            % Create ChannelPanel_530
            app.ChannelPanel_530 = uipanel(app.UIFigure);
            app.ChannelPanel_530.Title = '530 Channel';
            app.ChannelPanel_530.Position = [20 327 322 97];
            app.ChannelPanel_530.Enable = 'off';
            
            % Create min530Label
            app.min530Label = uilabel(app.ChannelPanel_530);
            app.min530Label.HorizontalAlignment = 'right';
            app.min530Label.Position = [97 31 34 22];
            app.min530Label.Text = 'min:';
            
            % Create min530
            app.min530 = uieditfield(app.ChannelPanel_530, 'numeric');
            app.min530.Position = [141 31 52 22];
            app.min530.Value = 570;
            
            % Create max530Label
            app.max530Label = uilabel(app.ChannelPanel_530);
            app.max530Label.HorizontalAlignment = 'right';
            app.max530Label.Position = [215 31 34 22];
            app.max530Label.Text = 'max:';
            
            % Create max530
            app.max530 = uieditfield(app.ChannelPanel_530, 'numeric');
            app.max530.Position = [259 31 52 22];
            app.max530.Value = 650;
            
            % Create ButtonGroup_2
            app.ButtonGroup_2 = uibuttongroup(app.ChannelPanel_530);
            app.ButtonGroup_2.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup_2.Position = [3 1 100 74];
            
            % Create tdTomatoButton
            app.tdTomatoButton = uiradiobutton(app.ButtonGroup_2);
            app.tdTomatoButton.Text = 'tdTomato';
            app.tdTomatoButton.Position = [11 47 72 22];
            app.tdTomatoButton.Value = true;
            
            app.OF550Button = uiradiobutton(app.ButtonGroup_2);
            app.OF550Button.Text = 'OF550';
            app.OF550Button.Position = [11 25 72 22];
            
            app.RhdButton = uiradiobutton(app.ButtonGroup_2);
            app.RhdButton.Text = 'Rhd';
            app.RhdButton.Position = [11 3 72 22];
            
            % Create ChannelPanel_635
            app.ChannelPanel_635 = uipanel(app.UIFigure);
            app.ChannelPanel_635.Title = '635 Channel';
            app.ChannelPanel_635.Position = [18 200 322 125];
            app.ChannelPanel_635.Enable = 'off';
            
            % Create min635Label
            app.min635Label = uilabel(app.ChannelPanel_635);
            app.min635Label.HorizontalAlignment = 'right';
            app.min635Label.Position = [97 31 34 22];
            app.min635Label.Text = 'min:';
            
            % Create min635
            app.min635 = uieditfield(app.ChannelPanel_635, 'numeric');
            app.min635.ValueChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.min635.Position = [141 31 52 22];
            app.min635.Value = 680;
            
            % Create max635Label
            app.max635Label = uilabel(app.ChannelPanel_635);
            app.max635Label.HorizontalAlignment = 'right';
            app.max635Label.Position = [215 31 34 22];
            app.max635Label.Text = 'max:';
            
            % Create max635
            app.max635 = uieditfield(app.ChannelPanel_635, 'numeric');
            app.max635.ValueChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.max635.Position = [259 31 52 22];
            app.max635.Value = 720;
            
            % Create ButtonGroup_3
            app.ButtonGroup_3 = uibuttongroup(app.ChannelPanel_635);
            app.ButtonGroup_3.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup_3.Position = [3 1 100 100];
            
            % Create IR680Button
            app.IR680Button = uiradiobutton(app.ButtonGroup_3);
            app.IR680Button.Text = 'IR680';
            app.IR680Button.Position = [11 75 54 22];
            app.IR680Button.Value = true;
            
            % Create AF680Button
            app.AF680Button = uiradiobutton(app.ButtonGroup_3);
            app.AF680Button.Text = 'AF680';
            app.AF680Button.Position = [11 50 65 22];
            
            % Create AF647Button
            app.AF647Button = uiradiobutton(app.ButtonGroup_3);
            app.AF647Button.Text = 'AF647';
            app.AF647Button.Position = [11 25 65 22];
            
            app.OF650Button = uiradiobutton(app.ButtonGroup_3);
            app.OF650Button.Text = 'OF650';
            app.OF650Button.Position = [11 1 72 22];
            
            
            %             % Create PPIXButton
            %             app.PPIXButton = uiradiobutton(app.ButtonGroup_3);
            %             app.PPIXButton.Text = 'PPIX';
            %             app.PPIXButton.Position = [11 3 65 22];
            
            
            % Create ChannelPanel_760
            app.ChannelPanel_760 = uipanel(app.UIFigure);
            app.ChannelPanel_760.Title = '760 Channel';
            app.ChannelPanel_760.Position = [18 100 322 97];
            app.ChannelPanel_760.Enable = 'off';
            
            % Create min760Label
            app.min760Label = uilabel(app.ChannelPanel_760);
            app.min760Label.HorizontalAlignment = 'right';
            app.min760Label.Position = [97 31 34 22];
            app.min760Label.Text = 'min:';
            
            % Create min760
            app.min760 = uieditfield(app.ChannelPanel_760, 'numeric');
            app.min760.Position = [141 31 52 22];
            app.min760.Value = 780;
            app.min760.Enable = 'off';
            
            % Create max530Label
            app.max760Label = uilabel(app.ChannelPanel_760);
            app.max760Label.HorizontalAlignment = 'right';
            app.max760Label.Position = [215 31 34 22];
            app.max760Label.Text = 'max:';
            
            % Create max530
            app.max760 = uieditfield(app.ChannelPanel_760, 'numeric');
            app.max760.Position = [259 31 52 22];
            app.max760.Value = 850;
            app.max760.Enable = 'off';
            
            % Create ButtonGroup_4
            app.ButtonGroup_4 = uibuttongroup(app.ChannelPanel_760);
            app.ButtonGroup_4.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup_4.Position = [3 1 100 74];
            
            % Create IR800Button
            app.IR800Button = uiradiobutton(app.ButtonGroup_4);
            app.IR800Button.Text = 'IR800';
            app.IR800Button.Position = [11 47 54 22];
            app.IR800Button.Value = true;
            
            % Create ICGButton
            app.ICGButton = uiradiobutton(app.ButtonGroup_4);
            app.ICGButton.Text = 'ICG';
            app.ICGButton.Position = [11 25 65 22];
            
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
            
            % Iterate through the channels and enable the channels that
            % should be turned on.
            for c = 1:length(app.channels)
                switch(app.channels(c))
                    case "470"
                        app.ChannelPanel.Enable = 'on';
                        app.agent470 = "FITC";
                    case "530"
                        app.ChannelPanel_530.Enable = 'on';
                        app.agent530 = "tdTomato";
                    case "635"
                        app.ChannelPanel_635.Enable = 'on';
                        app.agent635 = "IR680";
                    case "760"
                        app.ChannelPanel_760.Enable = 'on';
                        app.agent760 = "IR800";
                    case "WL"
                        app.minWL = 450;
                        app.maxWL = 510;
                    otherwise
                        disp('Channel not recognized');
                end
            end
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = setRange_app(channels)
            
            % give the app the correct channels to display
            app.channels = channels;
            
            % Create UIFigure and components
            createComponents(app)
            
            % Register the app with App Designer
            registerApp(app, app.UIFigure)
            
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