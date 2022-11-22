function[channelRange] =  setchannelRange(channels)

fig = uifigure('Position', [ 100 100 330 270]);
lblTitle = uilabel(fig,'Position',[75 (230)  200 22],'Text','Select Channel Bandwidths');

fig.Resize = 'Off';
global channelRange;

% Set up the initial field condition
for c = 1:length(channels)
    switch channels(c)
        case '470'
            channelRange(c,1) = 510;
            channelRange(c,2) = 550;
        case '532'
            channelRange(c,1) = 570;
            channelRange(c,2) = 650;
        case '530'
            channelRange(c,1) = 570;
            channelRange(c,2) = 650;
        case '635'
            channelRange(c,1) = 680;
            channelRange(c,2) = 720;
        case '760'
            channelRange(c,1) = 0;
            channelRange(c,2) = 0;
        case 'WL'
            channelRange(c,1) = 450;
            channelRange(c,2) = 510;
        otherwise
            error(strcat('Channel Unrecognized: ', channels(c)));
    end
    
end

% set the figure's property for the channel range to default values
fig.Parent.UserData.channelRange = channelRange;
%% SET UP INITIAL VALUES
for c = 1:length(channels)
    pnl = uipanel(fig, 'Position', [5 (195 - ((c-1) * 30)) 300 35]);
    
    switch channels(c)
        case '470'
            
            min470 = uieditfield(fig,'numeric',...
                'Position',[130 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,1),'ValueChangedFcn',@(min470,event) textChanged(min470,channelRange,c,1));
            max470 = uieditfield(fig,'numeric',...
                'Position',[225 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,2),'ValueChangedFcn',@(max470,event) textChanged(max470,channelRange,c,2));
            
        case '532'
            min532 = uieditfield(fig,'numeric',...
                'Position',[130 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,1),'ValueChangedFcn',@(min532,event) textChanged(min532,channelRange,c,1));
            max532 = uieditfield(fig,'numeric',...
                'Position',[225 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,2),'ValueChangedFcn',@(max532,event) textChanged(max532,channelRange,c,2));
            
        case '530'
            min530 = uieditfield(fig,'numeric',...
                'Position',[130 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,1),'ValueChangedFcn',@(min530,event) textChanged(min530,channelRange,c,1));
            max530 = uieditfield(fig,'numeric',...
                'Position',[225 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,2),'ValueChangedFcn',@(max530,event) textChanged(max530,channelRange,c,2));

        case '635'
            
            min635 = uieditfield(fig,'numeric',...
                'Position',[130 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,1),'ValueChangedFcn',@(min635,event) textChanged(min635,channelRange,c,1));
            max635 = uieditfield(fig,'numeric',...
                'Position',[225 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,2),'ValueChangedFcn',@(max635,event) textChanged(max635,channelRange,c,2));
            
        case '760'
            
            min760 = uieditfield(fig,'numeric',...
                'Position',[130 (200 - ((c-1) * 30)) 40 22],'Value',780);
            max760 = uieditfield(fig,'numeric',...
                'Position',[225 (200 - ((c-1) * 30)) 40 22],'Value',900);
            
            min760.Enable = 'Off';
            max760.Enable = 'Off';
            
        case 'WL'
            minWL = uieditfield(fig,'numeric',...
                'Position',[130 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,1),'ValueChangedFcn',@(minWL,event) textChanged(minWL,channelRange,c,1));
            maxWL = uieditfield(fig,'numeric',...
                'Position',[225 (200 - ((c-1) * 30)) 40 22],'Value',channelRange(c,2),'ValueChangedFcn',@(maxWL,event) textChanged(maxWL,channelRange,c,2));
            
        otherwise
            error(strcat('Channel Unrecognized: ', channels(c)));
            
    end
    
    lblmin = uilabel(fig,'Position',[105 (200 - ((c-1) * 30))  40 22],'Text','Min:');
    lbmax= uilabel(fig,'Position',[200 (200 - ((c-1) * 30))  40 22],'Text','Max:');
    lblChannel = uilabel(fig,'Position',[10 (200 - ((c-1) * 30)) 100 22],'Text',strcat('Channel ',{' '}, channels(c),':'));
    
end

acceptButton = uibutton(fig, 'push', ....
    'Position', [50 (150 - ((c-1) * 30)) 80 35],...
    'Text', 'Accept','ButtonPushedFcn', @(acceptButton,event) acceptchannelRanges(acceptButton, channelRange));

cancelButton = uibutton(fig, 'push', ....
    'Position', [150 (150 - ((c-1) * 30)) 80 35],...
    'Text', 'Cancel','ButtonPushedFcn', @(cancelButton,event) acceptchannelRanges(cancelButton));

uiwait();  % # <--------------------- Output will be available when the GUI closes (or ui is resumed)
channelRange = fig.Parent.UserData.channelRange;
close(fig);

% Code the callback function.
    function[channelRange] = textChanged(hdl,channelRange,c,ind)
        % gather the channel range from the fig parent properties
        channelRange = fig.Parent.UserData.channelRange;
        channelRange(c,ind) = hdl.Value;
        uiresume(); % resume gui to change the channel range
        fig.Parent.UserData.channelRange = channelRange;
        disp('channelRange successfully changed');
        uiwait(); % pause the gui to allow the user to continue to select data
    end


    function[channelRange] =  acceptchannelRanges(btn, channelRange)
        if (strcmp(btn.Text,'Accept'))
            pause(0.2);
            channelRange = fig.Parent.UserData.channelRange
            uiresume();  % resume gui to finalize the values
            disp('Finalized channelRange');
        else
            channelRange = [100 100];
        end
    end

    function [channelRange] = returnVals()
        channelRange = fig.Parent.UserData.channelRange;
    end

end


