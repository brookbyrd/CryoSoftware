function reportMacroError(config)
% Function created to check on running macrotome status by looking at how
% many slices there are in the newest folder every 60 sec.
% Created by Brook Byrd 8-22-21

% Emails list when ever there is a disturbance.

addpath(pwd);

stat=true;
time0 = tic;
timeLimit = 60*15*1; % 15 min time limit

% Initially check how many slices have been acquired today
numSlices_old = checkSlices(config);
numSlices_current = numSlices_old;

% loop to check # of slices acquired
while(stat==true)
    
    
    % Still on the same slice
    if(numSlices_current == numSlices_old)
        numSlices_old = numSlices_current;
        numSlices_current = checkSlices(config); % Continously check the number of slices
        numSlices_current
        
        % if we reach our time limit (15 min), email the warning
        if(toc(time0)>timeLimit)
            
            if(numSlices_current >= config.numberOfSlices(1)) % EXPERIMENT COMPLETE!
                 disp(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,4), {' '},'of todays experiment'));
                sendolmail('f00355n@dartmouth.edu',char(strcat('MACROTOME FINISHED:',{' '},config.studyName,',', {' '}, string(datetime('now')))),char(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has completed the ', {' '}, config.studyName, {' '}, ' experiment on slice #',num2str(numSlices_current,4), {' '})));
             
            else
                disp(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,4), {' '},'of todays experiment'));
                sendolmail('f00355n@dartmouth.edu',char(strcat('MACROTOME DISTURBANCE:',{' '},config.studyName,',', {' '}, string(datetime('now')))),char(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,4), {' '},'of todays experiment')));
                sendolmail('f005cgp@dartmouth.edu',char(strcat('MACROTOME DISTURBANCE:',{' '},config.studyName,',' ,{' '}, string(datetime('now')))),char(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,2), {' '},'of todays experiment')));
                %             sendolmail('d12184@dartmouth.edu',char(strcat('MACROTOME DISTURBANCE:',{' '},config.studyName,',' ,{' '}, string(datetime('now')))),char(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,2), {' '},'of todays experiment')));
                %             sendolmail('d12184@dartmouth.edu',char(strcat('MACROTOME DISTURBANCE:',{' '},config.studyName,',' ,{' '}, string(datetime('now')))),char(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,2), {' '},'of todays experiment')));
                %             sendolmail('f002p4s@dartmouth.edu',char(strcat('MACROTOME DISTURBANCE:',{' '},config.studyName,',' ,{' '}, string(datetime('now')))),char(strcat('Beginning at',{' '}, string(datetime('now')),' , the Davis macrotome has been inactive for 15 min. while slicing ', {' '}, config.studyName, {' '}, 'on slice #',num2str(numSlices_current,2), {' '},'of todays experiment')));
            end
            
            % break the for loop
            stat = false;
        else
            pause(30); % only run this once a minute to save on processing load
        end
        
    else % we've jumped to a new slice! record the data and email milestones
        % record the amount of time it took
        timePerSlice = toc(time0);
        disp(strcat('Processed Slice #',{' '},num2str(numSlices_current) ,{' '}, 'in', {' '}, num2str(timePerSlice/60,4), {' '}, 'min'));
        
        if(mod( numSlices_current, 50) == 0)
            sendolmail('f00355n@dartmouth.edu',char(strcat('MACROTOME UPDATES:', {' '},config.studyName,',',{' '}, string(datetime('now')))), char(strcat('Today we have processed ',num2str(numSlices_current) ,{' '},...
                'The most recent slice took', {' '}, num2str(timePerSlice/60,4), {' '}, 'minutes to complete.')));
        end
        
        % Reset the number of slices to be the new #
        numSlices_old = numSlices_current;
        time0 = tic% Restart the clock
    end
end

