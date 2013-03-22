function supfig(figfilename1,figfilename2)

% SUPFIG(figfilename1,figfilename2)
% opens fig files and places the contents into a new figure.

% Import fig1 info
ImportFig1 = hgload(figfilename1); % Open fig1 file
AxesFig1 = get(ImportFig1,'CurrentAxes'); % Get the axes
LinesFig1 = get(AxesFig1,'child'); % Get the lines (ie the data)

% Import fig2 info
ImportFig2 = hgload(figfilename2); % Open fig2 file
AxesFig2 = get(ImportFig2,'CurrentAxes'); % Get the axes
LinesFig2 = get(AxesFig2,'child'); % Get the lines (ie the data)

%for ii = 1:length(LinesFig2)
%set(LinesFig2(ii),'Color','y'); % Bonus
%end

% Concatenate data
LinesFigs = [ LinesFig2 ; LinesFig1 ];

% Create new figure
figure;
NewFig = gcf;
copyobj(AxesFig1,NewFig);
copyobj(AxesFig2,NewFig);
ChildNewFig = get(NewFig,'child');
copyobj(LinesFigs,ChildNewFig);

% Delete imported figures
delete(ImportFig1);
delete(ImportFig2);
