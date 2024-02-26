clc; clear

% Display a message box to guide users on selecting image files. Instructs on multi-selection and directory constraint. 
waitfor(msgbox({'Please select images to process.' 'Use shift + left mouse to select more than one file.',...
    'Files may only be selected from one directory at a time.'},'Select Image Help'));

% Open a file selection dialog for the user to choose image files. Supports multiple selections and various image formats
[fname1,pname1]=uigetfile({'*.*';'*.jpg';'*.png';'*.tiff';'*.czi'},'Select Image(s)','MultiSelect', 'on');

%Convert the selected file names into a cell array if not already.
fname1 = cellstr(fname1);

%Initialize variables to store processing results. 
%Include variables for image number, cell count, areas, intensities, circularity, eccentricity, diameter, and solidity.
imagenum = [];
cellnum = [];
allarea = [];
truegreenintall = [];
trueredintall = [];
allCirc = [];
allEccen = [];
allDia = [];
allSolid = [];
save = 0;

% Ask the user whether they want to save images with circled cells through a dialog box.
        choice = questdlg('Would you like to save the images with circled cells','Save Images?','Yes','No','No');
        % Handle the user's choice for saving images.
        switch choice
            case 'Yes'
                save=1;
            case 'No'
                save=0;
        end          

% Main Program
% Iterates through each selected file for processing.
for p = 1:length(fname1)

    % Concatenate path and file name to obtain full path to the image.
    filename = fullfile ( pname1, fname1{p} );
    % Extract the file name (without extension) to use as an identifier.
    [~,name,~] = fileparts(filename);  
    ID = name;

% Read the current image file using Bio-Formats.
r = bfGetReader(filename);

% Separate and define RGB color channels. Assumes the first three channels are Red, Green, and Blue respectively.
redchannel = bfGetPlane(r, 1);
greenchannel = bfGetPlane(r, 2);
bluechannel = bfGetPlane(r, 3);

% Enhance each channel for better visualization and apply median filter to reduce noise.
R1= imadjust(redchannel);
G1= imadjust(greenchannel); 
B1 = imadjust(bluechannel);
R1filt = medfilt2(R1);
G1filt = medfilt2(G1);
B1filt = medfilt2(B1);

% Binarize filtered channels using predefined thresholds to separate foreground from background.
bwR1 = imbinarize(R1filt,0.1);
bwG1 = imbinarize(G1filt,0.4);
bwB1 = imbinarize(B1filt,0.1);

% Combine filtered channels for visualization. This is not used for analysis, just for visual inspection.
visual = imadd(R1filt, B1filt);

% Process blue channel (nuclear stain) to identify objects, filling holes and removing small objects
Bchannel = bwB1;
se1=strel("diamond",5);
bwB2 = imopen(Bchannel,se1);
bwB3 = imfill(bwB2,'holes');
bwB5 = bwareaopen(bwB3, 200);
Boverlay= labeloverlay(visual,bwB5);

% Similar processing for the red channel (cell body), with adjustments for object size and shape.
Rchannel = bwR1;
se2=strel("diamond",4);
bwR2 = imdilate(Rchannel,se2);
bwR3 = imfill(bwR2,'holes');
bwR4 = bwareaopen(bwR3, 400);
se3 = strel("disk",6);
bwR40 = imerode(bwR4,se3);
bwR5 = bwareafilt(bwR40, [10000 100000000]);
Roverlay= labeloverlay(visual,bwR5);

% Background processing to differentiate foreground objects from the background.
Rbkg = bwR5;
bkgse = strel('diamond',40);
Rbkg1 = imdilate(Rbkg,bkgse);
Rbkg2 = imfill(Rbkg1, 'holes');
bkg = ~Rbkg2; %inverts mask to quantify background instead of foreground
bkgoverlay = labeloverlay(visual, bkg);

%Apply Watershed algorithm to segment closely touching objects.
D = -bwdist(~bwB5);
ld = watershed(D);

bw2 = bwB5;
bw2(ld == 0) = 0;
mask = imextendedmin(D,2);

D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bwR5;
bw3(Ld2 == 0) = 0;
bw4 = imclearborder(bw3); % Remove objects touching the border.
bw5 = bwareafilt(bw4, [2000 100000000]); % Filter objects by size.

% Measure properties of labeled regions.
objects = regionprops(bw5, "Area","Centroid","Circularity","Eccentricity","EquivDiameter","PixelList","Solidity");

% Draw the perimeter around objects for visualization.
bw_perim = bwperim(bw5); 
bw_perim = imdilate(bw_perim,strel('disk',4)); 
overlay1 = imoverlay(visual, bw_perim, 'magenta'); %Overlay the perimeter of objects on image

% Background Quantification for both Green and Red channels.
BkgData = regionprops(bkg, 'all'); 

   % Convert pixel locations from the background mask to linear indices for green and red channels
   BkgInds = sub2ind(size(greenchannel), BkgData(1).PixelList(:,2),BkgData(1).PixelList(:,1)); %indexes pixel list in original image
   AvgGreenBkg = mean(greenchannel(BkgInds)); % Mean intensity of background on green channel
   AvgRedBkg = mean(redchannel(BkgInds)); % Mean intensity of background on red channel

% Foreground Quantification for both Green and Red channels.
% Initialize vectors for storing quantified data of detected objects.
v=zeros(1,length(objects)); % Vector for mean green intensity of objects.
rev=zeros(1,length(objects)); % Vector for mean red intensity of objects.
varea=zeros(1,length(objects)); % Vector for area of objects.
cellcirc = zeros(1,length(objects)); % Vector for circularity of objects.
celleccen = zeros(1,length(objects)); % Vector for eccentricity of objects.
cellDia = zeros(1,length(objects)); % Vector for equivalent diameter of objects.
cellSolid = zeros(1,length(objects)); % Vector for solidity of objects.

% Loop through each detected object to quantify its properties.
for n=1:1:length(objects)

   % Convert object pixel locations to linear indices for intensity measurement.
   ind = sub2ind(size(greenchannel), objects(n).PixelList(:,2),objects(n).PixelList(:,1)); 
   v(n) = mean(greenchannel(ind)); % Store mean green intensity for the object.
   varea(n) = objects(n).Area; % Store area of the object.

   % Repeat for red channel.
   ind = sub2ind(size(redchannel), objects(n).PixelList(:,2),objects(n).PixelList(:,1)); 
   rev(n) = mean(redchannel(ind)); % Store mean red intensity for the object.
  
   % Store geometric properties of the object.
   cellcirc(n) = objects(n).Circularity;
   celleccen(n) = objects(n).Eccentricity;
   cellDia(n) = objects(n).EquivDiameter;
   cellSolid(n) = objects(n).Solidity;

end

% Aggregate quantified data for all objects.
imagenum = [imagenum, p+zeros(1,n)]; % Track image numbers for data association.
cellnum = [cellnum, 1:1:n]; % Sequential numbers of detected objects.
allarea = [allarea varea]; % Accumulate areas of all objects.
allCirc = [allCirc, cellcirc]; % Accumulate circularity measurements.
allEccen = [allEccen, celleccen]; % Accumulate eccentricity measurements.
allDia = [allDia, cellDia]; % Accumulate equivalent diameters.
allSolid = [allSolid, cellSolid]; % Accumulate solidity measurements.

% Calculate corrected total cellular fluorescence (CTCF) for green channels.
GreenIntDen = v.*varea; % Multiply mean intensity by area for green channel.
Greencorv = varea.*AvgGreenBkg; % Background correction for green channel.
truegreenint = GreenIntDen - Greencorv; % Corrected total green fluorescence.
truegreenintall = [truegreenintall truegreenint]; % Accumulate for all objects.

% Calculate corrected total cellular fluorescence (CTCF) for red channel
RedIntDen = rev.*varea; % Multiply mean intensity by area for red channel.
Redcorv = varea.*AvgRedBkg; % Background correction for red channel.
trueredint = RedIntDen - Redcorv; % Corrected total red fluorescence.
trueredintall = [trueredintall trueredint]; % Accumulate for all objects.

 
% Save processed images if the user opted to, with perimeter overlays and object number annotations
if save == 1
     fg= figure('visible','off'); % Invisible figure for displaying and saving the overlay image with annotations.
     imshow(overlay1) 

        % Annotate each detected object with its number.
        for k = 1:1:length(objects)
            itemCentroid = objects(k).Centroid;	% Retrieve object centroid
            text(itemCentroid(1), itemCentroid(2), num2str(k), 'Color', [1 0.2 0.2], 'fontsize',25); % Write object number over its centroid.
        end
        
        title("Total Count: " +length(objects)); % Display total count of detected objects.

        % Save the annotated image to file and close figure.
         saveas(gcf,fullfile(pname1, [ID '_Analysis.png']), 'png');
         close(fg);
end

close(figure(1)); % Close the figure window after processing each file.
    
end

% Data Aggregation and Export
% Convert collected data into columns for table creation.
colimagenum = imagenum';
colcellnum = cellnum';
colallarea = allarea';
colCirc = allCirc';
colEccen = allEccen';
colDia = allDia';
colSolid = allSolid';
coltruegreenintall = truegreenintall';
coltrueredintall = trueredintall';

%Combine collected data into a matrix and create a table with descriptive headers.
    datamat = horzcat(colimagenum, colcellnum, colallarea, colCirc, colEccen, colDia, colSolid, coltrueredintall, coltruegreenintall);

    T = array2table(datamat,'VariableNames',{'File', 'Cell', 'Area','Circularity','Eccentricity','EquivDiameter','Solidity', 'Red CTCF', 'GreenCTCF'});

% Save the table to an Excel file in the selected directory.
    fullFileName=  fullfile(pname1,'CellAnalysis.xls');
    writetable(T, fullFileName);