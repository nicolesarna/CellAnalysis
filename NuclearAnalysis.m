clc;clear;
% Display a message box to guide users on selecting image files. Instructs on multi-selection and directory constraint.
waitfor(msgbox({'Please select images to process.' 'Use shift + left mouse to select more than one file.',...
    'Files may only be selected from one directory at a time.'},'Select Image Help'));

% Open a file selection dialog for the user to choose image files. Supports multiple selections and various image formats.
[fname1,pname1]=uigetfile({'*.*';'*.jpg';'*.png';'*.tiff';'*.czi'},'Select Image(s)','MultiSelect', 'on');

% Convert the selected file names into a cell array if not already.
fname1 = cellstr(fname1);

% Initialize variables to store processing results. 
% Include variables for image number, cell count, areas, intensities, circularity, eccentricity, diameter, and solidity.
imagenum = [];
cellnum = [];
nucGCTCF = [];
allarea = [];
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
bwR1 = imbinarize(R1filt,0.15);
bwG1 = imbinarize(G1filt,0.1);
bwB1 = imbinarize(B1filt,0.1);

% Combine filtered channels for visualization. This is not used for analysis, just for visual inspection.
visual = imadd(R1,B1);
allchannels = imadd(visual,G1);
allchannelsfilt = medfilt2(allchannels);
bwall = imbinarize(allchannelsfilt,0.35);

% Process blue channel (nuclear stain) to identify nuclei, filling holes and removing small objects
Bchannel = bwB1;
se1=strel("diamond",5);
bwB2 = imopen(Bchannel,se1);
bwB3 = imfill(bwB2,'holes');
bwB4 = bwareaopen(bwB3, 200);
bwB5 = imclearborder(bwB4);
Boverlay= labeloverlay(B1,bwB5);

% Process the combined channel image for background segmentation.
allbkg = bwall;
se3 = strel("diamond",2);
allbkg2 = imerode(allbkg,se3);
allbkg3 = bwareaopen(allbkg2, 20);
bkgse = strel('diamond',40);
allbkg4 = imdilate(allbkg3,bkgse);
allbkg5 = imfill(allbkg4, 'holes');
bkgfin = allbkg5;
bkgoverlay = labeloverlay(allchannelsfilt, bkgfin);

% Apply Watershed algorithm to segment closely touching nuclei.
D = -bwdist(~bwB5);
ld = watershed(D);
bw2 = bwB5;
bw2(ld == 0) = 0;
mask = imextendedmin(D,2);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
nucbw3 = bwB5;
nucbw3(Ld2 == 0) = 0;
nucbw4 = imclearborder(nucbw3); % Remove objects touching the border.
finnuc = bwareafilt(nucbw4, [2000 100000000]); % Filter objects by size.

% Quantify properties of segmented nuclei.
nuc = regionprops(finnuc, 'all');

% Draw the perimeter around objects for visualization.
nucbw_perim = bwperim(finnuc); 
nucbw_perim=imdilate(nucbw_perim,strel('disk',4)); 
nucoverlay = imoverlay(visual, nucbw_perim, 'magenta'); %Overlay the perimeter of objects on image

%Background Quantification
bkg = ~bkgfin; % Inverts binary objects mask to quantify background instead of foreground
BkgData = regionprops(bkg, 'all'); % Collect properties of the background
  
    % Convert pixel locations from the background mask to linear indices
   BkgInds = sub2ind(size(greenchannel), BkgData(1).PixelList(:,2),BkgData(1).PixelList(:,1)); 
   AvgGreenBkg = mean(greenchannel(BkgInds)); % Mean intensity of background on green channel

% Foreground Quantification
% Initialize vectors for storing nuclear properties.
nucgreenint=zeros(1,length(nuc)); % For average intensity of each nucleus in the green channel.
nucarea=zeros(1,length(nuc)); % For area of each nucleus.

% Loop through each nucleus to calculate its properties.
for n=1:1:length(nuc)

   % Convert nuclear pixel coordinates to linear indices.
   nucind = sub2ind(size(greenchannel), nuc(n).PixelList(:,2),nuc(n).PixelList(:,1)); % Indexes pixel list on original green channel
   nucgreenint(n) = mean(greenchannel(nucind)); % Calculate and store the average green intensity.
   nucarea(n) = nuc(n).Area; % Store the area.

end

% Accumulate data for all nuclei and images processed.
imagenum = [imagenum, p+zeros(1,n)]; % Track which image each nucleus belongs to.
cellnum = [cellnum, 1:1:n]; % Assign a unique identifier to each nucleus.
allarea = [allarea nucarea]; % Collect areas of all nuclei.

% Calculate corrected total cellular fluorescence (CTCF) for each nucleus.
GreenIntDen = nucgreenint.*nucarea; % Multiply mean intensity by area for green channel.
Greencorv = nucarea.*AvgGreenBkg; % Background correction for green channel.
Gintall = GreenIntDen - Greencorv; % Corrected total green fluorescence.
nucGCTCF = [nucGCTCF Gintall]; % Accumulate for all nuclei.


    % Save processed images if the user opted to, with perimeter overlays and object number annotations
    if save == 1
        fg= figure('visible','off'); % Invisible figure for displaying and saving the overlay image with annotations.
        imshow(nucoverlay) 

            % Annotate each detected object with its number.
            for j = 1:1:length(nuc) % Loop through all counted cells
                itemCentroid = nuc(j).Centroid; % Retrieve object centroid
                text(itemCentroid(1), itemCentroid(2), num2str(j), 'Color', [1 0.2 0.2], 'fontsize',20); % Write object number over its centroid.
            end
                title("Total Count: " +length(nuc)); % Display total count of detected objects.
        
                % Save the annotated image to file and close figure.
                saveas(gcf,fullfile(pname1, [ID '_nuc.png']), 'png');
                close(fg);
    end
    
    close(figure(1)); % Close the figure window after processing each file.

end
 
% Data Aggregation and Export
% Convert collected data into columns for table creation.
colimagenum = imagenum';
colcellnum = cellnum';
colallarea = allarea';
colnucgreenintall = nucGCTCF';

% Combine collected data into a matrix and create a table with descriptive headers.
    datamat = horzcat(colimagenum, colcellnum, colallarea, colnucgreenintall);
    T = array2table(datamat,'VariableNames',{'File', 'Cell', 'NuclearArea', 'NuclearGreenCTCF'});

% Save the table to an Excel file in the selected directory.
    fullFileName=  fullfile(pname1,'NucleusAnalysis.xls');
    writetable(T, fullFileName);