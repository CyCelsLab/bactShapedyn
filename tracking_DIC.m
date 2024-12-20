%% Modified by Tanvi Kale, October 2023

outFolder = '/Users/tanvikale/Desktop/CA/Bac_cell_size/mg1655_a22_27092023/'; 
filename  =  '/Users/tanvikale/Desktop/CA/Bac_cell_size/mg1655_a22_27092023/2002_85frames-singlecell_adjusted.tif';
all_data.Area = ones(height(all_data), 1); 
[tracks, obno, outmat] = DICOT_tracking(outFolder, all_data, 10, 1, 'um', 's'); 
finfo = imfinfo(filename); 
makemovie(outFolder, 1:length(finfo), filename, outmat)
%% Show trajectories 

%%
function [tracks,objno,outmat] = DICOT_tracking(outfolder, savestats, micron_search_radius, scal_fact, interval, distUnit, timeUnit)
%DICOT_TRACKING tracking code for DICOT, modified version of
% trackfil_dicot
%{
By default function removes all tracks that are less than 3 units in
tracking times
%}
Centroid_coordinates = cat(1, savestats.Centroid);
pixel_search_radius= micron_search_radius;  %/scal_fact;
x = Centroid_coordinates(:,2)';
y = Centroid_coordinates(:,1)';
frame = cat(1,savestats.Frame);
len= cat(1, savestats.Length); 
width = cat(1, savestats.Width); 
beadlabel=zeros(size(x)); % vector of bead labels.
i=min(frame); spanA=find(frame==i); % initialize w/ first frame.
beadlabel(1:length(spanA))=1:length(spanA); % refers to absolute indexing of x,y,frame,etc.
lastlabel=length(spanA); % start off with unique bead labels for all the beads in the first frame

w =waitbar(0, 'Tracking objects..');
p=0;
fidx= min(frame):max(frame)-1;
tic; 
for i= fidx
    p=p+1;
    waitbar(p/numel(fidx));
    spanA= find(frame==i); % ARC 'find'-here, gives positions 
    spanB= find(frame==i+1);
    dx = ones(length(spanA),1)*x(spanB) - x(spanA)'*ones(1,length(spanB));
    dy = ones(length(spanA),1)*y(spanB) - y(spanA)'*ones(1,length(spanB));
    
    dr2 = sqrt(dx.^2 + dy.^2); 
     [from, to, orphan] = beadsorterMod(dr2, pixel_search_radius); 
    from=spanA(from);
    to=spanB(to);
    orphan=spanB(orphan);
    beadlabel(to)= beadlabel(from); 
    % modified ARC 28/5/2016
    if ~isempty(orphan)
        beadlabel(orphan)=lastlabel+(1:length(orphan));
        lastlabel=lastlabel+length(orphan); 
    end
end
t = toc; 
disp(t)
emptybead.x=0; emptybead.y=0; emptybead.area=0; emptybead.frame=0;

%% Initialize for purposes of speed and memory management.
ALLtracks=cell(lastlabel,1); % modified ARC 28/5/2016
re=0;

for i=1:lastlabel
    p=p+1;
    waitbar(p/(numel(fidx)+lastlabel));% Y [mod]% reassemble beadlabel into a structured array 'tracks' containing all info
    
    beadi=find(beadlabel==i);
    if  numel(x(beadi))>=3 % ARC minimum no. of data points in a trajectory
    re=re+1;%ARC renumbering
    tracks(re).x=x(beadi);
    tracks(re).y=y(beadi);
    tracks(re).len=len(beadi)*scal_fact;
    tracks(re).width=width(beadi)*scal_fact;
    tracks(re).frame=frame(beadi);
    % modified ARC 28/5/2016
    ALLtracks{re}= [(re*ones(size(tracks(re).x)))',(tracks(re).frame), (tracks(re).x)',...
        (tracks(re).y)', (tracks(re).len), (tracks(re).width),((tracks(re).frame))*interval];
    %ALLtracks{re}(:,6)=ALLtracks{re}(:,6)-ALLtracks{re}(1,7);
    % object number, frame number, centroid(x,y)(um),...
    % length of skeleton(um), time
    else
        continue
    end
end
%ARC 25/2/2018:
%% write data to file
outmat=cat(1,ALLtracks{:});
objno=re;
if exist(outfolder, 'dir')==0
   mkdir(outfolder)
end

if exist([outfolder, '/trajectories.txt'], 'file') 
    delete([outfolder,'/trajectories.txt']);
end

%fid =fopen([outfolder,'/trajectories.txt'], 'w');
%fprintf(fid, ['ObjID	Frame	X	Y	Length (',distUnit,')	Time (', timeUnit, ')\r\n']);
%writematrix([outfolder, '/trajectories.txt'], outmat,'-append',...
%    'delimiter', '\t','newline', 'pc', 'precision', '%.2f');
%fclose(fid);

delete(w);

end

%%-------------------------------------------------------------------------
function [from,to, orphan]= beadsorterMod(dr2, pixel_search_radius)
% All bead tracking is done here.  Everything else is bookkeeping.  NOT ROBUST.  Look here first for problems!
% find indices of minimum value in each row 
% Get all hits within the search radius 
[i,j]  = min(dr2,[],2); % j is the row index
orphan=find(i > pixel_search_radius); 
ColArray = 1:size(dr2,1);
ColArray(orphan) = [];
j(orphan) = []; 
from = ColArray'; 
to = j; 
orphan=setdiff(1:size(dr2,2),to);

end
%%
%% 24/2/2018 Anushree, iiser pune
%% Modified 7/3/2018
%% Overlay tracks on images in real time i.e. make a tracking movie
%% Yash MOD
function makemovie(infolder, Findex,imagename, outmat)

%% INPUT:
% Findex = frame nos. start:end
% imagename = tif time-series
% ALLtracks = tracked output matrix from 'link_cells4.m' 
% (1:obj ID, 2:x (pixels), 3:y (pixels), 4:time (specific units), 5:frame, 6:length (specific units))
%% OUTPUT:
% trackingmovie.tif = output tif movie
%% CODE:

outfile=[infolder,'/trackingmovie.tif'];
if exist(outfile, 'file')
    delete(outfile);
end
ALLtracks=outmat;
% ALLtracks:  obj no, frame, x, y, time, frame, length 
% MIGHT NEED TO CHANGE COLUMN NUMBERS, depending on output from Dhruv's
% (FluoreT)code
Null_Alltracks = [ALLtracks(:,1) ALLtracks(:,4) ALLtracks(:,3),...
    ALLtracks(:,6), ALLtracks(:,2)];
ALLtracks = Null_Alltracks; 
objno=unique(ALLtracks(:,1));
objwise=cell(1, max(objno));
frmwise=cell(1,Findex(end));
p=0;
w =waitbar(0, 'Saving movie..');
for j=Findex
    p=p+1;
    waitbar(p/numel(Findex));
    r=find(ALLtracks(:,5)==j);
    frmwise{j}=[ALLtracks(r,1:3),ALLtracks(r,5)]; %obj no, x, y, frame
    fiG= imread(imagename,j);
    fj= figure(j); 
    set(fj,'visible', 'off'), 
    imshow(fiG, 'Border', 'tight')
    
    for i=min(objno):max(objno)
        rr= find(frmwise{j}(:,1)==i);
        objwise{i}=[objwise{i}; frmwise{j}(rr,2:3)];% x,y
        
        if j==frmwise{j}(rr,4)
            set(gca, 'visible', 'off')
                hold on, plot(objwise{i}(1:end,1),objwise{i}(1:end,2),...
                '-r.',...
                'Linewidth', 1.2,...
                'MarkerSize', 2)

            hold on, plot(objwise{i}(end,1),objwise{i}(end,2),'.b','MarkerSize', 1)
            %hold on, text(objwise{i}(end,1),objwise{i}(end,2),sprintf('%i', i),...
                %'Color', 'y','FontSize', 0.0)
        else
            continue
        end
    end
    f=getframe(gcf); % converting frame to matrix
    imwrite(f.cdata, outfile, 'tif', 'Compression', 'none',...
        'WriteMode', 'append');%, 'Resolution', 1/Scaling_factor);
    delete(fj);
end
delete(w);
end