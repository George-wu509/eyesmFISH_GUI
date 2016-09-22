function varargout = eyesmFISH_GUI(varargin)
% EYESMFISH_GUI MATLAB code for eyesmFISH_GUI.fig
%      EYESMFISH_GUI, by itself, creates a new EYESMFISH_GUI or raises the existing
%      singleton*.
%
%      H = EYESMFISH_GUI returns the handle to a new EYESMFISH_GUI or the handle to
%      the existing singleton*.
%
%      EYESMFISH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EYESMFISH_GUI.M with the given input arguments.
%
%      EYESMFISH_GUI('Property','Value',...) creates a new EYESMFISH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eyesmFISH_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eyesmFISH_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eyesmFISH_GUI

% Last Modified by GUIDE v2.5 30-Aug-2016 04:48:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eyesmFISH_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @eyesmFISH_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end
function p=p_setting(handles,p,lsm_stack)

% 1.1 general p main setting
p.lsm_type = '*.lsm';
p.figure_tail_tif = '.tif';
if isfield(p,'chal_num')~=1
    p.chal_num=size(lsm_stack(1).data,2); % total number of image channel
end


% 1.2 rawImage parameters
p.savetiff=1;
p.savefinalresult=1;
p.rawImage16=1;
p.closeallfigure=1;
p.DAPI_chal=1;
p.FISH_chal=3;


% 1.3 gaussain filter parameters
p.gau_type='gaussian'; %gaussian filter type 'gaussian' or 'log'
p.gau_hsize=3;    %gaussian hsize
p.gau_sigma=1;    %gaussian sigma
p.init_stack =11;
p.final_stack = 20;


% 1.4 background estimation and substration
p.back_th_num = 100;


% 1.5 image filter parameters
p.open_size = 5;


% gaussian fitting parameters
p.remove_bou_back=0;
p.candi_start=3;
p.candi_end=8;
p.nearby_radius_max=200;



end

%% GUI unit callback functions
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles,out]=openfile(handles,hObject);
if out==1
    handles=image_conversion(handles,hObject);
end

% Update handles
guidata(hObject, handles);

end
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameter_setting;
end
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=eyesmFISH(handles);

% Update handles
guidata(hObject, handles);
end
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% show axes image
%set(handles.uibuttongroup2,'Visible','on');
%tiff_numshow=['/',num2str(size(tiff_images,2))];
%set(handles.text5,'String',tiff_numshow);

    load([handles.datafolder 'data.mat'],'p','tiff_images','max_image');
    %Va = get(hObject,'Value');
    %set(handles.slider1,'String',num2str(round(Va)));

    tiff_num=str2num(get(handles.slider1,'String'));%tifftxt=[num2str(tiff_num) ' stack'];
    tiff_nummax=size(tiff_images,2)+1;
    if tiff_num ~=tiff_nummax
        tifftxt=[num2str(tiff_num) ' stack'];
    else
        tifftxt='Maximage';
    end
    
    set(handles.text6,'String',tifftxt);
    tiff_chal_t=get(handles.popupmenu1,'Value');
    if tiff_chal_t==1
        tiff_chal=p.DAPI_chal;
    elseif tiff_chal_t==2
        tiff_chal=p.FISH_chal;
    end
    
    %myImage = tiff_images{1,tiff_num}(:,:,tiff_chal);
    if tiff_num ~=tiff_nummax
        myImage = tiff_images{1,tiff_num}(:,:,tiff_chal);
    else
        myImage = max_image(:,:,tiff_chal);
    end
    
    set(handles.axes1,'Units','pixels');
    resizePos = get(handles.axes1,'Position');
    myImage= imresize(myImage, [resizePos(3) resizePos(3)]);
    axes(handles.axes1);
    imagesc(myImage);
    set(handles.axes1,'Units','normalized');
    
    
    guidata(hObject, handles);
end
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    %if rem(Va,1)~=0
    %    set(hObject,'Value',round(Va))
    %end
    load([handles.datafolder 'data.mat'],'p','tiff_images','max_image');
    Va = get(hObject,'Value');
    set(handles.slider1,'String',num2str(round(Va)));

    tiff_num=str2num(get(handles.slider1,'String'));%tifftxt=[num2str(tiff_num) ' stack'];
    tiff_nummax=size(tiff_images,2)+1;
    if tiff_num ~=tiff_nummax
        tifftxt=[num2str(tiff_num) ' stack'];
    else
        tifftxt='Maximage';
    end
    
    set(handles.text6,'String',tifftxt);
    tiff_chal_t=get(handles.popupmenu1,'Value');
    if tiff_chal_t==1
        tiff_chal=p.DAPI_chal;
    elseif tiff_chal_t==2
        tiff_chal=p.FISH_chal;
    end
    
    %myImage = tiff_images{1,tiff_num}(:,:,tiff_chal);
    if tiff_num ~=tiff_nummax
        myImage = tiff_images{1,tiff_num}(:,:,tiff_chal);
    else
        myImage = max_image(:,:,tiff_chal);
    end
    
    set(handles.axes1,'Units','pixels');
    resizePos = get(handles.axes1,'Position');
    myImage= imresize(myImage, [resizePos(3) resizePos(3)]);
    axes(handles.axes1);
    imagesc(myImage);
    set(handles.axes1,'Units','normalized');
    
    
    guidata(hObject, handles);
end

%% Main eyesmFISH function
function [handles,out]=openfile(handles,hObject)

    % open Lsm image file
    if isempty(handles)==1
        handles = guihandles(gcf);
    end
    set(handles.edit1,'String','');pause(0.1)
    [FileName,folder_name] = uigetfile('*.lsm','Select the MATLAB code file');
    if FileName==0
        out=0;
    else
    handles.folder_name=folder_name;
    handles.FileName=FileName;

    % path setting
    handles.root_folder=cd;
    if ismac==1
        handles.functions_folder=[cd '/functions'];
    elseif ispc==1
        handles.functions_folder=[cd '\functions'];
    end
    addpath(genpath(handles.functions_folder));
    addpath(genpath(handles.folder_name));

    % LSM file loading/resave
    lsmfile=[handles.folder_name handles.FileName];
    % Lsm file loading
    set(handles.edit1,'String','loading.....');pause(0.1)
    try
        lsm_stack = tiffread(lsmfile);
    catch
        lsm_stack = readlsm(lsmfile);
    end
    % import parameters file and get image information
    p=[];
    p=p_setting(handles,p,lsm_stack);    
    
    lsmfolder=strrep(lsmfile,p.lsm_type(2:end),'');
    warning('off','all');
    if exist(lsmfolder,'dir') ~= 7
        if ismac||isunix==1
            mkdir(strrep(lsmfile,p.lsm_type(2:end),'/'));
        elseif ispc==1
            mkdir(strrep(lsmfile,p.lsm_type(2:end),'\'));
        end
    end
       
    % parameter file create
    if ismac||isunix==1
        datafolder=strrep(lsmfile,p.lsm_type(2:end),'/data/');
        handles.datafolder=datafolder;
        if exist(handles.datafolder,'dir') ~= 7
            mkdir(handles.datafolder);
        end
        save([handles.datafolder 'data.mat'],'p','lsm_stack','-v7.3');
        save([handles.functions_folder '/root.mat'],'datafolder');
    elseif ispc==1
        datafolder=strrep(lsmfile,p.lsm_type(2:end),'\data\');
        handles.datafolder=datafolder;
        if exist(handles.datafolder,'dir') ~= 7
            mkdir(handles.datafolder);
        end
        save([handles.datafolder 'data.mat'],'p','lsm_stack','-v7.3');
        save([handles.functions_folder '/root.mat'],'datafolder');
    end
    
    set(handles.edit1,'String',lsmfile);pause(0.1);
    set(handles.pushbutton2,'Enable','on');
    set(handles.pushbutton3,'Enable','on');
    set(handles.pushbutton4,'Enable','on');
    set(handles.slider1,'Value',1);
    guidata(hObject, handles);
    out=1;
    end
    
    function stack = tiffread(filename, indices)

% tiffread, version 2.9 May 10, 2010
%
% stack = tiffread;
% stack = tiffread(filename);
% stack = tiffread(filename, indices);
%
% Reads 8,16,32 bits uncompressed grayscale and (some) color tiff files,
% as well as stacks or multiple tiff images, for example those produced
% by metamorph, Zeiss LSM or NIH-image.
%
% The function can be called with a file name in the current directory,
% or without argument, in which case it pops up a file opening dialog
% to allow for a manual selection of the file.
% If the stacks contains multiples images, reading can be restricted by
% specifying the indices of the desired images (eg. 1:5), or just one index (eg. 2).
%
% The returned value 'stack' is a vector struct containing the images 
% and their meta-data. The length of the vector is the number of images.
% The image pixels values are stored in a field .data, which is a simple
% matrix for gray-scale images, or a cell-array of matrices for color images.
%
% The pixels values are returned in their native (usually integer) format,
% and must be converted to be used in most matlab functions.
%
% Example:
% im = tiffread('spindle.stk');
% imshow( double(im(5).data) );
%
% Only a fraction of the TIFF standard is supported, but you may extend support
% by modifying this file. If you do so, please return your modification to us,
% such that the added functionality can be redistributed to everyone.
%
% Copyright (C) 1999-2010 Francois Nedelec, 
% with contributions from:
%   Kendra Burbank for the waitbar
%   Hidenao Iwai for the code to read floating point images,
%   Stephen Lang to be more compliant with PlanarConfiguration
%   Jan-Ulrich Kreft for Zeiss LSM support
%   Elias Beauchanp and David Kolin for additional Metamorph support
%   Jean-Pierre Ghobril for requesting that image indices may be specified
%   Urs Utzinger for the better handling of color images, and LSM meta-data
%   O. Scott Sands for support of GeoTIFF tags
%   
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details:
% <http://www.gnu.org/licenses/>.
%
% Francois Nedelec
% nedelec (at) embl.de
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
% http://www.embl.org
% http://www.cytosim.org




%Optimization: join adjacent TIF strips: this results in faster reads
consolidateStrips = 1;

%without argument, we ask the user to choose a file:
if nargin < 1
    [filename, pathname] = uigetfile('*.tif;*.stk;*.lsm', 'select image file');
    filename = [ pathname, filename ];
end

if (nargin<=1);  indices = 1:10000; end


% not all valid tiff tags have been included, as they are really a lot...
% if needed, tags can easily be added to this code
% See the official list of tags:
% http://partners.adobe.com/asn/developer/pdfs/tn/TIFF6.pdf
%
% the structure IMG is returned to the user, while TIF is not.
% so tags usefull to the user should be stored as fields in IMG, while
% those used only internally can be stored in TIF.

global TIF;
TIF = [];

%counters for the number of images read and skipped
img_skip  = 0;
img_read  = 1;
hWaitbar  = [];

%% set defaults values :
TIF.SampleFormat     = 1;
TIF.SamplesPerPixel  = 1;
TIF.BOS              = 'ieee-le';          %byte order string

if  isempty(findstr(filename,'.'))
    filename = [filename,'.tif'];
end

TIF.file = fopen(filename,'r','l');
if TIF.file == -1
    stkname = strrep(filename, '.tif', '.stk');
    TIF.file = fopen(stkname,'r','l');
    if TIF.file == -1
        error(['File "',filename,'" not found.']);
    else
        filename = stkname;
    end
end
[s, m] = fileattrib(filename);

% obtain the full file path:
filename = m.Name;

% find the file size in bytes:
% m = dir(filename);
% filesize = m.bytes;


%% read header
% read byte order: II = little endian, MM = big endian
byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
    TIF.BOS = 'ieee-le';                                % Intel little-endian format
elseif ( strcmp(byte_order','MM') )
    TIF.BOS = 'ieee-be';
else
    error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies file as TIFF format
tiff_id = fread(TIF.file,1,'uint16', TIF.BOS);
if (tiff_id ~= 42)
    error('This is not a TIFF file (missing 42).');
end

%% ---- read the byte offset for the first image file directory (IFD)
TIF.img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);

while  TIF.img_pos ~= 0 

    clear IMG;
    IMG.filename = filename;
    % move in the file to the first IFD
    status = fseek(TIF.file, TIF.img_pos, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end

    %disp(strcat('reading img at pos :',num2str(TIF.img_pos)));

    %read in the number of IFD entries
    num_entries = fread(TIF.file,1,'uint16', TIF.BOS);
    %disp(strcat('num_entries =', num2str(num_entries)));

    %read and process each IFD entry
    for i = 1:num_entries

        % save the current position in the file
        file_pos  = ftell(TIF.file);

        % read entry tag
        TIF.entry_tag = fread(TIF.file, 1, 'uint16', TIF.BOS);
        % read entry
        entry = readIFDentry;
        %disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));

        switch TIF.entry_tag
            case 254
                TIF.NewSubfiletype = entry.val;
            case 256         % image width - number of column
                IMG.width          = entry.val;
            case 257         % image height - number of row
                IMG.height         = entry.val;
                TIF.ImageLength    = entry.val;
            case 258         % BitsPerSample per sample
                TIF.BitsPerSample  = entry.val;
                TIF.BytesPerSample = TIF.BitsPerSample / 8;
                IMG.bits           = TIF.BitsPerSample(1);
                %fprintf('BitsPerSample %i %i %i\n', entry.val);
            case 259         % compression
                if ( entry.val ~= 1 )
                    error(['Compression format ', num2str(entry.val),' not supported.']);
                end
            case 262         % photometric interpretation
                TIF.PhotometricInterpretation = entry.val;
                if ( TIF.PhotometricInterpretation == 3 )
                    warning('tiffread2:LookUp', 'Ignoring TIFF look-up table');
                end
            case 269
                IMG.document_name  = entry.val;
            case 270         % comments:
                IMG.info           = entry.val;
            case 271
                IMG.make           = entry.val;
            case 273         % strip offset
                if ~exist('previous_offset','var')
                    previous_offset = zeros(size(entry.val));
                    current_offset = zeros(size(entry.val));
                    wrap_offset = zeros(size(entry.val));
                end
                current_offset = entry.val;

                if ~exist('stack','var') || (IMG.width == stack(1).width && IMG.height == stack(1).height)
                    wrap_offset = wrap_offset+2^32.*(current_offset < previous_offset);
                    previous_offset = entry.val;
                    TIF.StripOffsets   = current_offset+wrap_offset;
                else
                    TIF.StripOffsets   = current_offset;
                end

                TIF.StripNumber    = entry.cnt;
                %fprintf('StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
            case 277         % sample_per pixel
                TIF.SamplesPerPixel  = entry.val;
                %fprintf('Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
            case 278         % rows per strip
                TIF.RowsPerStrip   = entry.val;
            case 279         % strip byte counts - number of bytes in each strip after any compressio
                TIF.StripByteCounts= entry.val;
            case 282         % X resolution
                IMG.x_resolution   = entry.val;
            case 283         % Y resolution
                IMG.y_resolution   = entry.val;
            case 284         %planar configuration describe the order of RGB
                TIF.PlanarConfiguration = entry.val;
            case 296         % resolution unit
                IMG.resolution_unit= entry.val;
            case 305         % software
                IMG.software       = entry.val;
            case 306         % datetime
                IMG.datetime       = entry.val;
            case 315
                IMG.artist         = entry.val;
            case 317        %predictor for compression
                if (entry.val ~= 1); error('unsuported predictor value'); end
            case 320         % color map
                IMG.cmap           = entry.val;
                IMG.colors         = entry.cnt/3;
            case 339
                TIF.SampleFormat   = entry.val;
            case 33550       % GeoTIFF ModelPixelScaleTag
                IMG.ModelPixelScaleTag    = entry.val;
            case 33628       %metamorph specific data
                IMG.MM_private1    = entry.val;
            case 33629       %this tag identify the image as a Metamorph stack!
                TIF.MM_stack       = entry.val;
                TIF.MM_stackCnt    = entry.cnt;
            case 33630       %metamorph stack data: wavelength
                TIF.MM_wavelength  = entry.val;
            case 33631       %metamorph stack data: gain/background?
                TIF.MM_private2    = entry.val;
            case 33922       % GeoTIFF ModelTiePointTag
                IMG.ModelTiePointTag    = entry.val;
            case 34412       % Zeiss LSM data
                LSM_info           = entry.val;
            case 34735       % GeoTIFF GeoKeyDirectory
                IMG.GeoKeyDirTag       = entry.val;
            case 34737       % GeoTIFF GeoASCIIParameters
                IMG.GeoASCII       = entry.val;
            case 42113       % GeoTIFF GDAL_NODATA
                IMG.GDAL_NODATA    = entry.val;
            otherwise
                fprintf( 'Ignored TIFF entry with tag %i (cnt %i)\n', TIF.entry_tag, entry.cnt);
        end
        
        % calculate bounding box  if you've got the stuff
        if isfield(IMG, 'ModelPixelScaleTag') && isfield(IMG, 'ModelTiePointTag') && isfield(IMG, 'height')&& isfield(IMG, 'width'),
            IMG.North=IMG.ModelTiePointTag(5)-IMG.ModelPixelScaleTag(2)*IMG.ModelTiePointTag(2);
            IMG.South=IMG.North-IMG.height*IMG.ModelPixelScaleTag(2);
            IMG.West=IMG.ModelTiePointTag(4)+IMG.ModelPixelScaleTag(1)*IMG.ModelTiePointTag(1);
            IMG.East=IMG.West+IMG.width*IMG.ModelPixelScaleTag(1);
        end

        % move to next IFD entry in the file
        status = fseek(TIF.file, file_pos+12, -1);
        if status == -1
            error('invalid file offset (error on fseek)');
        end
    end

    %Planar configuration is not fully supported
    %Per tiff spec 6.0 PlanarConfiguration irrelevent if SamplesPerPixel==1
    %Contributed by Stephen Lang
    if (TIF.SamplesPerPixel ~= 1) && ( ~isfield(TIF, 'PlanarConfiguration') || TIF.PlanarConfiguration == 1 )
        error('PlanarConfiguration = 1 is not supported');
    end

    %total number of bytes per image:
    PlaneBytesCnt = IMG.width * IMG.height * TIF.BytesPerSample;

    %% try to consolidate the TIFF strips if possible
    
    if consolidateStrips
        %Try to consolidate the strips into a single one to speed-up reading:
        BytesCnt = TIF.StripByteCounts(1);

        if BytesCnt < PlaneBytesCnt

            ConsolidateCnt = 1;
            %Count how many Strip are needed to produce a plane
            while TIF.StripOffsets(1) + BytesCnt == TIF.StripOffsets(ConsolidateCnt+1)
                ConsolidateCnt = ConsolidateCnt + 1;
                BytesCnt = BytesCnt + TIF.StripByteCounts(ConsolidateCnt);
                if ( BytesCnt >= PlaneBytesCnt ); break; end
            end

            %Consolidate the Strips
            if ( BytesCnt <= PlaneBytesCnt(1) ) && ( ConsolidateCnt > 1 )
                %fprintf('Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
                TIF.StripByteCounts = [BytesCnt; TIF.StripByteCounts(ConsolidateCnt+1:TIF.StripNumber ) ];
                TIF.StripOffsets = TIF.StripOffsets( [1 , ConsolidateCnt+1:TIF.StripNumber] );
                TIF.StripNumber  = 1 + TIF.StripNumber - ConsolidateCnt;
            end
        end
    end

    %% read the next IFD address:
    TIF.img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);
    %if (TIF.img_pos) disp(['next ifd at', num2str(TIF.img_pos)]); end

    if isfield( TIF, 'MM_stack' )

        sel = ( indices <= TIF.MM_stackCnt );
        indices = indices(sel);
        
        if numel(indices) > 1
            hWaitbar = waitbar(0,'Reading images...','Name','TiffRead');
        end

        %this loop reads metamorph stacks:
        for ii = indices

            TIF.StripCnt = 1;
            offset = PlaneBytesCnt * (ii-1);

            %read the image channels
            for c = 1:TIF.SamplesPerPixel
                IMG.data{c} = read_plane(offset, IMG.width, IMG.height, c);
            end

            % print a text timer on the main window, or update the waitbar
            % fprintf('img_read %i img_skip %i\n', img_read, img_skip);
            if ~isempty( hWaitbar )
                waitbar(img_read/numel(indices), hWaitbar);
            end
            
            [ IMG.MM_stack, IMG.MM_wavelength, IMG.MM_private2 ] = splitMetamorph(ii);
            
            stack(img_read) = IMG;
            img_read = img_read + 1;

        end
        break;

    else

        %this part reads a normal TIFF stack:
        
        read_img = any( img_skip+img_read == indices );
        if exist('stack','var')
            if IMG.width ~= stack(1).width || IMG.height ~= stack(1).height
                %setting read_it=0 will skip dissimilar images:
                %comment-out the line below to allow dissimilar stacks
                read_img = 0;
            end
        end
        
        if read_img
            TIF.StripCnt = 1;
            %read the image channels
            for c = 1:TIF.SamplesPerPixel
                IMG.data{c} = read_plane(0, IMG.width, IMG.height, c);
            end

            try
                stack(img_read) = IMG;  % = orderfields(IMG);
                img_read = img_read + 1;
            catch
                fprintf('Tiffread skipped dissimilar image %i\n', img_read+img_skip);
                img_skip = img_skip + 1;
             end
            
            if  all( img_skip+img_read > indices )
                break;
            end

        else
            img_skip = img_skip + 1;
        end

    end
end

%% remove the cell structure if there is always only one channel
flat = 1;
for i = 1:numel(stack)
    if numel(stack(i).data) ~= 1
        flat = 0;
        break;
    end
end

if flat
    for i = 1:numel(stack)
        stack(i).data = stack(i).data{1};
    end
end


%% distribute the MetaMorph info
if isfield(TIF, 'MM_stack') && isfield(IMG, 'info') && ~isempty(IMG.info)
    MM = parseMetamorphInfo(IMG.info, TIF.MM_stackCnt);
    for i = 1:numel(stack)
        stack(i).MM = MM(i);
    end
end

%% duplicate the LSM info
if exist('LSM_info', 'var')
    for i = 1:numel(stack)
        stack(i).lsm = LSM_info;
    end
end


%% return

if ~ exist('stack', 'var')
    stack = [];
end

%clean-up
fclose(TIF.file);
if ~isempty( hWaitbar )
    delete( hWaitbar );
end


end
    function plane = read_plane(offset, width, height, plane_nb)

global TIF;

%return an empty array if the sample format has zero bits
if ( TIF.BitsPerSample(plane_nb) == 0 )
    plane=[];
    return;
end

%fprintf('reading plane %i size %i %i\n', plane_nb, width, height);

%determine the type needed to store the pixel values:
switch( TIF.SampleFormat )
    case 1
        classname = sprintf('uint%i', TIF.BitsPerSample(plane_nb));
    case 2
        classname = sprintf('int%i', TIF.BitsPerSample(plane_nb));
    case 3
        if ( TIF.BitsPerSample(plane_nb) == 32 )
            classname = 'single';
        else
            classname = 'double';
        end
    otherwise
        error('unsuported TIFF sample format %i', TIF.SampleFormat);
end

% Preallocate a matrix to hold the sample data:
try
    plane = zeros(width, height, classname);
catch
    %compatibility with older matlab versions:
    eval(['plane = ', classname, '(zeros(width, height));']);
end

% Read the strips and concatenate them:
line = 1;
while ( TIF.StripCnt <= TIF.StripNumber )

    strip = read_strip(offset, width, plane_nb, TIF.StripCnt, classname);
    TIF.StripCnt = TIF.StripCnt + 1;

    % copy the strip onto the data
    plane(:, line:(line+size(strip,2)-1)) = strip;

    line = line + size(strip,2);
    if ( line > height )
        break;
    end

end

% Extract valid part of data if needed
if ~all(size(plane) == [width height]),
    plane = plane(1:width, 1:height);
    warning('tiffread2:Crop','Cropping data: found more bytes than needed');
end

% transpose the image (otherwise display is rotated in matlab)
plane = plane';

    end
    function strip = read_strip(offset, width, plane_nb, stripCnt, classname)

global TIF;

%fprintf('reading strip at position %i\n',TIF.StripOffsets(stripCnt) + offset);
StripLength = TIF.StripByteCounts(stripCnt) ./ TIF.BytesPerSample(plane_nb);

%fprintf( 'reading strip %i\n', stripCnt);
status = fseek(TIF.file, TIF.StripOffsets(stripCnt) + offset, 'bof');
if status == -1
    error('invalid file offset (error on fseek)');
end

bytes = fread( TIF.file, StripLength, classname, TIF.BOS );

if any( length(bytes) ~= StripLength )
    error('End of file reached unexpectedly.');
end

strip = reshape(bytes, width, StripLength / width);

    end
    function [nbBytes, matlabType] = convertType(tiffType)
switch (tiffType)
    case 1
        nbBytes=1;
        matlabType='uint8';
    case 2
        nbBytes=1;
        matlabType='uchar';
    case 3
        nbBytes=2;
        matlabType='uint16';
    case 4
        nbBytes=4;
        matlabType='uint32';
    case 5
        nbBytes=8;
        matlabType='uint32';
    case 7
        nbBytes=1;
        matlabType='uchar';
    case 11
        nbBytes=4;
        matlabType='float32';
    case 12
        nbBytes=8;
        matlabType='float64';
    otherwise
        error('tiff type %i not supported', tiffType)
end
    end
    function  entry = readIFDentry()

global TIF;
entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.BOS);
entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.BOS);
%disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);

[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32', TIF.BOS);
    %disp(strcat('offset = ', num2str(offset)));
    status = fseek(TIF.file, offset, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end

end


if TIF.entry_tag == 33629   % metamorph 'rationals'
    entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.BOS);
elseif TIF.entry_tag == 34412  %TIF_CZ_LSMINFO
    entry.val = readLSMinfo;
else
    if entry.tiffType == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);
    else
        entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.BOS);
    end
end

if ( entry.tiffType == 2 );
    entry.val = char(entry.val');
end

    end
    function [MMstack, MMwavelength, MMprivate2] = splitMetamorph(imgCnt)

global TIF;

MMstack = [];
MMwavelength = [];
MMprivate2 = [];

if TIF.MM_stackCnt == 1
    return;
end

left  = imgCnt - 1;

if isfield( TIF, 'MM_stack' )
    S = length(TIF.MM_stack) / TIF.MM_stackCnt;
    MMstack = TIF.MM_stack(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_wavelength' )
    S = length(TIF.MM_wavelength) / TIF.MM_stackCnt;
    MMwavelength = TIF.MM_wavelength(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_private2' )
    S = length(TIF.MM_private2) / TIF.MM_stackCnt;
    MMprivate2 = TIF.MM_private2(S*left+1:S*left+S);
end

    end
    function mm = parseMetamorphInfo(info, cnt)

info   = regexprep(info, '\r\n|\o0', '\n');
parse  = textscan(info, '%s %s', 'Delimiter', ':');
tokens = parse{1};
values = parse{2};

first = char(tokens(1,1));

k = 0;
mm = struct('Exposure', zeros(cnt,1));
for i=1:size(tokens,1)
    tok = char(tokens(i,1));
    val = char(values(i,1));
    %fprintf( '"%s" : "%s"\n', tok, val);
    if strcmp(tok, first)
        k = k + 1;
    end
    if strcmp(tok, 'Exposure')
        [v, c, e, pos] = sscanf(val, '%i');
        unit = val(pos:length(val));
        %return the exposure in milli-seconds
        switch( unit )
            case 'ms'
                mm(k).Exposure = v;
            case 's'
                mm(k).Exposure = v * 1000;
            otherwise
                warning('tiffread2:Unit', ['Exposure unit "',unit,'" not recognized']);
                mm(k).Exposure = v;
        end
    else
        switch tok
            case 'Binning'
                % Binning: 1 x 1 -> [1 1]
                mm(k).Binning = sscanf(val, '%d x %d')';
            case 'Region'
                mm(k).Region = sscanf(val, '%d x %d, offset at (%d, %d)')';
            otherwise
                field  = regexprep(tok, ' ', '');
                if strcmp(val, 'Off')
                    eval(['mm(k).',field,'=0;']);
                elseif strcmp(val, 'On')
                    eval(['mm(k).',field,'=1;']);
                elseif isstrprop(val,'digit')
                    eval(['mm(k).',field,'=str2num(val)'';']);
                else
                    eval(['mm(k).',field,'=val;']);
                end
        end
    end
end

    end
    function R = readLSMinfo()

% Read part of the LSM info table version 2
% this provides only very partial information, since the offset indicate that
% additional data is stored in the file
global TIF;

R.MagicNumber            = sprintf('0x%09X',fread(TIF.file, 1, 'uint32', TIF.BOS));
StructureSize          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionX             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionY             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionZ             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionChannels      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionTime          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.IntensityDataType      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ThumbnailX             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ThumbnailY             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.VoxelSizeX             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeY             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeZ             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginX                = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginY                = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginZ                = fread(TIF.file, 1, 'float64', TIF.BOS);
R.ScanType               = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.SpectralScan           = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.DataType               = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetVectorOverlay    = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetInputLut         = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetOutputLut        = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetChannelColors    = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.TimeInterval           = fread(TIF.file, 1, 'float64', TIF.BOS);
OffsetChannelDataTypes = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetScanInformation  = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetKsData           = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetTimeStamps       = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetEventList        = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetRoi              = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetBleachRoi        = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetNextRecording    = fread(TIF.file, 1, 'uint32', TIF.BOS);

% There are more information stored in this table, which is not read here


%read real acquisition times:
if ( OffsetTimeStamps > 0 )
    
    status = fseek(TIF.file, OffsetTimeStamps, -1);
    if status == -1
        error('error on fseek');
    end
    
    StructureSize          = fread(TIF.file, 1, 'int32', TIF.BOS);
    NumberTimeStamps       = fread(TIF.file, 1, 'int32', TIF.BOS);
    for i=1:NumberTimeStamps
        R.TimeStamp(i)       = fread(TIF.file, 1, 'float64', TIF.BOS);
    end
    
    %calculate elapsed time from first acquisition:
    R.TimeOffset = R.TimeStamp - R.TimeStamp(1);
    
end


    end

end
function handles=image_conversion(handles,hObject)

    % Load image information: lsm_stack
    %guidata(hObject, handles);
    set(handles.edit1,'String','running.....');pause(0.1)
    load([handles.datafolder 'data.mat'],'p','lsm_stack');
    tiff_num=size(lsm_stack,2); % total number of image pages
    if lsm_stack(1).bits==8&&p.rawImage16==1
        lsm_stack=convert8to16lsm(lsm_stack);   % convert unit8 lsm into unit16 lsm
    end

    % Transfrom lsm image into tif files
    lsmfile=[handles.folder_name handles.FileName];
    if ismac||isunix==1
        tiffolder=strrep(lsmfile,p.lsm_type(2:end),'/tif/');
    elseif ispc==1
        tiffolder=strrep(lsmfile,p.lsm_type(2:end),'\tif\');
    end
    if exist(tiffolder,'dir') ~= 7
        mkdir(tiffolder);
    end
    handles.tiffolder=tiffolder;
    for I_layer = 1:tiff_num
        stack_raw = lsm_stack(I_layer).data;    % stack_raw: 4 channel tiff images
        tiff_image = zeros(0);
        if iscell(stack_raw)
            for I_color = 1:length(stack_raw)
                tiff_image = cat(3,tiff_image,stack_raw{I_color});
            end
        else
            tiff_image = stack_raw;
        end

        if size(tiff_image,3) == 1
            tiff_image(:,:,2) = tiff_image(:,:,1);
        end
        if size(tiff_image,3) == 2
            tiff_image(:,:,3) = tiff_image(:,:,2);
        end
        
        % output tif files
        tiff_images{I_layer}=tiff_image;
        if p.savetiff==1
            imwrite(tiff_image,[tiffolder 'stack' num2str(I_layer) '.tif']);
        end
    end
    
    % calculate max iamges
    max_image = maximage(tiff_images,p);
    
    % show axes image
    set(handles.uibuttongroup2,'Visible','on');
    tiff_nummax=size(tiff_images,2)+1;
       
    set(handles.slider1,'Max',tiff_nummax);set(handles.slider1,'Min',1);    
    tiff_num=get(handles.slider1,'Value');
    if tiff_num ~=tiff_nummax
        tifftxt=[num2str(tiff_num) ' stack'];
    else
        tifftxt='Maximage';
    end
    set(handles.slider1, 'SliderStep', [1/(tiff_nummax-1) , 1/(tiff_nummax-1) ]);
    set(handles.slider1,'String',num2str(1));
    guidata(hObject, handles);
    
    set(handles.text6,'String',tifftxt);
    tiff_chal_t=get(handles.popupmenu1,'Value');
    if tiff_chal_t==1
        tiff_chal=p.DAPI_chal;
    elseif tiff_chal_t==2
        tiff_chal=p.FISH_chal;
    end
    
    if tiff_num ~=tiff_nummax
        myImage = tiff_images{1,tiff_num}(:,:,tiff_chal);
    else
        myImage = max_image(:,:,tiff_chal);
    end
    
    set(handles.axes1,'Units','pixels');
    resizePos = get(handles.axes1,'Position');
    myImage= imresize(myImage, [resizePos(3) resizePos(3)]);
    axes(handles.axes1);
    imagesc(myImage);
    set(handles.axes1,'Units','normalized');
    
    p_original = p;
    save([handles.datafolder 'data.mat'],'tiff_images','max_image','p','p_original','-v7.3');
    set(handles.edit1,'String','image conversion FINISH!!');pause(0.1)
    
    function lsm_new=convert8to16lsm(lsm_old)
    [~,nmax]=size(lsm_old);
    [~,chal_max]=size(lsm_old(1).data);
    lsm_new=lsm_old;
    for ni=1:nmax
        lsm_new(ni).bits=16;
        for chal=1:chal_max
            lsm_new(ni).data{1,chal} = im2uint16(lsm_old(ni).data{1,chal});
        end
    end
    end 
end
function handles=eyesmFISH(handles)

% load data

set(handles.edit1,'String','running.....');pause(0.1)
load([handles.datafolder 'data.mat'],'p','tiff_images','max_image');
%p=p_setting(handles,p);
tiffolder=handles.tiffolder;
tiff_num=str2num(get(handles.slider1,'String'));

%max_image
if tiff_num<=size(tiff_images,2)
    fishimage(tiff_images,tiff_num,max_image,p,handles);
else
    fishimage(tiff_images,0,max_image,p,handles);
end

end
function fishimage(tiff_im,tiff_num,max_image,p,handles)
% tiff_im = 1x29 cell[512x512x3]: all z-stacks image with all channels
% tiff_num = N:select stack numer  or maximage=0
% max_image: [512x512x3]:  maximage from selected z-stacks images


% ==== Step 0. Raw images collect ====
if tiff_num==0
    %raw_dapi = max_image(:,:,p.DAPI_chal);
    raw_fish = max_image(:,:,p.FISH_chal);
    raw_maxdapi = max_image(:,:,p.DAPI_chal);
    %raw_maxfish = max_image(:,:,p.FISH_chal);
else
    %raw_dapi = tiff_im{1,tiff_num}(:,:,p.DAPI_chal);
    raw_fish = tiff_im{1,tiff_num}(:,:,p.FISH_chal);
    raw_maxdapi = max_image(:,:,p.DAPI_chal);
    %raw_maxfish = max_image(:,:,p.FISH_chal);
end


% ===== Step 1. Identify background and nucleus masks =====
H = fspecial(p.gau_type,p.gau_hsize,p.gau_sigma);
im1 = imfilter(raw_fish,H);
im1(im1<0) = 0;
[propA,back_mask,nucle_mask,overlay1,back_blob_overlay]=conu_theh(im1,raw_maxdapi,p,handles);  % use raw_maxdapi to identify masks
if p.closeallfigure~=1
    figure;surf(raw_fish,'Edgecolor','none');colorbar;title('1. raw fish image');view(0, -90);xlim([0 512]);ylim([0 512]);
end



% ===== Step 2. Gaussian smooth and image background substraction=====
[raw_fishall,raw_fish]=image_background(tiff_im,tiff_num,max_image,back_mask,p);
H = fspecial(p.gau_type,p.gau_hsize,p.gau_sigma);
im1 = imfilter(raw_fish,H);
im1(im1<0) = 0;
im1all=zeros(size(raw_fishall));
for i=1:size(raw_fishall,3)
    im1all(:,:,i) = max(imfilter(raw_fishall(:,:,i),H),0);
end
if p.closeallfigure~=1
    figure;surf(raw_fish,'EdgeColor','none');colorbar;title('2. after gaussian filter and subtract background: im1');view(0, -90);xlim([0 512]);ylim([0 512]);
end

% ===== Step 3. Ommatadia identification algorithm =====
STATS_nuc = regionprops(nucle_mask,'Centroid','Area','Solidity','BoundingBox','Image');
STATS_back = regionprops(back_mask,'Centroid','Area','Solidity','BoundingBox','Image');
nucle_blob_overlay = make_overlaymap(STATS_nuc,raw_fish);

%
[overlay2,bn_matrix_est,back_blob_overlay,nucle_blob,back_mask] = nuc_back_integrate(p,STATS_back,STATS_nuc,raw_fish,back_mask,nucle_blob_overlay,back_blob_overlay);
% = output =
%if back# = N1, nucle# = N2, n nucle in each back
% overlay1: [512x512x3 unit8] overlay show
% bn_matrix_est: [1xN1 cell](nx6 double) = [nucle no, x, y, Area, Solidity, distance to back]
% back_blob_overlay: [1xN1 cell](512x512 log): nucle image
% nucle_blob_overlay: [1xN2 cell](512x512 log): nucle image
% nucle_blob: [N2x7 double]: [nucle no, x,y, Area, Solidity, belong to back, distance to back]
% back_mask: [512x512 log]: new back_mask only target region

%{
% back_area matrix
back_area=zeros(size(STATS_back,1),5);
for i=1:size(STATS_back,1)
    back_area(i,1) = i;
    back_area(i,2) = STATS_back(i,1).Centroid(1);
    back_area(i,3) = STATS_back(i,1).Centroid(2);
    back_area(i,4) = STATS_back(i,1).Area;
    back_area(i,5) = STATS_back(i,1).Solidity;
end

% nucle_lob matrix
nucle_blob=zeros(size(STATS_nuc,1),7);  %nucle_blob=[no., x, y, A, solidity, Backno., distance]
for i=1:size(STATS_nuc,1)
    nucle_blob(i,1) = i;
    nucle_blob(i,2) = STATS_nuc(i,1).Centroid(1);
    nucle_blob(i,3) = STATS_nuc(i,1).Centroid(2);
    nucle_blob(i,4) = STATS_nuc(i,1).Area;
    nucle_blob(i,5) = STATS_nuc(i,1).Solidity;
end

% Sort back_area and find corresponding nucle
back_area = sortrows(back_area,-4);
for ii=1:size(back_area,1)
    bn_dist = zeros(size(nucle_blob,1),1);
    for jj=1:size(nucle_blob,1)
        bn_dist(jj,1) = sqrt((nucle_blob(jj,2)-back_area(ii,2))^2+(nucle_blob(jj,3)-back_area(ii,3))^2);
    end
    % find nearest nuclear blobs
    bn_matrix = [nucle_blob(:,1:5) bn_dist]; %bn_matrix=[no., x, y, A, solidity, distance]
    bn_matrix = sortrows(bn_matrix,6);
    d=diff(bn_matrix(:,6));
    bn_matrix_est = bn_matrix(1:min(find(d>20)),:);
    
    for jj=1:size(bn_matrix_est,1)
        d=bn_matrix_est(jj,6);
        if d<nucle_blob(bn_matrix_est(jj,1),7)||nucle_blob(bn_matrix_est(jj,1),7)==0
            nucle_blob(bn_matrix_est(jj,1),6)=ii;
            nucle_blob(bn_matrix_est(jj,1),7)=d;
        end
    end
    
end

%if p.closeallfigure~=1
    
    % calculate overlay figures
    nucle_blob2 = sortrows(nucle_blob,6);
    %nucle_blob2=nucle_blob2(nucle_blob2(:,6)>0,:);
    ma=max(nucle_blob2(:,6));
    ma_color=rand(ma,3);
    overlay1 = imoverlay(raw_fish./max(max(raw_fish)),bwperim(back_mask),[0,1,0]); % background
    for oo=1:size(nucle_blob2,1)
        if nucle_blob2(oo,6)==0
            overlay1 = imoverlay(overlay1,bwperim(nucle_blob_overlay{1,nucle_blob2(oo,1)}),[0,0,1]); % nucleus
        else
            overlay1 = imoverlay(overlay1,nucle_blob_overlay{1,nucle_blob2(oo,1)},ma_color(nucle_blob2(oo,6),:)); % nucleus
        end
    end
%}    
    %figure;
    axes(handles.axes2);
    imshow(overlay2);pause(0.1);

%end


% ===== Step 4. nuclear region identification algorithm =====
%{
overlay3 = overlay2;
for i=1:size(bn_matrix_est,2)
    temp=false(size(back_blob_overlay{1,1}));
    for j=1:size(bn_matrix_est{1,i},1)
        % draw lines on temp image
    end
    overlay3 = imoverlay(overlay3,temp,[1,1,1]); % nucleus
end
%}

% ===== Step 5. Advanced Object detection using multi z-stacks =====

zmax=size(tiff_im,2);
for z=1:zmax
    [propA,back_mask,nucle_mask,overlay1,back_blob_overlay]=conu_theh(tiff_im{1,z}(:,:,p.FISH_chal),tiff_im{1,z}(:,:,p.DAPI_chal),p,handles);
end




% ===== Step 5. Object detection algorithm =====
im_object=im3;
im_object_n=im_object;
%im_object_n(30:220,150:360)=0;
newBoxPolygon=ObjectDetection(p,im_object_n,im_object(30:220,150:360));
%newBoxPolygon=ObjectDetection(im_object_n,im_object(160:220,180:230));


% ===== Step 6. Waterahed inside target =====
im_targets=im_object;
%figure;imagesc(im5);colorbar;title('6. im5 = im4 > 5: im5')
nucleus_prop = regionprops(im_targets,'Area','Centroid','Image','BoundingBox');
nucleus_center=zeros(size(nucleus_prop,1),6);
for n=1:size(nucleus_prop,1)
    nucleus_center(n,1:2)=nucleus_prop(n,1).Centroid;
    nucleus_center(n,3:4)=[nucleus_prop(n,1).BoundingBox(1)+0.5,nucleus_prop(n,1).BoundingBox(2)+0.5];
    nucleus_center(n,5:6)=[nucleus_center(n,3)+nucleus_prop(n,1).BoundingBox(3)-1,nucleus_center(n,4)+nucleus_prop(n,1).BoundingBox(4)-1];
end
for newB=1:size(newBoxPolygon,1)
    [nucleus_prop,new_perim]=DivideTargetRegion(raw_im,nucleus_prop,nucleus_center,newBoxPolygon{newB});
end

%D = bwdist(im3);D=-D;D(~im3)=-Inf;
%L=watershed(D);rgb=label2rgb(L);
%figure;imshow(rgb);


% ===== Step 7. Gaussian fitting interaged intensity value =====
[raw_parameter,xyI,ex_cell] = ptrack_par_nogauss(p,double(raw_im),nucleus_prop);
% xyI: have 1~(n)~N peak area
% xyI{n,1} = [fitting intensity, real intensity, x,y,z, xcenter, ycenter]
% xyI{n,2} = Xdatapoint, xyI{n,3} = Ydatapoint, xyI{n,4} = fitting intensity, xyI{n,5} = real intensity
% raw_parameter = [Ip, a, c, 2b, Ib, x0, y0, sigma1, sigma2, total fluorescence intensity, spot#]

%figure;surf(im3,'EdgeColor','none');colorbar;title('after watershed');view(0, 90);xlim([0 512]);ylim([0 512]);title('4. after imopen');



% ===== Step 8. Final overlay figures =====
overlay = imoverlay(raw_im./max(max(raw_im)),bwperim(im_targets),[0,0,1]); % nucleus
for se=1:size(new_perim,2)
    overlay = imoverlay(overlay,new_perim{se},[0,1/size(new_perim,2)*se,0]); % nucleus
end

%figure;
axes(handles.axes2);
imshow(overlay);hold on;
for j=1:size(newBoxPolygon,1)
    line(newBoxPolygon{j,1}(:, 1), newBoxPolygon{j,1}(:, 2), 'Color', 'y');hold on
end
hold off
for i=1:size(raw_parameter,1)
    text(raw_parameter(i,7),raw_parameter(i,6),num2str(round(raw_parameter(i,11))),'Color','white','FontSize',11);
end



end
function fishimage_old(tiff_im,tiff_num,max_image,p,handles)
% tiff_im = 1x29 cell[512x512x3]: all z-stacks image with all channels
% tiff_num = N:select stack numer  or maximage=0
% max_image: [512x512x3]:  maximage from selected z-stacks images


% ==== Step 0. Raw images collect ====
if tiff_num==0
    %raw_dapi = max_image(:,:,p.DAPI_chal);
    raw_fish = max_image(:,:,p.FISH_chal);
    raw_maxdapi = max_image(:,:,p.DAPI_chal);
    %raw_maxfish = max_image(:,:,p.FISH_chal);
else
    %raw_dapi = tiff_im{1,tiff_num}(:,:,p.DAPI_chal);
    raw_fish = tiff_im{1,tiff_num}(:,:,p.FISH_chal);
    raw_maxdapi = max_image(:,:,p.DAPI_chal);
    %raw_maxfish = max_image(:,:,p.FISH_chal);
end


% ===== Step 1. Identify background and nucleus masks =====
[propA,back_mask,nucle_mask,overlay1]=conu_theh(raw_fish,raw_maxdapi,p,handles);  % use raw_maxdapi to identify masks
if p.closeallfigure~=1
    figure;surf(raw_fish,'Edgecolor','none');colorbar;title('1. raw fish image');view(0, -90);xlim([0 512]);ylim([0 512]);
end



% ===== Step 2. Gaussian smooth and image background substraction=====
[raw_fishall,raw_fish]=image_background(tiff_im,tiff_num,max_image,back_mask,p);
H = fspecial(p.gau_type,p.gau_hsize,p.gau_sigma);
im1 = imfilter(raw_fish,H);
im1(im1<0) = 0;
im1all=zeros(size(raw_fishall));
for i=1:size(raw_fishall,3)
    im1all(:,:,i) = max(imfilter(raw_fishall(:,:,i),H),0);
end
if p.closeallfigure~=1
    figure;surf(raw_fish,'EdgeColor','none');colorbar;title('2. after gaussian filter and subtract background: im1');view(0, -90);xlim([0 512]);ylim([0 512]);
end

% ===== Step 3. Morphological filter =====
p.open_size=5;p.close_size=10;

se1 = strel('disk', p.open_size);
im2 = imopen(im1, se1);
if p.closeallfigure~=1
    figure;imagesc(im2);colorbar;title('3. Open filter(im1):  im2');
end

se2 = strel('disk', p.close_size);
im3 = imerode(im2, se2);
if p.closeallfigure~=1
    figure;imagesc(im3);colorbar;title('4. Close filter(im2):  im3');
end

%im3 = watershed(im1);
%g_prop = regionprops(im2,'Area'); 
%se1 = strel('disk', 5);
%im3 = imopen(im1, se1);
%figure;imagesc(im3);colorbar;title('4. after imopen: im1 --> im3');


% ===== Step 4. imfill calculation =====
im2=im2>0;
im3 = imfill(im2,'holes');
%se2 = strel('disk', 40);
%ie = imdilate(im3, se2);
%im4 = imreconstruct(ie, im3);
if p.closeallfigure~=1
    figure;imagesc(im3);title('4. imfill(im2): im3');
end
%im4 = imregionalmax(im3);


% ===== Step 5. Object detection algorithm =====
im_object=im3;
im_object_n=im_object;
%im_object_n(30:220,150:360)=0;
newBoxPolygon=ObjectDetection(p,im_object_n,im_object(30:220,150:360));
%newBoxPolygon=ObjectDetection(im_object_n,im_object(160:220,180:230));


% ===== Step 6. Waterahed inside target =====
im_targets=im_object;
%figure;imagesc(im5);colorbar;title('6. im5 = im4 > 5: im5')
nucleus_prop = regionprops(im_targets,'Area','Centroid','Image','BoundingBox');
nucleus_center=zeros(size(nucleus_prop,1),6);
for n=1:size(nucleus_prop,1)
    nucleus_center(n,1:2)=nucleus_prop(n,1).Centroid;
    nucleus_center(n,3:4)=[nucleus_prop(n,1).BoundingBox(1)+0.5,nucleus_prop(n,1).BoundingBox(2)+0.5];
    nucleus_center(n,5:6)=[nucleus_center(n,3)+nucleus_prop(n,1).BoundingBox(3)-1,nucleus_center(n,4)+nucleus_prop(n,1).BoundingBox(4)-1];
end
for newB=1:size(newBoxPolygon,1)
    [nucleus_prop,new_perim]=DivideTargetRegion(raw_im,nucleus_prop,nucleus_center,newBoxPolygon{newB});
end

%D = bwdist(im3);D=-D;D(~im3)=-Inf;
%L=watershed(D);rgb=label2rgb(L);
%figure;imshow(rgb);


% ===== Step 7. Gaussian fitting interaged intensity value =====
[raw_parameter,xyI,ex_cell] = ptrack_par_nogauss(p,double(raw_im),nucleus_prop);
% xyI: have 1~(n)~N peak area
% xyI{n,1} = [fitting intensity, real intensity, x,y,z, xcenter, ycenter]
% xyI{n,2} = Xdatapoint, xyI{n,3} = Ydatapoint, xyI{n,4} = fitting intensity, xyI{n,5} = real intensity
% raw_parameter = [Ip, a, c, 2b, Ib, x0, y0, sigma1, sigma2, total fluorescence intensity, spot#]

%figure;surf(im3,'EdgeColor','none');colorbar;title('after watershed');view(0, 90);xlim([0 512]);ylim([0 512]);title('4. after imopen');



% ===== Step 8. Final overlay figures =====
overlay = imoverlay(raw_im./max(max(raw_im)),bwperim(im_targets),[0,0,1]); % nucleus
for se=1:size(new_perim,2)
    overlay = imoverlay(overlay,new_perim{se},[0,1/size(new_perim,2)*se,0]); % nucleus
end

%figure;
axes(handles.axes2);
imshow(overlay);hold on;
for j=1:size(newBoxPolygon,1)
    line(newBoxPolygon{j,1}(:, 1), newBoxPolygon{j,1}(:, 2), 'Color', 'y');hold on
end
hold off
for i=1:size(raw_parameter,1)
    text(raw_parameter(i,7),raw_parameter(i,6),num2str(round(raw_parameter(i,11))),'Color','white','FontSize',11);
end



end
function [raw_fishall,raw_fish]=image_background(tiff_im,tiff_num,max_image,back_mask,p)

% Estimate background intensity for all stacks
back_all = zeros(size(tiff_im,2),1);
for ii=1:size(tiff_im,2)
    imall(:,:,ii)=tiff_im{1,ii}(:,:,p.FISH_chal);
    back_all(ii,1) = sum(double(imall(:,:,ii)).*double(back_mask))/sum(double(back_mask));
    raw_fishall(:,:,ii)=max(imall(:,:,ii) - back_all(ii,1),0);
end

if tiff_num==0
    im=max_image(:,:,p.FISH_chal);
    back = sum(double(im).*double(back_mask))/sum(double(back_mask));
    raw_fish = max(im-back,0);
else
    raw_fish = max(raw_fishall(:,:,tiff_num),0);
end


end
function [new_nucleus_prop,new_perim]=DivideTargetRegion(raw_im,nucleus_prop,nucleus_center,newBoxPolygon)

IN = inpolygon(nucleus_center(:, 1),nucleus_center(:, 2),newBoxPolygon(:, 1),newBoxPolygon(:, 2));
IN=find(IN);
ii=1;
for i=1:size(nucleus_prop,1)
    if isempty(find(i-IN==0, 1))==1
        new_nucleus_prop(ii,1)=nucleus_prop(i);
        ii=ii+1;
    end
end

TargetRegion = false(size(raw_im));
for i=1:size(IN,1)
    TargetRegion(nucleus_center(IN(i,1),4):nucleus_center(IN(i,1),6),nucleus_center(IN(i,1),3):nucleus_center(IN(i,1),5))=nucleus_prop(IN(i,1),1).Image;
end
[Targetx,Targety]=find(TargetRegion);


% cell1
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon(2:4, 1),newBoxPolygon(2:4, 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{1}=cell1;

% cell2
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon([2,1,12], 1),newBoxPolygon([2,1,12], 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{2}=cell1;

% cell3
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon(10:12, 1),newBoxPolygon(10:12, 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{3}=cell1;

% cell4
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon(8:10, 1),newBoxPolygon(8:10, 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{4}=cell1;

% cell5
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon(7:8, 1),newBoxPolygon(7:8, 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{5}=cell1;

% cell6
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon(6:7, 1),newBoxPolygon(6:7, 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{6}=cell1;

% cell7
poly1=[[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))];newBoxPolygon(4:6, 1),newBoxPolygon(4:6, 2);[mean(newBoxPolygon(:,1)),mean(newBoxPolygon(:,2))]];
IN1 = inpolygon(Targety,Targetx,poly1(:, 1),poly1(:, 2));
cell1=false(size(raw_im));IN1=find(IN1);
for ii=1:size(IN1,1)
    cell1(Targetx(IN1(ii)),Targety(IN1(ii)))=1;
end
new_nucleus_prop=[new_nucleus_prop;regionprops(cell1,'Area','Centroid','Image','BoundingBox')];
new_perim{7}=cell1;

end
function nucle_blob_overlay = make_overlaymap(STATS,raw)

for i=1:size(STATS,1)
    nucle_blob_overlay{i}=false(size(raw));
    [xx,yy]=find(STATS(i,1).Image==1);
    x0=STATS(i,1).BoundingBox(1)-0.5;y0=STATS(i,1).BoundingBox(2)-0.5;
    xxx=yy+x0;yyy=xx+y0;
    for ii=1:size(xxx,1)
        nucle_blob_overlay{i}(yyy(ii),xxx(ii))=1;
    end
end
end
function STATS2=STAT_subBoundary(STATS,max_image)

n=size(STATS,1);
ii=1;
for i=1:n
    A=STATS(i,1).BoundingBox;
    A(3)=A(3)+A(1);A(4)=A(4)+A(2);
    if A(1)<1||A(2)<1||A(3)>size(max_image,1)||A(4)>size(max_image,2)
    else
        STATS2(ii,1)=STATS(i,1);
        ii=ii+1;
    end
end
if ii==1
   STATS2 = []; 
end


end
function back_mask = Intemask(several_mask)

if size(several_mask,2)>=1
    back_mask=zeros(size(several_mask{1,1}));
    for i=1:size(several_mask,2)
        back_mask = back_mask|several_mask{1,i};
    end
else
    back_mask=[];
end

end
function [overlay1,bn_matrix_est,back_blob_overlay,nucle_blob,back_mask] = nuc_back_integrate(p,STATS_back,STATS_nuc,raw_fish,back_mask,nucle_blob_overlay,back_blob_overlay)

% back_area matrix
back_area=zeros(size(STATS_back,1),5);
for i=1:size(STATS_back,1)
    back_area(i,1) = i;
    back_area(i,2) = STATS_back(i,1).Centroid(1);
    back_area(i,3) = STATS_back(i,1).Centroid(2);
    back_area(i,4) = STATS_back(i,1).Area;
    back_area(i,5) = STATS_back(i,1).Solidity;
end

% nucle_lob matrix
nucle_blob=zeros(size(STATS_nuc,1),7);  %nucle_blob=[no., x, y, A, solidity, Backno., distance]
for i=1:size(STATS_nuc,1)
    nucle_blob(i,1) = i;
    nucle_blob(i,2) = STATS_nuc(i,1).Centroid(1);
    nucle_blob(i,3) = STATS_nuc(i,1).Centroid(2);
    nucle_blob(i,4) = STATS_nuc(i,1).Area;
    nucle_blob(i,5) = STATS_nuc(i,1).Solidity;
end

% Sort back_area and find corresponding nucle
back_area = sortrows(back_area,-4);back_drange=zeros(size(back_area,1),2);
for ii=1:size(back_area,1)
    bn_dist = zeros(size(nucle_blob,1),1);
    for jj=1:size(nucle_blob,1)
        bn_dist(jj,1) = sqrt((nucle_blob(jj,2)-back_area(ii,2))^2+(nucle_blob(jj,3)-back_area(ii,3))^2);
    end
    % find nearest nuclear blobs
    bn_matrix = [nucle_blob(:,1:5) bn_dist]; %bn_matrix=[no., x, y, A, solidity, distance]
    bn_matrix = sortrows(bn_matrix,6);
    d=diff(bn_matrix(:,6));
    
    % nearby nucleus identification method 1: diff peak find method 
    LocMin1=[1:100]';thresh=0;
    while size(LocMin1,1)>1&&thresh<500
        [LocMin1] = peakfinder(d(p.candi_start:p.candi_end),thresh,5,1);
        % peakfinder(inputs, threshold, smaller than to be minima, -1 if minima)
        % threshold: Larger values mean the algorithm is more selective in finding peaks
        thresh=thresh+1;LocMin1=LocMin1+2;
    end
    % nearby nucleus identification method 2: threshold radius method
    LocMin2 = max(find(bn_matrix(:,6)<p.nearby_radius_max));
    % combine method 1,2
    bn_matrix_est{ii} = bn_matrix(1:min(LocMin1,LocMin2),:);   
    back_drange(ii,1) = bn_matrix_est{ii}(1,6);
    back_drange(ii,2) = bn_matrix_est{ii}(end,6);
    temp=[];   
    for jj=1:size(bn_matrix_est{ii},1)
        dis=bn_matrix_est{ii}(jj,6);
        if nucle_blob(bn_matrix_est{ii}(jj,1),7)==0 %||d<nucle_blob(bn_matrix_est(jj,1),7)
            nucle_blob(bn_matrix_est{ii}(jj,1),6)=ii;
            nucle_blob(bn_matrix_est{ii}(jj,1),7)=dis;
        else
            temp=[temp jj];
        end
    end
    otherlist{ii}=temp;
end
    
    % calculate overlay figures
    iiii=1;nucle_blob2=nucle_blob;
    for iii=1:size(back_area,1)
        if isempty(find(nucle_blob(:,6)==iii, 1))==1||(isempty(find(otherlist{1,iii}==1))~=1&&isempty(find(otherlist{1,iii}==2))~=1&&isempty(find(otherlist{1,iii}==3))~=1)
            nucle_blob2(find(nucle_blob2(:,6)==iii),6)=0;
        else
            back_blob_overlay2{iiii}=back_blob_overlay{back_area(iii,1)};
            bn_matrix_est2{iiii}=bn_matrix_est{iii};
            nucle_blob2(find(nucle_blob2(:,6)==iii),6)=iiii;
            iiii=iiii+1;
        end       
    end
    back_mask2 = Intemask(back_blob_overlay2);
    
    nucle_blob3 = sortrows(nucle_blob2,6);
    %nucle_blob2=nucle_blob2(nucle_blob2(:,6)>0,:);
    ma=max(nucle_blob3(:,6));
    ma_color=rand(ma,3);
    overlay1 = imoverlay(raw_fish./max(max(raw_fish)),bwperim(back_mask2),[0,1,0]); % background
    for oo=1:size(nucle_blob3,1)
        if nucle_blob3(oo,6)==0
            overlay1 = imoverlay(overlay1,bwperim(nucle_blob_overlay{1,nucle_blob3(oo,1)}),[0,0,1]); % nucleus
        else
            overlay1 = imoverlay(overlay1,nucle_blob_overlay{1,nucle_blob3(oo,1)},ma_color(nucle_blob3(oo,6),:)); % nucleus
        end
    end
    
    % output
    nucle_blob = nucle_blob3;
    back_mask = back_mask2;
    bn_matrix_est = bn_matrix_est2;
    back_blob_overlay = back_blob_overlay2;

end

%% Gaussian fitting 
function [raw_parameter,xyI,ex_cell] = ptrack_par_nogauss(p,max_image0noback,nucleus_prop)
% function peak_parameter = ptrack_par(imstack0,imstack,mask_out)
% For each integrated peak(mRNA spot)
% peak_parameter[N,11] = [Ip, sigma1, sigma2, ?, Ib, x0, y0, ?, ?, total fluorescence intensity]
% xyI{N} = [I(predicted), I(data), x,y,z]  for every locations in each peark
% raw_cell{Ifit} = [Ip, a, x0, c, y0, 2b, xf1(:,2:3), Ib, resum(?)];
tic

% Gaussian filter smoothing
H = fspecial(p.gau_type,p.gau_hsize,p.gau_sigma);
%H = -fspecial('log',p.gfit_hsize,p.gfit_sigma);
%H = fspecial('gaussian',p.RnaSe_hsize);
max_image0f = imfilter(max_image0noback,H,'replicate');
max_image0f(max_image0f<0) = 0;
mask_out=zeros(size(max_image0f));
for ii=1:size(nucleus_prop,1)
    mask_out(round(nucleus_prop(ii).Centroid(2)),round(nucleus_prop(ii).Centroid(1)))=1;
end

imstack0 = max_image0noback;
imstack = max_image0f;

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xrange = 3;  %%% Peak fitting range on x direction 6 -->3
yrange = 3;  %%% Peak fitting range on y direction 6 -->3
zrange = 0;  %%% Peak fitting range on z direction
sigmax0 = 2;   %%% Initial value of sigma_x
sigmay0 = 2;   %%% Initial value of sigma_y
sigmaz0 = 1;   %%% Initial value of sigma_z
crossxy = 0;   %%% Initial value of cross coefficient for "xy" terms
gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator
%gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))+x(',num2str(n),',4)'];   %%% 1D single-gaussian model function text generator
%gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-xdata(:,4)).^2+x(',num2str(n),',3)*(xdata(:,2)-xdata(:,5)).^2)+x(',num2str(n),',4)*(xdata(:,1)-xdata(:,4)).*(xdata(:,2)-xdata(:,5)))+x(',num2str(n),',5)'];   %%% 2D single-gaussian model function text generator
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Coordinate matrix generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim0 = size(imstack0);
pdim0 = prod(dim0);
%[XC,YC,ZC] = ndgrid(1:dim0(1),1:dim0(2),1:dim0(3));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak sorting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mask_out(:,:,[1,end]) = false;
%Ipeak = find(double(imstack).*mask_out);
Ipeak = find(double(imstack).*mask_out);
[~,Itr] = sort(imstack(Ipeak),'descend');
Ipeak = Ipeak(Itr);
[X0,Y0,Z0] = ind2sub(dim0,Ipeak);   %%% Peak coordinates list

Nucleus0=zeros(size(X0,1),5);Cen_xy=zeros(size(X0,2));
for iii=1:size(X0,1)
    Cen_xy(iii,1:2)=nucleus_prop(iii,1).Centroid;
end
for iii=1:size(X0,1)
    [~,aa]=min((Cen_xy(:,1)-Y0(iii)).^2+(Cen_xy(:,2)-X0(iii)).^2);
    Nucleus0(iii,1)=aa;
    Nucleus0(iii,2)=nucleus_prop(aa).BoundingBox(2)-0.5;Nucleus0(iii,3)=nucleus_prop(aa).BoundingBox(2)+nucleus_prop(aa).BoundingBox(4)-0.5;
    Nucleus0(iii,4)=nucleus_prop(aa).BoundingBox(1)-0.5;Nucleus0(iii,5)=nucleus_prop(aa).BoundingBox(1)+nucleus_prop(aa).BoundingBox(3)-0.5;
end

% X0 = XC(Ipeak);   %%% Peak coordinates x list
% Y0 = YC(Ipeak);   %%% Peak coordinates y list
% Z0 = ZC(Ipeak);   %%% Peak coordinates z list
vpeak = imstack0(Ipeak);   %%% Peak value list
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = [];%zeros(length(Ipeak),10);
%exstatus = false(length(Ipeak),1);
all_offset = zeros(pdim0,1);

lb = [-inf,-inf,-inf,-inf,-inf,-inf,-inf];
ub = [inf,inf,inf,inf,inf,inf,inf];

options = optimset('Display','off');

%%% Arrange fitting groups: %%% ===========================================
% Dx = pdist(X0);   %%%%% distance matrix on x dimension
% Dy = pdist(Y0);   %%%%% distance matrix on y dimension
% Dz = pdist(Z0);   %%%%% distance matrix on z dimension
% Dxyz = squareform((pdist(X0) <= 2*xrange)&(pdist(Y0) <= 2*yrange)&(pdist(Z0) <= 2*zrange)); %%%%% neighbor matrix
isfit = false(size(X0));   %%% Peak fitting status
If_peak = cell(0);
fmax = 1;

for Ifit = 1:length(Ipeak)
    if ~isfit(Ifit)   %%% check whether the peak has been fitted
        Jfit = Ifit;   %%% instant peak in process
        Nfit = [];   %%% total peaks in fitting
        %%%%% Fitting range/peak arrangement: %%%%% -----------------------
        while Jfit
            Nfit = union(Nfit,Jfit);
            Dxyz = (pdist2(X0(Jfit),X0) <= 2*xrange)&(pdist2(Y0(Jfit),Y0) <= 2*yrange)&(pdist2(Z0(Jfit),Z0) <= 2*zrange);
            [~,Jfit] = find(Dxyz);
            Jfit = setdiff(Jfit,Nfit);
        end
        %%%%%% ------------------------------------------------------------
        isfit(Nfit) = true;
        If_peak = cat(1,If_peak,{Nfit});
        fmax = max(fmax,length(Nfit));
    end
end
clear Dxyz
%%% =======================================================================

        
%%% Fitting preparation: %%% ==============================================
%%%%% 2D multi-gaussian fit function generation:
gau2 = cell(1,fmax);
for I_gau = 1:fmax
    textfun = '@(x,xdata) ';
    for n = 1:I_gau
        textfun = [textfun,gau2_gen(n),'+'];
    end
    textfun = [textfun(1:(end-1)),';'];
    gau2{I_gau} = eval(textfun);
end
        
%%%%% Data points collection: 
%raw_cell = cell(length(If_peak),1);
ex_cell = cell(length(If_peak),1);

dim0(3)=1;
Ifit_xyI=1;
for Ifit = 1:length(If_peak)
    Nfit = If_peak{Ifit};
    Nfit = reshape(Nfit,numel(Nfit),1);
    xdata = zeros(0,5);
    xy0 = [];xminmax=[10000 0];yminmax=[10000 0];
    for I_temp = 1:length(Nfit)
        %[Xtemp,Ytemp,Ztemp] = ndgrid(max((X0(Nfit(I_temp))-xrange),1):min((X0(Nfit(I_temp))+xrange),dim0(1)),max((Y0(Nfit(I_temp))-yrange),1):min((Y0(Nfit(I_temp))+yrange),dim0(2)),max((Z0(Nfit(I_temp))-zrange),1):min((Z0(Nfit(I_temp))+zrange),dim0(3)));
        [Xtemp,Ytemp,Ztemp] = ndgrid(max(Nucleus0(Nfit(I_temp),2),1):min(Nucleus0(Nfit(I_temp),3),dim0(1)),max(Nucleus0(Nfit(I_temp),4),1):min(Nucleus0(Nfit(I_temp),5),dim0(2)),max((Z0(Nfit(I_temp))-zrange),1):min((Z0(Nfit(I_temp))+zrange),dim0(3)));
        xdata = union(xdata,[Xtemp(:),Ytemp(:),Ztemp(:),X0(Nfit(I_temp))*ones(size(Xtemp(:))),Y0(Nfit(I_temp))*ones(size(Xtemp(:)))],'rows');
        xy0=[xy0;[X0(Nfit(I_temp)) Y0(Nfit(I_temp))]];
        xminmax=[min(min(min(Xtemp)),xminmax(1)) max(max(max(Xtemp)),xminmax(2))];
        yminmax=[min(min(min(Ytemp)),yminmax(1)) max(max(max(Ytemp)),yminmax(2))];
    end
    ydata =imstack0(sub2ind(dim0,xdata(:,1),xdata(:,2),xdata(:,3)));
%%%%% Fitting parameter initialization:
    %xf0 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), X0(Nfit), 1./2./sigmay0.^2.*ones(size(Nfit)), Y0(Nfit), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit))];
    %xf0 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), 1./2./sigmay0.^2.*ones(size(Nfit)), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit))];
    %xf10 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit)), Z0(Nfit), min(ydata)*ones(size(Nfit))];

        raw_cell = [];
        raw_cell(:,1:5) = zeros(1,5);
        raw_cell(:,6:7) = Cen_xy(Ifit,[2,1]);
        raw_cell(:,8) = 1./sqrt(raw_cell(:,2)+raw_cell(:,3)+sqrt((raw_cell(:,2)-raw_cell(:,3)).^2+raw_cell(:,4).^2));
        raw_cell(:,9) = 1./sqrt(raw_cell(:,2)+raw_cell(:,3)-sqrt((raw_cell(:,2)-raw_cell(:,3)).^2+raw_cell(:,4).^2));
        raw_cell(:,10) = sum(ydata); %2*pi*raw_cell(:,1).*raw_cell(:,8).*raw_cell(:,9);
        raw_cell(:,11) = sum(ydata); %Ifit_xyI;
        xyI{Ifit_xyI,1} = [ydata ydata xdata];
        raw_parameter=[raw_parameter;real(raw_cell)];
        
        [xyI{Ifit_xyI,2},xyI{Ifit_xyI,3}] = ndgrid(xminmax(1):xminmax(2), yminmax(1):yminmax(2));
        xyI{Ifit_xyI,4}=zeros(xminmax(2)-xminmax(1)+1,yminmax(2)-yminmax(1)+1);
        xyI{Ifit_xyI,5}=zeros(xminmax(2)-xminmax(1)+1,yminmax(2)-yminmax(1)+1);
        xyI{Ifit_xyI,6}=0; %resnorm;
        for ii=1:size(xyI{Ifit_xyI,1},1)
            xyI{Ifit_xyI,4}(xyI{Ifit_xyI,1}(ii,3)-xminmax(1)+1,xyI{Ifit_xyI,1}(ii,4)-yminmax(1)+1)=xyI{Ifit_xyI,1}(ii,1);
            xyI{Ifit_xyI,5}(xyI{Ifit_xyI,1}(ii,3)-xminmax(1)+1,xyI{Ifit_xyI,1}(ii,4)-yminmax(1)+1)=xyI{Ifit_xyI,1}(ii,2);
        end
        Ifit_xyI=Ifit_xyI+1;
    %end
%%% =======================================================================
end
raw_parameter=distn(raw_parameter,nucleus_prop);

toc
end
function raw_parameter_new=distn(raw_parameter,nucleus_prop)
% raw_parameter: [nx11]double = [Ip, a, c, 2b, Ib, x0, y0, sigma1, sigma2, total fluorescence intensity, spot#]
% nucleus_prop: [mx1]struct

m=size(nucleus_prop,1);B=zeros(m,2);
n=size(raw_parameter,1);A=zeros(n,1);
for i=1:m
    B(i,:)=nucleus_prop(i).Centroid;
end

for j=1:n
    [~,A(j,1)]=min((raw_parameter(j,7)-B(:,1)).^2+(raw_parameter(j,6)-B(:,2)).^2);
end
raw_parameter_new=[raw_parameter A];
end
function out = imoverlay(in, mask, color)
%IMOVERLAY Create a mask-based image overlay.
%   OUT = IMOVERLAY(IN, MASK, COLOR) takes an input image, IN, and a binary
%   image, MASK, and produces an output image whose pixels in the MASK
%   locations have the specified COLOR.
%
%   IN should be a grayscale or an RGB image of class uint8, uint16, int16,
%   logical, double, or single.  If IN is double or single, it should be in
%   the range [0, 1].  If it is not in that range, you might want to use
%   mat2gray to scale it into that range.
%
%   MASK should be a two-dimensional logical matrix.
%
%   COLOR should be a 1-by-3 vector of values in the range [0, 1].  [0 0 0]
%   is black, and [1 1 1] is white.
%
%   OUT is a uint8 RGB image.
%
%   Examples
%   --------
%   Overlay edge detection result in green over the original image.
%       
%       I = imread('cameraman.tif');
%       bw = edge(I, 'canny');
%       rgb = imoverlay(I, bw, [0 1 0]);
%       imshow(rgb)
%
%   Treating the output of peaks as an image, overlay the values greater than
%   7 in red.  The output of peaks is not in the usual grayscale image range
%   of [0, 1], so use mat2gray to scale it.
%
%       I = peaks;
%       mask = I > 7;
%       rgb = imoverlay(mat2gray(I), mask, [1 0 0]);
%       imshow(rgb, 'InitialMagnification', 'fit')

%   Steven L. Eddins, The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2007/08/15 13:18:08 $

% If the user doesn't specify the color, use white.
DEFAULT_COLOR = [1 1 1];
if nargin < 3
    color = DEFAULT_COLOR;
end

% Make the uint8 the working data class.  The output is also uint8.
in_uint8 = im2uint8(in);
color_uint8 = im2uint8(color);

% Initialize the red, green, and blue output channels.
if ndims(in_uint8) == 2
    % Input is grayscale.  Initialize all output channels the same.
    out_red   = in_uint8;
    out_green = in_uint8;
    out_blue  = in_uint8;
else
    % Input is RGB truecolor.
    out_red   = in_uint8(:,:,1);
    out_green = in_uint8(:,:,2);
    out_blue  = in_uint8(:,:,3);
end

% Replace output channel values in the mask locations with the appropriate
% color value.
out_red(mask)   = color_uint8(1);
out_green(mask) = color_uint8(2);
out_blue(mask)  = color_uint8(3);

% Form an RGB truecolor image by concatenating the channel matrices along
% the third dimension.
out = cat(3, out_red, out_green, out_blue);
end
function max_image = maximage(tiff_images,p)

% size of original tiff cells
zn=size(tiff_images,2);
xysize=size(tiff_images{1,1});

% choice init and final stake considered
init_stack=p.init_stack;
final_stack=p.final_stack;

%max_image=zeros(xysize(1),xysize(2),xysize(3));
%max_image=max(tiff_images{1,:}(1,1,1));
max_image1=zeros(xysize(1),xysize(2),zn);
max_image2=zeros(xysize(1),xysize(2),zn);
max_image3=zeros(xysize(1),xysize(2),zn);
for i=init_stack:final_stack
    max_image1(:,:,i-init_stack+1)=tiff_images{1,i}(:,:,1);
    max_image2(:,:,i-init_stack+1)=tiff_images{1,i}(:,:,2);
    max_image3(:,:,i-init_stack+1)=tiff_images{1,i}(:,:,3);
end
max_image(:,:,1)=max(max_image1,[],3);
max_image(:,:,2)=max(max_image2,[],3);
max_image(:,:,3)=max(max_image3,[],3);

end
function [propA,back_mask,nucle_mask,overlay1,back_blob_overlay]=conu_theh(im,max_image,p,handles)

% parameter setting
th_num=p.back_th_num;
max_value=max(max(max_image));
min_value=min(min(max_image));
d=floor((max_value-min_value)/th_num);
Area=size(max_image,1)*size(max_image,2);
propA=zeros(th_num,7);   
% propA = [threshold, nuclear#, nuclear area std, nuclear weighted solidity,...
%          back#, back area std, back weighted solidity]

% try different threshold values
for i=1:th_num
    th=min_value+d*(i-1);
    temp_nuc = max_image > th;
    temp_back = max_image < th;
    temp_nuc1 = bwareaopen(temp_nuc,50);
    temp_back1 = bwareaopen(temp_back,500);
    STATS_nuc = regionprops(temp_nuc1,'Centroid','Area','Solidity','BoundingBox');
    STATS_back = regionprops(temp_back1,'Centroid','Area','Solidity','BoundingBox');
    STATS_back=STAT_subBoundary(STATS_back,max_image);
    
    temp_nuc=zeros(size(STATS_nuc,1),2);
    for ii=1:size(STATS_nuc,1)
        temp_nuc(ii,1)=STATS_nuc(ii,1).Area;
    end
    for ii=1:size(STATS_nuc,1)
        temp_nuc(ii,2)=temp_nuc(ii,1)/sum(temp_nuc(:,1))*STATS_nuc(ii,1).Solidity;
    end
      
    temp_back=zeros(size(STATS_back,1),2);
    for ii=1:size(STATS_back,1)
        temp_back(ii,1)=STATS_back(ii,1).Area;
    end
    for ii=1:size(STATS_back,1)
        temp_back(ii,2)=temp_back(ii,1)/sum(temp_back(:,1))*STATS_back(ii,1).Solidity;
    end
    
    propA(i,1)=th;             % threshold value
    propA(i,2)=size(STATS_nuc,1);  % nuclear  number 
    propA(i,3)=std(temp_nuc(:,1));  % nuclear std(area)/mean(area)
    propA(i,4)=sum(temp_nuc(:,2));        % nuclear weighted averaged Solidity
    propA(i,5)=size(STATS_back,1);  %background  number
    propA(i,6)=std(temp_back(:,1));  % background std(area)/mean(area)
    propA(i,7)=sum(temp_back(:,2));       % background weighted averaged Solidity
end


    % Create background mask
    % method 1
    %{
    peakLoc=(1:100)';thresh=0;
    while size(peakLoc(peakLoc(:,1)>10&peakLoc(:,1)<90),1)>2&&thresh<100
        [peakLoc] = peakfinder(propA(:,2),thresh,max(propA(:,2))/5,-1);
        thresh=thresh+1;
    end
    back_mask = max_image < propA(peakLoc(end,1),1);n=peakLoc(end,1);
    %}
    % method 2
    [~,n] = min(propA(:,4));
    back_mask = max_image < propA(n,1);
    
    se = strel('disk', 10);
    back_mask = imerode(back_mask, se);
    back_mask= imfill(back_mask,'holes');
    back_mask = bwareaopen(back_mask,500);
    
    % delete boundary background
    STATS_back1 = regionprops(back_mask,'Centroid','Area','BoundingBox','Image');
    if p.remove_bou_back==1
        STATS_back2=STAT_subBoundary(STATS_back1,max_image);   %remove boundary backgroind blobs
    else
        STATS_back2=STATS_back1;
    end
    back_blob_overlay = make_overlaymap(STATS_back2,max_image);   %create back blob maps
    back_mask = Intemask(back_blob_overlay);   %Integrate several blob maps into one


    % Create nuclear mask
    %method 1
    %{
    try_th = 0.83;
    nu_temp = find(propA(:,4)>try_th);
    while isempty(nu_temp)
        nu_temp = find(propA(:,4)>try_th);
        try_th = try_th -0.05;
    end
    nuc_threshold = nu_temp(min(find(nu_temp>50)));
    nucle_mask = max_image > propA(nuc_threshold,1);
    %}
    %method 2
    LocMin=[1:100]';thresh=0;
    while size(LocMin,1)>3&&thresh<500
        [LocMin] = peakfinder(propA(51:100,2),thresh,200,-1);
        LocMin = LocMin+50;
        % peakfinder(inputs, threshold, smaller than to be minima, -1 if minima)
        % threshold: Larger values mean the algorithm is more selective in finding peaks
        thresh=thresh+1;
    end    
    nucle_mask = max_image > propA(LocMin(2),1);
    
    se = strel('disk', 5);
    nucle_mask = imopen(nucle_mask, se);
    nucle_mask = imfill(nucle_mask,'hole');
    nucle_mask = bwareaopen(nucle_mask,50);
    
    
    
    % calculate overlay figures
    overlay1 = imoverlay(im./max(max(im)),bwperim(back_mask),[0,1,0]); % background
    overlay1 = imoverlay(overlay1,nucle_mask,[0,0,1]); % nucleus
    
    %figure;
    axes(handles.axes2);
    imshow(overlay1);pause(0.1);
   
if p.closeallfigure~=1
    %figure;AX=plotyy(propA(:,1),propA(:,2),propA(:,1),propA(:,4));hold on
    figure;plot(propA(:,1),propA(:,4));hold on
    %axis(AX(1)); axis(AX(2)); 
    xlabel('threshold value');
    ylabel('nuclear weighted averaged Solidity');
    %ylabel(AX(1),'background  spot number'); ylabel(AX(2),'nuclear weighted averaged Solidity');
    plot([propA(n,1),propA(n,1)],[min(propA(:,4)),max(propA(:,4))],'b-.');hold on;
    %plot([propA(peakLoc(end,1),1),propA(peakLoc(end,1),1)],[0,20000],'b-.');hold on;
    %plot([0.8*max_value,0.8*max_value],[0,20000],'b-.');hold off;
    %plot([propA(n,1),propA(n,1)],[0,20000],'g-.');hold on;
    plot([propA(nuc_threshold,1),propA(nuc_threshold,1)],[min(propA(:,4)),max(propA(:,4))],'r-.');hold off;
    
    %figure;imshow(overlay1);
end
end
function varargout = peakfinder(x0, sel, thresh, extrema)
        %PEAKFINDER Noise tolerant fast peak finding algorithm
        %   INPUTS:
        %       x0 - A real vector from the maxima will be found (required)
        %       sel - The amount above surrounding data for a peak to be
        %           identified (default = (max(x0)-min(x0))/4). Larger values mean
        %           the algorithm is more selective in finding peaks.
        %       thresh - A threshold value which peaks must be larger than to be
        %           maxima or smaller than to be minima.
        %       extrema - 1 if maxima are desired, -1 if minima are desired
        %           (default = maxima, 1)
        %   OUTPUTS:
        %       peakLoc - The indicies of the identified peaks in x0
        %       peakMag - The magnitude of the identified peaks
        %
        %   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
        %       are at least 1/4 the range of the data above surrounding data.
        %
        %   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
        %       that are at least sel above surrounding data.
        %
        %   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
        %       maxima that are at least sel above surrounding data and larger
        %       (smaller) than thresh if you are finding maxima (minima).
        %
        %   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
        %       data if extrema > 0 and the minima of the data if extrema < 0
        %
        %   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
        %       local maxima as well as the magnitudes of those maxima
        %
        %   If called with no output the identified maxima will be plotted along
        %       with the input data.
        %
        %   Note: If repeated values are found the first is identified as the peak
        %
        % Ex:
        % t = 0:.0001:10;
        % x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
        % x(1250:1255) = max(x);
        % peakfinder(x)
        %
        % Copyright Nathanael C. Yoder 2011 (nyoder@gmail.com)

        % Perform error checking and set defaults if not passed in
        error(nargchk(1,4,nargin,'struct'));
        error(nargoutchk(0,2,nargout,'struct'));

        s = size(x0);
        flipData =  s(1) < s(2);
        len0 = numel(x0);
        if len0 ~= s(1) && len0 ~= s(2)
            error('PEAKFINDER:Input','The input data must be a vector')
        elseif isempty(x0)
            varargout = {[],[]};
            return;
        end
        if ~isreal(x0)
            warning('PEAKFINDER:NotReal','Absolute value of data will be used')
            x0 = abs(x0);
        end

        if nargin < 2 || isempty(sel)
            sel = (max(x0)-min(x0))/4;
        elseif ~isnumeric(sel) || ~isreal(sel)
            sel = (max(x0)-min(x0))/4;
            warning('PEAKFINDER:InvalidSel',...
                'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
        elseif numel(sel) > 1
            warning('PEAKFINDER:InvalidSel',...
                'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
            sel = sel(1);
        end

        if nargin < 3 || isempty(thresh)
            thresh = [];
        elseif ~isnumeric(thresh) || ~isreal(thresh)
            thresh = [];
            warning('PEAKFINDER:InvalidThreshold',...
                'The threshold must be a real scalar. No threshold will be used.')
        elseif numel(thresh) > 1
            thresh = thresh(1);
            warning('PEAKFINDER:InvalidThreshold',...
                'The threshold must be a scalar.  The first threshold value in the vector will be used.')
        end

        if nargin < 4 || isempty(extrema)
            extrema = 1;
        else
            extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
            if extrema == 0
                error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
            end
        end

        x0 = extrema*x0(:); % Make it so we are finding maxima regardless
        thresh = thresh*extrema; % Adjust threshold according to extrema.
        dx0 = diff(x0); % Find derivative
        dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
        ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign

        % Include endpoints in potential peaks and valleys
        x = [x0(1);x0(ind);x0(end)];
        ind = [1;ind;len0];

        % x only has the peaks, valleys, and endpoints
        len = numel(x);
        minMag = min(x);


        if len > 2 % Function with peaks and valleys

            % Set initial parameters for loop
            tempMag = minMag;
            foundPeak = false;
            leftMin = minMag;

            % Deal with first point a little differently since tacked it on
            % Calculate the sign of the derivative since we taked the first point
            %  on it does not neccessarily alternate like the rest.
            signDx = sign(diff(x(1:3)));
            if signDx(1) <= 0 % The first point is larger or equal to the second
                ii = 0;
                if signDx(1) == signDx(2) % Want alternating signs
                    x(2) = [];
                    ind(2) = [];
                    len = len-1;
                end
            else % First point is smaller than the second
                ii = 1;
                if signDx(1) == signDx(2) % Want alternating signs
                    x(1) = [];
                    ind(1) = [];
                    len = len-1;
                end
            end

            % Preallocate max number of maxima
            maxPeaks = ceil(len/2);
            peakLoc = zeros(maxPeaks,1);
            peakMag = zeros(maxPeaks,1);
            cInd = 1;
            % Loop through extrema which should be peaks and then valleys
            while ii < len
                ii = ii+1; % This is a peak
                % Reset peak finding if we had a peak and the next peak is bigger
                %   than the last or the left min was small enough to reset.
                if foundPeak
                    tempMag = minMag;
                    foundPeak = false;
                end

                % Make sure we don't iterate past the length of our vector
                if ii == len
                    break; % We assign the last point differently out of the loop
                end

                % Found new peak that was lager than temp mag and selectivity larger
                %   than the minimum to its left.
                if x(ii) > tempMag && x(ii) > leftMin + sel
                    tempLoc = ii;
                    tempMag = x(ii);
                end

                ii = ii+1; % Move onto the valley
                % Come down at least sel from peak
                if ~foundPeak && tempMag > sel + x(ii)
                    foundPeak = true; % We have found a peak
                    leftMin = x(ii);
                    peakLoc(cInd) = tempLoc; % Add peak to index
                    peakMag(cInd) = tempMag;
                    cInd = cInd+1;
                elseif x(ii) < leftMin % New left minima
                    leftMin = x(ii);
                end
            end

            % Check end point
            if x(end) > tempMag && x(end) > leftMin + sel
                peakLoc(cInd) = len;
                peakMag(cInd) = x(end);
                cInd = cInd + 1;
            elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
                peakLoc(cInd) = tempLoc;
                peakMag(cInd) = tempMag;
                cInd = cInd + 1;
            end

            % Create output
            peakInds = ind(peakLoc(1:cInd-1));
            peakMags = peakMag(1:cInd-1);
        else % This is a monotone function where an endpoint is the only peak
            [peakMags,xInd] = max(x);
            if peakMags > minMag + sel
                peakInds = ind(xInd);
            else
                peakMags = [];
                peakInds = [];
            end
        end

        % Apply threshold value.  Since always finding maxima it will always be
        %   larger than the thresh.
        if ~isempty(thresh)
            m = peakMags>thresh;
            peakInds = peakInds(m);
            peakMags = peakMags(m);
        end



        % Rotate data if needed
        if flipData
            peakMags = peakMags.';
            peakInds = peakInds.';
        end



        % Change sign of data if was finding minima
        if extrema < 0
            peakMags = -peakMags;
            x0 = -x0;
        end
        % Plot if no output desired
        if nargout == 0
            if isempty(peakInds)
                disp('No significant peaks found')
            else
                figure;
                plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
            end
        else
            varargout = {peakInds,peakMags};
        end
    end

%% Object detection algorithm
function newBoxPolygon=ObjectDetection(p,sceneImage,boxImage)


newBoxPolygon{1}=ObjectDetection_SURF(p,sceneImage,boxImage);
%ObjectDetection_MSER(sceneImage,boxImage);
%ObjectDetection_Cascade(sceneImage,boxImage);
end
function newBoxPolygon=ObjectDetection_SURF(p,sceneImage,boxImage)

% Step 1: Load images
%load matlab.mat;


% Step 2: Detect Feature Points
boxPoints = detectSURFFeatures(boxImage);
scenePoints = detectSURFFeatures(sceneImage);

if p.closeallfigure~=1
    figure;
    imshow(boxImage);
    title('100 Strongest Feature Points from Box Image');
    hold on;
    plot(selectStrongest(boxPoints, 100));

    figure;
    imshow(sceneImage);
    title('300 Strongest Feature Points from Scene Image');
    hold on;
    plot(selectStrongest(scenePoints, 300));
end

% Step 3: Extract Feature Descriptors
[boxFeatures, boxPoints] = extractFeatures(boxImage, boxPoints);
[sceneFeatures, scenePoints] = extractFeatures(sceneImage, scenePoints);


% Step 4: Find Putative Point Matches
boxPairs = matchFeatures(boxFeatures, sceneFeatures);
matchedBoxPoints = boxPoints(boxPairs(:, 1), :);
matchedScenePoints = scenePoints(boxPairs(:, 2), :);
if p.closeallfigure~=1
    figure;
    showMatchedFeatures(boxImage, sceneImage, matchedBoxPoints, ...
        matchedScenePoints, 'montage');
    title('Putatively Matched Points (Including Outliers)');
end

% Step 5: Locate the Object in the Scene Using Putative Matches
[tform, inlierBoxPoints, inlierScenePoints] = ...
    estimateGeometricTransform(matchedBoxPoints, matchedScenePoints, 'affine');
if p.closeallfigure~=1
    figure;
    showMatchedFeatures(boxImage, sceneImage, inlierBoxPoints, ...
        inlierScenePoints, 'montage');
    title('Matched Points (Inliers Only)');
end

%{
boxPolygon = [1, 1;...                           % top-left
        size(boxImage, 2), 1;...                 % top-right
        size(boxImage, 2), size(boxImage, 1);... % bottom-right
        1, size(boxImage, 1);...                 % bottom-left
        1, 1];                   % top-left again to close the polygon
%}
boxPolygon=[50,-5;105,-5;160,-5;190,45;220,95;190,145;160,195;85,195;50,195;25,145;0,95;25,45;50,-5];
    
newBoxPolygon = transformPointsForward(tform, boxPolygon);
if p.closeallfigure~=1
figure;
    imshow(sceneImage);
    hold on;
    line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'y');
    title('Detected Box');
end
% Website link:
% http://www.mathworks.com/help/vision/examples/object-detection-in-a-cluttered-scene-using-point-feature-matching.html?prodcode=VP&language=en

end
function ObjectDetection_MSER(sceneImage,boxImage)

% Step 1: Load images
%oad matlab.mat;


% Step 2: Detect Feature Points

boxregions = detectMSERFeatures(boxImage);
sceneregions = detectMSERFeatures(sceneImage);
%[boxregions,boxcomp] = detectMSERFeatures(boxImage);
%[sceneregions,scenecomp] = detectMSERFeatures(sceneImage);
%boxStats = regionprops(boxregions, 'BoundingBox', 'Eccentricity', ...
%    'Solidity', 'Extent', 'Euler', 'Image');
%sceneStats = regionprops(sceneregions, 'BoundingBox', 'Eccentricity', ...
%    'Solidity', 'Extent', 'Euler', 'Image');

figure;
imshow(boxImage);
title('100 Strongest Feature Points from Box Image');
hold on;
plot(boxregions);

figure;
imshow(sceneImage);
title('300 Strongest Feature Points from Scene Image');
hold on;
plot(sceneregions);

%figure
%imshow(boxImage)
%hold on
%plot(boxregions, 'showPixelList', true,'showEllipses',false)
%title('MSER regions')
%hold off


% Step 3: Extract Feature Descriptors
[boxFeatures, boxregions] = extractFeatures(boxImage, boxregions);
[sceneFeatures, sceneregions] = extractFeatures(sceneImage, sceneregions);


% Step 4: Find Putative Point Matches
boxPairs = matchFeatures(boxFeatures, sceneFeatures);
matchedBoxPoints = boxPoints(boxPairs(:, 1), :);
matchedScenePoints = scenePoints(boxPairs(:, 2), :);
figure;
showMatchedFeatures(boxImage, sceneImage, matchedBoxPoints, ...
    matchedScenePoints, 'montage');
title('Putatively Matched Points (Including Outliers)');



end
function ObjectDetection_Cascade(sceneImage,boxImage)

close all; 
clear all; 
clc; 

load('stopSigns.mat'); 
imDir = fullfile(matlabroot,'toolbox','vision','visiondemos','stopSignImages'); 
addpath(imDir); 

negativeFolder = fullfile(matlabroot,'toolbox','vision','visiondemos','non_stop_signs'); 
trainCascadeObjectDetector('stopSignDetector.xml',data,negativeFolder,'FalseAlarmRate',0.2,'NumCascadeStages',5); 

detector = vision.CascadeObjectDetector('stopSignDetector.xml'); 
img = imread('stopSignTest.jpg'); 
bbox = step(detector,img); 
detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign'); 
figure;
imshow(detectedImg); 
end


%% GUI unit open functions
% --- Executes just before eyesmFISH_GUI is made visible.
function eyesmFISH_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eyesmFISH_GUI (see VARARGIN)

% Choose default command line output for eyesmFISH_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes eyesmFISH_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%setappdata(My_GUI_handle, 'IgnoreCloseAll', 1)
end

% --- Outputs from this function are returned to the command line.
function varargout = eyesmFISH_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% show axes image
%tiff_numshow=['/',num2str(size(tiff_images,2))];
%set(handles.text5,'String',tiff_numshow);

load([handles.datafolder 'data.mat'],'p','tiff_images');

tiff_num=str2num(get(handles.edit2,'String'));
tiff_chal=get(handles.popupmenu1,'Value');

myImage = tiff_images{1,tiff_num}(:,:,tiff_chal);
set(handles.axes1,'Units','pixels');
resizePos = get(handles.axes1,'Position');
myImage= imresize(myImage, [resizePos(3) resizePos(3)]);
axes(handles.axes1);
imagesc(myImage);
set(handles.axes1,'Units','normalized');
% Update handles
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.

%load([handles.datafolder 'data.mat'],'p','tiff_images');
%tiff_nummax=size(tiff_images,2);

%set(handles.slider1,'max',tiff_nummax);set(handles.slider1,'min',1);    
%tiff_num=get(handles.slider1,'Value');tifftxt=[num2str(tiff_num) ' stack'];
    
%set(hObject,'Max',tt-1,'Min',1)
%set(hObject,'SliderStep',[1/(tt-1) 0.1])
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);   
end
end
