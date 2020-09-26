classdef TimeSeriesAnalyzer < handle
    %Flexible and performant time series viewer and analysis tool.
    %
    %   -------------------------------------------------------------------
    %   Usage
    %   -------------------------------------------------------------------
    %   tsa = TimeSeriesAnalyzer;
    %   tsa.setData(...);
    %   Now you have command line access to everything via tsa.Data,
    %   tsa.Measurement, etc.
    %
    %   -------------------------------------------------------------------
    %   Data for one or more time series.
    %   -------------------------------------------------------------------
    %   Data -> struct array of time series (x,y) data
    %       Data.xdata -> Nx1 column vector
    %       Data.ydata -> Nx1 column vector
    %       Data.xlabel -> string, e.g. "Time, s"
    %       Data.ylabel -> string, e.g. "Current, pA"
    %       ...
    %   setData(...) <- !!! use this to be safe
    %
    %   -------------------------------------------------------------------
    %   Visualize multiple time series with arbitrary grouping (e.g. for
    %   multiple channels in a recording). Each group is shown in a
    %   separate plot axes. You can display a selected subset of groups.
    %   e.g. setVisibleGroups([2 5]) -> show 2nd and 5th groups only.
    %   -------------------------------------------------------------------
    %   groups() -> group indexes of each time series in Data
    %   setGroups(...)
    %
    %   visibleGroups() -> groups visible in the UI
    %   setVisibleGroups(...)
    %
    %   -------------------------------------------------------------------
    %   Display a selected subset of time series within each group.
    %   e.g. setVisibleSweeps([1 3]) -> show 1st and 3rd time series in each group
    %   -------------------------------------------------------------------
    %   visibleSweeps() -> subset of time series to show in each group
    %   setVisibleSweeps(...)
    %
    %   -------------------------------------------------------------------
    %   Each time series in Data can have multiple associated time series
    %   stored in the fields Data.yNAME that are refered to by NAME.
    %   The basic required field is ydata, but any number of other yNAME
    %   fields can exist. For example, yideal may refer to an idealization
    %   of ydata, and it is logical that it is grouped together in the
    %   same struct as the ydata it refers to. You can display a selected
    %   subset of NAMEs for each time series. Each plotted NAME will have
    %   its own color.
    %   e.g. setVisibleNames(["raw","data"]) -> show Data.yraw and Data.ydata
    %   -------------------------------------------------------------------
    %   names() -> ordered list of NAME from all Data.yNAME
    %   visibleNames() -> subset of NAMEs to display
    %   setVisibleNames(...)
    %
    %   -------------------------------------------------------------------
    %   Updating the UI
    %   -------------------------------------------------------------------
    %   .refresh() -> updates entire UI (e.g. if groups are changed)
    %   .replot() -> updates all displayed plots (e.g. if visible time series are changed)
    %
    %   -------------------------------------------------------------------
    %   X-axis range selections as interactive graphical ROIs
    %   -------------------------------------------------------------------
    %   Clicking the XROIsButton [] will toggle the active state of the
    %   XROIMgr (instance of XAxisROIManager). While active you can click
    %   and drag in the plot axes to draw x-axis range ROIs or drag/resize
    %   current ROIs. Each ROI has its own context menu. These ROIs define
    %   a set of selected ranges wrt the x-axis. Measurements, curve
    %   fitting and operatiions (see below) will pay attention to these
    %   ROIs if they are active.
    %
    %   -------------------------------------------------------------------
    %   Measurements (see plot axes context menu)
    %   -------------------------------------------------------------------
    %   Measurement -> table of results for the last measurement
    %   measure(axes, method) -> performs the measurement defined by
    %       method on all visible plots in axes. If XROIMgr is active,
    %       measurement will be restricted to each ROI for each plot.
    %
    %   -------------------------------------------------------------------
    %   Curve fitting (see plot axes context menu)
    %   -------------------------------------------------------------------
    %   curveFit(axes, methodOrExpr) -> fits each Data.ydata visible
    %       in axes according to methodOrExpr (e.g. mean, line, polynomial,
    %       spline, a*exp(-x/b)+c). The resulting fits are stored in Data.yfit.
    %       If XROIMgr is active, the fit error will only include data
    %       points within the ROIs. There is also the option to limit the
    %       fit to within a single ROI. !!! This will overwrite Data.yfit.
    %   Fit plot objects have a context menu for deleting the fit or using
    %   it to baseline or normalize its associated ydata.
    %
    %   -------------------------------------------------------------------
    %   Mask, zero, interpolate, set, add, subtract, multiply, divide...
    %   (see plot axes context menu)
    %   -------------------------------------------------------------------
    %   operation(axes, operation) -> apply operation to all visible
    %       plots in axes (except for yraw). A dialog will ask whether or
    %       not to save the operation. If XROIMgr is active, the operation
    %       will only be applied to data points within the ROIs.
    %
    %   -------------------------------------------------------------------
    %   Really large time series
    %   -------------------------------------------------------------------
    %   Fast performance for extremely large plots via Jim Hokanson's
    %   plotBig_Matlab (find it in the Add-On Explorer).
    %
    %   Buttons exist for paging left/right while zoomed.
    %
    %   -------------------------------------------------------------------
    %   Author: Marcel Goldschen-Ohm
    %   -------------------------------------------------------------------
    
    properties
        % Time series data (struct array over multiple time series)
        %{
        ///////////////////////////////////////////////////////////////////
        //                     REQUIRED FIELDS
        ///////////////////////////////////////////////////////////////////
        
        -------------------------------------------------------------------
        Time series (x,y) data
        -------------------------------------------------------------------
        .xdata            Nx1 column data
        .ydata            Nx1 column data
        
        .xlabel           string, e.g. "Time, s"
        .ylabel           string, e.g. "Current, pA"
        
        ///////////////////////////////////////////////////////////////////
        //                     OPTIONAL FIELDS
        //            (these have functionality in the UI)
        ///////////////////////////////////////////////////////////////////
        
        -------------------------------------------------------------------
        Group index for arbitrary grouping of time series.
        A default option is to group according to ylabel.
        For example, data might alternate between multiple channels.
        Each group has its own plot in the viewer.
        -------------------------------------------------------------------
        .group            index
        
        -------------------------------------------------------------------
        Associated time series that you want to plot overlaid on .ydata.
        e.g. An idealization.
        !!! Replace NAME with whatever => .ymydata, .yMyOtherData
        In the viewer you can select by NAME which data to view.
        -------------------------------------------------------------------
        .xNAME            Nx1 (optional, defaults to .xdata)
        .yNAME            Nx1 OR cfit object OR piecewise polynomial struct
        
        e.g.
        .yraw             Nx1 copy of original unmodified .ydata
        .yideal           Nx1 idealization of .ydata
        .yfit             Nx1 OR cfit object OR piecewise polynomial struct
        .ybaseline        Nx1 baseline for .ydata
        
        ///////////////////////////////////////////////////////////////////
        //                      OTHER FIELDS
        ///////////////////////////////////////////////////////////////////
        
        It's just a struct array, so put whatever you want in it.
        But it's up to you to add functionality in the UI for any of it.
        
        %}
        % !!! WARNING !!! You should use setData() to set this unless you
        %                 know what you are doing.
        Data = struct('xdata', [], 'ydata', [], 'xlabel', "", 'ylabel', "");
        
        % User interface graphics objects
        ui = struct();
        
        % Attempt to use plotBig instead of plot if >= BigN points.
        % See plotBig_Matlab by Jim Hokanson in Add-On Explorer or at
        % https://github.com/JimHokanson/blog/tree/master/2018_01_PlotBig_Matlab
        BigN = 1e6;
        
        % Manager for x-axis range ROIs (instance of XAxisROIManager)
        XROIMgr
        
        % Set by measure method (table of measurements)
        Measurement = [];
    end
    
    properties (Dependent)
        % parent graphics object of .UI.MainPanel
        Parent
    end

    properties (Access = protected)
        % Visible group indices, [] -> all visible
        % e.g. [1,3] => show groups 1 and 3
        % !!! Do NOT access directly, get via .visibleGroups()
        VisibleGroups = [];
        
        % Visible within group time series (sweeps) indices, [] -> all visible
        % e.g. [1,3] => show the 1st and 3rd time series in each group
        % !!! Do NOT access directly, get via .visibleWithinGroupTraces()
        VisibleSweeps = [];
        
        % Visible names refer to .yNAME fields, [] -> all visible
        % e.g. ["data","ideal"] => show .ydata and .yideal
        % !!! Do NOT access directly, get via .visibleNames()
        VisibleNames = string.empty;
    end
    
    methods
        % constructor/destructor
        function obj = TimeSeriesAnalyzer(parent)
            % Parent container for the main panel.
            % If none is supplied, create a figure to hold the main panel.
            if ~exist('parent', 'var') || ~isvalid(parent)
                parent = figure( ...
                    'Name', 'Time Series Analyzer', ...
                    'NumberTitle', 0, ...
                    'Menubar', 'none', ...
                    'Toolbar', 'none', ...
                    'AutoResizeChildren', 0);
                didCreateOuterFigure = true;
            else
                didCreateOuterFigure = false;
            end
            
            % Get handle to outer figure container.
            obj.ui.Figure = ancestor(parent, {'Figure', 'matlab.ui.Figure'}, 'toplevel');
            
            % Everyting will be inside the main panel.
            % This way we can reparent the whole UI by just reparenting the
            % main panel. Similar for toggling visiblity, etc.
            % This also makes it easy to place this UI inside another UI.
            obj.ui.MainPanel = uipanel(parent, ...
                'BorderType', 'none', ...
                'Units', 'normalized', ...
                'Position', [0 0 1 1], ...
                'AutoResizeChildren', 0, ...
                'DeleteFcn', @(varargin) delete(obj));
            
            % plot axes
            [obj.ui.TsAxes, obj.ui.HistAxes] = obj.createAxes();
            
            % x-axis ROI manager
            try
                obj.XROIMgr = XAxisROIManager(obj.ui.TsAxes);
            catch
                obj.XROIMgr = [];
            end
            
            % toolbar above axes
            obj.ui.ToolbarPanel = obj.createToolbarPanel();
            
            % add figure menubar only if the UI is not part of another UI
            if didCreateOuterFigure
                obj.createMenubar();
            end
            
            % plot big?
            try
                h = plotBig(obj.ui.TsAxes(1), 1:obj.BigN, 1:obj.BigN);
                delete(h);
                obj.ui.HasPlotBig = true;
            catch
                obj.ui.HasPlotBig = false;
            end
            
            % on main panel resize
            obj.ui.MainPanel.SizeChangedFcn = @(varargin) obj.resize();
            
            % although parent is already set, this sets up some listeners
            obj.Parent = obj.Parent;
            
            % refresh the UI
            obj.refresh();
        end
        function delete(obj)
            if isfield(obj.ui, 'menubar')
                delete(obj.ui.Menubar);
            end
            if ~obj.ui.MainPanel.BeingDeleted
                delete(obj.ui.MainPanel);
            end
        end
        
        %------------------------------------------------------------------
        % Data
        %------------------------------------------------------------------
        
        % Set the time series data
        function setData(obj, data)
            if isstruct(data)
                % expected to have appropriate format
                obj.Data = data;
            elseif iscell(data)
                % cell array of [x y] or [y]
                obj.Data = struct.empty;
                ndata = numel(data);
                for i = 1:ndata
                    if size(data{i}, 1) == 1 || size(data{i}, 2) == 1
                        % y or y'
                        obj.Data(i).xdata = reshape(1:length(data{i}), [], 1);
                        obj.Data(i).ydata = reshape(data{i}, [], 1);
                    else
                        % [x y]
                        obj.Data(i).xdata = data{i}(:,1);
                        obj.Data(i).ydata = data{i}(:,2);
                    end
                end
            else
                % [x y] or [y]
                obj.Data = struct;
                if size(data, 1) == 1 || size(data, 2) == 1
                    % y or y'
                    obj.Data.xdata = reshape(1:length(data), [], 1);
                    obj.Data.ydata = reshape(data, [], 1);
                else
                    % [x y]
                    obj.Data.xdata = data(:,1);
                    obj.Data.ydata = data(:,2);
                end
            end
            if ~isfield(obj.Data, 'xlabel')
                [obj.Data.xlabel] = deal("");
            end
            if ~isfield(obj.Data, 'ylabel')
                [obj.Data.ylabel] = deal("");
            end
            % group by ylabel
            obj.setGroups(obj.groupsByYLabel());
            % copy of original raw data
            if ~isfield(obj.Data, 'yraw')
                [obj.Data.yraw] = deal([]);
                for i = 1:numel(obj.Data)
                    obj.Data(i).yraw = obj.Data(i).ydata;
                end
            end
            % show only first trace in each group
            obj.setVisibleSweeps(1);
            % show all yNAME except yraw
            obj.setVisibleNames(string.empty);
            obj.showAllGroups(); % calls refresh()
            obj.autoscaleXY();
        end
        
        % Set the sample interval for all time series
        function setSampleInterval(obj, dx)
            if ~exist('dx', 'var')
                answer = inputdlg({'sample interval '}, 'Sample Interval', 1);
                if isempty(answer)
                    return
                end
                dx = str2num(answer{1});
            end
            for i = 1:numel(obj.Data)
                n = length(obj.Data(i).ydata);
                obj.Data(i).xdata = reshape((0:n-1) .* dx, size(obj.Data(i).ydata));
            end
            obj.replot();
        end
        
        % get (x,y) values for named time series in data struct
        function [x,y] = getXY(obj, data, name)
            if ~exist('name', 'var')
                name = "data";
            else
                name = string(name);
            end
            y = data.("y"+name);
            try
                x = data.("x"+name);
            catch
                x = data.xdata;
            end
            if isempty(x)
                x = data.xdata;
            end
            if class(y) == "cfit"
                y = y(x);
            elseif isstruct(y)
                if isfield(y, 'form') && y.form == "pp"
                    y = ppval(y, x);
                end
            end
        end
        
        % Edit table of Data struct array fields.
        % Only allows editing single value numeric or string fields.
        % Array fields can't be edited, but their size is displayed.
        % TODO: Open array fields in new table to edit the arrays.
        function editData(obj)
            if numel(obj.Data) == 0
                return
            elseif numel(obj.Data) == 1
                tdata = struct2table(obj.Data, 'AsArray', true);
            else
                tdata = struct2table(obj.Data);
            end
            fig = uifigure('Name', 'Data');
            uitable(fig, 'Data', tdata, ...
                'Position', [0 0 fig.Position(3:4)], ...
                'RowName', 'numbered', ...
                'ColumnEditable', true(1, size(tdata, 2)), ...
                'CellEditCallback', @obj.dataTableCellEdited);
            % align fig to upper left of main panel
            fbox = getpixelposition(obj.ui.Figure);
            pbox = getpixelposition(obj.ui.MainPanel);
            xy = fbox(1:2) + pbox(1:2) + [0 pbox(4)];
            fig.Position(1:2) = xy - [0 fig.Position(4)];
        end
        function dataTableCellEdited(obj, src, event)
            row = event.Indices(1);
            col = event.Indices(2);
            header = src.Data.Properties.VariableNames{col};
            tableColumn = src.Data.(header);
            if iscell(tableColumn)
                % if table column is a cell array:
                % - cell must itself be a collection or array
                % - don't edit and reset table cell to original data
                src.Data.(header){row} = obj.Data(row).(header);
                return
            end
            obj.Data(row).(header) = tableColumn(row);
            obj.refresh();
        end
        
        % Revert to original raw data (yraw)
        function revertToRaw(obj, tsi)
            if ~isfield(obj.Data, 'yraw')
                return
            end
            if ~exist('tsi', 'var') || isempty(tsi)
                tsi = 1:numel(obj.Data);
            end
            for i = 1:numel(tsi)
                j = tsi(i);
                [x,y] = obj.getXY(obj.Data(j), "raw");
                if ~isempty(y) && numel(x) == numel(y)
                    obj.Data(j).xdata = x;
                    obj.Data(j).ydata = y;
                end
            end
            obj.replot();
        end
        function revertToRawAllVisibleInAxes(obj, ax)
            tsi = obj.visibleTsInAxes(ax);
            obj.revertToRaw(tsi);
        end
        
        %------------------------------------------------------------------
        % File I/O
        %------------------------------------------------------------------
        
        % Load time series data from .mat file
        function loadData(obj, filepath)
            if ~exist('filepath', 'var') || isempty(filepath)
                [file, path] = uigetfile('*.mat', 'Open data file.');
                if isequal(file, 0)
                    return
                end
                filepath = fullfile(path, file);
            end
            try
                wb = waitbar(0, 'Loading data from file...');
                tmp = load(filepath);
                close(wb);
                if isfield(tmp, 'Data')
                    obj.setData(tmp.Data);
                elseif isfield(tmp, 'data')
                    % to cover my own naming convention change
                    obj.setData(tmp.data);
                end
                [path, name, ext] = fileparts(filepath);
                obj.ui.Figure.Name = name;
            catch err
                disp(err);
            end
        end
        
        % Save time series data to .mat file
        function saveData(obj, filepath)
            if ~exist('filepath', 'var') || isempty(filepath)
                [file, path] = uiputfile('*.mat', 'Save data to file.');
                if isequal(file, 0)
                    return
                end
                filepath = fullfile(path, file);
            end
            try
                Data = obj.Data;
                wb = waitbar(0, 'Saving data to file...');
                save(filepath, 'Data');
                close(wb);
                [path, name, ext] = fileparts(filepath);
                obj.ui.Figure.Name = name;
            catch err
                disp(err);
            end
        end
        
        % Import time series data from HEKA binary file
        function importHEKA(obj, filepath)
            if ~exist('filepath', 'var') || isempty(filepath)
                [file, path] = uigetfile('*.*', 'Open HEKA data file.');
                if isequal(file, 0)
                    return
                end
                filepath = fullfile(path, file);
            end
            try
                heka = HEKA_Importer(filepath);
            catch err
                disp(err);
                warndlg("!!! Requires package 'HEKA Patchmaster Importer' by Christian Keine. Find in MATLAB's Add-On Explorer.", ...
                    'HEKA file loader');
                return
            end
            % info for each recording are in the nonempty leaves of dataTree
            recdata = heka.trees.dataTree(:,end);
            i = 1;
            while i <= numel(recdata)
                if isempty(recdata{i})
                    recdata(i) = [];
                else
                    i = i + 1;
                end
            end
            numRecordings = size(heka.RecTable,1);
            if numel(recdata) ~= numRecordings
                warndlg('Unexpected data structure. Please report this error.');
                return
            end
            if numRecordings > 1
                % Ask which recordings to load.
                stimuli = {};
                for rec = 1:numRecordings
                    stimuli{rec} = heka.RecTable.Stimulus{rec};
                end
                selrec = listdlg('ListString', stimuli, ...
                    'PromptString', 'Select recordings to load:');
            else
                selrec = 1;
            end
            % load each recording
            data = struct.empty;
            t = 0;
            for i = 1:numel(selrec)
                rec = selrec(i);
                numChannels = numel(heka.RecTable.dataRaw{rec});
                numSweeps = size(heka.RecTable.dataRaw{rec}{1},2);
                numPts = size(heka.RecTable.dataRaw{rec}{1},1);
                for s = 1:numSweeps
                    for c = 1:numChannels
                        t = t + 1;
                        data(t).xdata = recdata{rec}.TrXInterval;
                        data(t).ydata = heka.RecTable.dataRaw{rec}{c}(:,s);
                        data(t).xlabel = "Time, " + string(heka.RecTable.TimeUnit{rec}{c});
                        data(t).ylabel = string(heka.RecTable.ChName{rec}{c}) + ", " + string(heka.RecTable.ChUnit{rec}{c});
                    end
                end
            end
            obj.setData(data);
            [path, name, ext] = fileparts(filepath);
            obj.ui.Figure.Name = name;
        end
        
        %------------------------------------------------------------------
        % Groups
        %------------------------------------------------------------------
        
        % Set/get group index of each time series
        function setGroups(obj, ind)
            if ~exist('ind', 'var')
                answer = inputdlg({'Group Indices '}, 'Groups', 1, {num2str(obj.groups())});
                if isempty(answer)
                    return
                end
                ind = str2num(answer{1});
            end
            nts = numel(obj.Data);
            nind = numel(ind);
            if isempty(ind)
                ind = obj.groupsByYLabel();
            elseif nind < nts
                % repeat for all data
                ind = repmat(reshape(ind, 1, []), 1, ceil(nts / nind));
                ind = ind(1:nts);
            elseif nind > nts
                % ignore extra
                ind = ind(1:nts);
            end
            assert(numel(ind) == nts);
            for i = 1:nts
                obj.Data(i).group = ind(i);
            end
        end
        function ind = groups(obj)
            if isempty(obj.Data)
                ind = [];
                return
            end
            if ~isfield(obj.Data, 'group')
                ind = obj.groupsByYLabel();
                return
            end
            ind = horzcat(obj.Data.group);
        end
        
        % Get number of groups
        function n = numGroups(obj)
            n = numel(unique(obj.groups()));
        end
        
        % Group time series by ylabel
        function groupByYLabel(obj)
            obj.setGroups(obj.groupsByYLabel());
            obj.setVisibleSweeps(1);
            obj.showAllGroups();
        end
        function ind = groupsByYLabel(obj)
            % group by data.ylabel
            nts = numel(obj.Data);
            if ~isfield(obj.Data, 'ylabel')
                ind = ones(1, nts);
                return
            end
            ylabels = arrayfun(@(ts) string(ts.ylabel), obj.Data);
            glabels = unique(ylabels);
            % sort glabels in order they first come in ts array
            idx = zeros(size(glabels));
            for i = 1:numel(glabels)
                idx(i) = find(ylabels == glabels(i), 1);
            end
            [~,idx] = sort(idx);
            glabels = glabels(idx);
            ngroups = numel(glabels);
            firstTraceInGroup = zeros(1, ngroups);
            for i = 1:ngroups
                firstTraceInGroup(i) = find(ylabels == glabels(i), 1);
            end
            [~,gorder] = sort(firstTraceInGroup);
            ind = zeros(1, nts);
            for i = 1:ngroups
                ind(ylabels == glabels(i)) = gorder(i);
            end
        end
        
        % Merge all time series into a single group
        function mergeAll(obj)
            obj.setGroups(1);
            obj.setVisibleSweeps(1);
            obj.showAllGroups();
        end
        
        % Separate each time series into its own group
        function separateAll(obj)
            obj.setGroups(1:numel(obj.Data));
            obj.showAllGroups();
            obj.showAllSweeps();
        end
        
        % Group every Nth time series
        function groupEveryN(obj, n)
            if ~exist('n', 'var')
                answer = inputdlg({'N: '}, 'Group Every N Traces', 1);
                if isempty(answer)
                    return
                end
                n = str2num(answer{1});
            end
            obj.setGroups(1:n);
            obj.setVisibleSweeps(1);
            obj.showAllGroups();
        end
        
        % Group blocks of N time series
        function groupBlocksOfN(obj, n)
            if ~exist('n', 'var')
                answer = inputdlg({'N: '}, 'Group Blocks of N Traces', 1);
                if isempty(answer)
                    return
                end
                n = str2num(answer{1});
            end
            nts = numel(obj.Data);
            obj.setGroups(reshape(repmat(1:nts, n, 1), 1, []));
            obj.setVisibleSweeps(1);
            obj.showAllGroups();
        end
        
        % Get cell array of index arrays of time series (Data) in each group
        function tsi = tsPerGroup(obj)
            groups = obj.groups(); % group index of each time series
            groupids = unique(groups);
            ngroups = numel(groupids);
            if ngroups == 0
                tsi = {};
                return
            end
            tsi{ngroups} = [];
            for i = 1:ngroups
                tsi{i} = find(groups == groupids(i));
            end
        end
        
        % Get array of number of time series in each group
        function n = numTsPerGroup(obj)
            n = cellfun(@(ind) numel(ind), obj.tsPerGroup());
        end
        
        % Get group index associated with axes object
        function group = axesGroupIndex(obj, ax)
            visGroups = obj.visibleGroups();
            axi = find(obj.ui.TsAxes == ax, 1);
            if isempty(axi)
                axi = find(obj.ui.HistAxes == ax, 1);
            end
            group = find(visGroups == axi, 1);
        end
        
        % Indices of time series in axes (each axes corresponds to a group)
        function tsi = tsInAxes(obj, ax)
            group = obj.axesGroupIndex(ax);
            tsi = obj.tsPerGroup();
            tsi = tsi{group};
        end
        
        %------------------------------------------------------------------
        % Labels
        %------------------------------------------------------------------
        
        % Set ylabel of all time series in each group to the group label
        function setGroupLabels(obj, labels)
            if iscell(labels)
                labels = cellfun(@(s) string(s), labels);
            end
            tsi = obj.tsPerGroup();
            ngroups = numel(tsi);
            for i = 1:ngroups
                [obj.Data(tsi{i}).ylabel] = deal(labels(i));
            end
        end
        
        % Get group labels as ylabel of first time series in each group
        function labels = groupLabels(obj)
            % string array of ylabel for the first trace in each group
            tsi = obj.tsPerGroup();
            ngroups = numel(tsi);
            labels = string.empty;
            for i = 1:ngroups
                try
                    labels(i,1) = string(obj.Data(tsi{i}(1)).ylabel);
                    if strlength(labels(i)) == 0
                        labels(i) = "Group " + string(i);
                    end
                catch
                    labels(i,1) = "Group " + string(i);
                end
            end
        end
        
        % Set ylabel of all time series in group to group label
        function setGroupLabel(obj, group, label)
            tsi = obj.tsPerGroup();
            tsi = tsi{group};
            firstTraceInGroup = tsi(1);
            if ~exist('label', 'var') %|| (~ischar(label) && ~isstring(label))
                try
                    default = {char(obj.Data(firstTraceInGroup).ylabel)};
                catch
                    default = {''};
                end
                answer = inputdlg({'ylabel '}, 'ylabel', 1, default);
                if isempty(answer)
                    return
                end
                label = string(answer{1});
            end
            [obj.Data(tsi).ylabel] = deal(string(label));
            obj.replot();
        end
        
        % Set xlabel of all time series
        function setXLabel(obj, label)
            if ~exist('label', 'var')
                try
                    default = {char(obj.Data(1).xlabel)};
                catch
                    default = {''};
                end
                answer = inputdlg({'xlabel '}, 'xlabel', 1, default);
                if isempty(answer)
                    return
                end
                label = string(answer{1});
            end
            [obj.Data.xlabel] = deal(string(label));
            obj.replot();
        end
        
        % Edit table of group labels
        function editGroupLabels(obj)
            Group = obj.groupLabels();
            tdata = table(Group);
            fig = uifigure('Name', 'Group Labels');
            uitable(fig, 'Data', tdata, ...
                'Position', [0 0 fig.Position(3:4)], ...
                'RowName', 'numbered', ...
                'ColumnEditable', true(1, size(tdata, 2)), ...
                'CellEditCallback', @obj.groupLabelsTableCellEdited);
            % align fig to upper left of main panel
            fbox = getpixelposition(obj.ui.Figure);
            pbox = getpixelposition(obj.ui.MainPanel);
            xy = fbox(1:2) + pbox(1:2) + [0 pbox(4)];
            fig.Position(1:2) = xy - [0 fig.Position(4)];
        end
        function groupLabelsTableCellEdited(obj, src, event)
            row = event.Indices(1);
            col = event.Indices(2);
            header = src.Data.Properties.VariableNames{col};
            if header == "Group"
                ind = obj.tsPerGroup();
                [obj.Data(ind{row}).ylabel] = deal(string(src.Data.Group(row)));
                obj.replot();
            end
        end
        
        %------------------------------------------------------------------
        % Visible Groups
        %------------------------------------------------------------------
        
        % Set/get array of visible group indices
        function setVisibleGroups(obj, ind)
            n = obj.numGroups();
            ind(ind < 1) = 1;
            ind(ind > n) = n;
            obj.VisibleGroups = unique(ind);
            obj.refresh();
        end
        function ind = visibleGroups(obj)
            ind = obj.VisibleGroups;
            n = obj.numGroups();
            ind(ind < 1) = [];
            ind(ind > n) = [];
            if isempty(ind)
                ind = 1:n;
            end
        end
        
        % Make all groups visible
        function showAllGroups(obj)
            obj.setVisibleGroups([]);
        end
        
        % Get the number of visible groups
        function n = numVisibleGroups(obj)
            n = numel(obj.visibleGroups());
        end
        
        % Create listbox for user selection of visible groups
        function lbox = createVisibleGroupsListBox(obj, parent)
            lbox = uicontrol(parent, ...
                'Style', 'listbox', ...
                'String', obj.groupLabels(), ...
                'Min', 1, ...
                'Max', obj.numGroups() + 1, ...
                'Value', obj.visibleGroups(), ...
                'Units', 'normalized', ...
                'Position', [0, 0, 1, 1], ...
                'FontSize', 12, ...
                'Callback', @(s,e) obj.setVisibleGroups(s.Value));
        end
        
        % Select visible groups in listbox
        function editVisibleGroups(obj)
            fig = figure('Name', 'Groups', ...
                'NumberTitle', 'off', ...
                'Menubar', 'none', ...
                'Toolbar', 'none', ...
                'Position', [0 0 150 150]);
            obj.createVisibleGroupsListBox(fig);
            % align fig to upper left of main panel
            fbox = getpixelposition(obj.ui.Figure);
            pbox = getpixelposition(obj.ui.MainPanel);
            xy = fbox(1:2) + pbox(1:2) + [0 pbox(4)];
            fig.Position(1:2) = xy - [0 fig.Position(4)];
        end
        
        %------------------------------------------------------------------
        % Visible Sweeps (Within Group Time Series)
        %------------------------------------------------------------------
        
        % Set/get array of within group visible sweep indices
        function setVisibleSweeps(obj, ind)
            if ~exist('ind', 'var')
                % get ind from edit uicontrol
                s = obj.ui.VisibleSweepsEdit.String;
                s = strsplit(s, 'of');
                ind = str2num(s{1});
            end
            n = obj.numTsPerGroup();
            n = n(obj.visibleGroups());
            n = max(n);
            ind(ind < 1) = 1;
            ind(ind > n) = n;
            ind = unique(ind);
            obj.VisibleSweeps = ind;
            if isempty(ind)
                obj.ui.VisibleSweepsEdit.String = sprintf('1:%d', n);
            elseif numel(ind) == 1
                obj.ui.VisibleSweepsEdit.String = sprintf('%d of %d', ind, n);
            elseif all(diff(ind) == 1)
                obj.ui.VisibleSweepsEdit.String = sprintf('%d:%d', ind(1), ind(end));
            else
                obj.ui.VisibleSweepsEdit.String = num2str(ind);
            end
            obj.replot();
        end
        function ind = visibleSweeps(obj)
            n = obj.numTsPerGroup();
            n = n(obj.visibleGroups());
            n = max(n);
            ind = obj.VisibleSweeps;
            ind(ind < 1) = [];
            ind(ind > n) = [];
            if isempty(ind)
                ind = 1:n;
            end
        end
        
        % Make all time series within each group visible
        function showAllSweeps(obj)
            obj.setVisibleSweeps([]);
        end
        
        % Incr/decr visible time series within each group
        function prevSweep(obj)
            ind = obj.visibleSweeps();
            if isempty(ind)
                ind = max(obj.numTsPerGroup());
            elseif numel(ind) > 1
                ind = ind(end);
            else
                ind = ind - 1;
            end
            obj.setVisibleSweeps(ind);
        end
        function nextSweep(obj)
            ind = obj.visibleSweeps();
            if isempty(ind)
                ind = 1;
            elseif numel(ind) > 1
                ind = ind(1);
            else
                ind = ind + 1;
            end
            obj.setVisibleSweeps(ind);
        end
        
        % Get cell array of arrays of indices for visible time series (Data) in each group
        function tsi = visibleTsPerGroup(obj)
            visSweeps = obj.visibleSweeps();
            tsi = obj.tsPerGroup();
            ngroups = numel(tsi);
            for i = 1:ngroups
                vis = intersect(1:numel(tsi{i}), visSweeps);
                tsi{i} = tsi{i}(vis);
            end
        end
        
        % Indices of visible time series in axes
        function tsi = visibleTsInAxes(obj, ax)
            group = obj.axesGroupIndex(ax);
            tsi = obj.visibleTsPerGroup();
            tsi = tsi{group};
        end
        
        %------------------------------------------------------------------
        % Visible Names (Data.yNAME)
        %------------------------------------------------------------------
        
        % Get string array of names for all yNAME fields except ylabel
        % e.g. ["raw" "data" "ideal" "mydata" ...]
        function names = names(obj)
            % .y*
            names = string.empty;
            fnames = fieldnames(obj.Data);
            for i = 1:numel(fnames)
                if startsWith(fnames{i}, 'y')
                    names(end+1) = string(fnames{i}(2:end));
                end
            end
            
            % not .ylabel
            names(names == "label") = [];
            
            % starts with: raw, data, ideal
            first = ["raw", "data", "ideal"];
            for i = numel(first):-1:1
                if ismember(first(i), names)
                    names = [first(i) names(names ~= first(i))];
                end
            end
            
            % ends with: fit, baseline
            last = ["fit", "baseline"];
            for i = 1:numel(last)
                if ismember(last(i), names)
                    names = [names(names ~= last(i)) last(i)];
                end
            end
        end
        
        % Set/get string array of visible names refering to yNAME fields
        function setVisibleNames(obj, names)
            if iscell(names)
                names = cellfun(@(s) string(s), names);
            end
            obj.VisibleNames = names;
            obj.replot();
        end
        function names = visibleNames(obj)
            if isempty(obj.VisibleNames)
                names = obj.names();
                % by default don't show the raw data
                names(names == "raw") = [];
            else
                names = intersect(obj.VisibleNames, obj.names(), 'stable');
            end
        end
        
        % Create listbox for user selection of visible names
        function lbox = createVisibleNamesListBox(obj, parent)
            names = obj.names();
            visNames = obj.visibleNames();
            ind = sort(arrayfun(@(visName) find(names == visName, 1), visNames));
            lbox = uicontrol(parent, ...
                'Style', 'listbox', ...
                'String', names, ...
                'Min', 1, ...
                'Max', numel(names) + 1, ...
                'Value', ind, ...
                'Units', 'normalized', ...
                'Position', [0, 0, 1, 1], ...
                'FontSize', 12, ...
                'Callback', @(s,e) obj.setVisibleNames(s.String(s.Value)));
        end
        
        % Select visible names in listbox
        function editVisibleNames(obj)
            fig = figure('Name', 'Names', ...
                'NumberTitle', 'off', ...
                'Menubar', 'none', ...
                'Toolbar', 'none', ...
                'Position', [0 0 150 150]);
            obj.createVisibleNamesListBox(fig);
            % align fig to upper left of main panel
            fbox = getpixelposition(obj.ui.Figure);
            pbox = getpixelposition(obj.ui.MainPanel);
            xy = fbox(1:2) + pbox(1:2) + [0 pbox(4)];
            fig.Position(1:2) = xy - [0 fig.Position(4)];
        end
        
        %------------------------------------------------------------------
        % Figure
        %------------------------------------------------------------------
        
        % Mouse press in figure
        function windowMousePress(obj, fig, varargin)
            cp = get(fig, 'CurrentPoint');
            x = cp(1,1);
            y = cp(1,2);
            if fig.SelectionType == "normal" % left
            elseif fig.SelectionType == "alt" % right
                ax = obj.ui.TsAxes(1);
                if x > ax.Position(1) - 72 && x < ax.Position(1)
                    for i = 1:numel(obj.ui.TsAxes)
                        ax = obj.ui.TsAxes(i);
                        if y > ax.Position(2) && y < sum(ax.Position([2 4]))
                            % y axis left
                            obj.popupContextMenu(obj.yaxisContextMenu(ax));
                        end
                    end
                elseif x > ax.Position(1) && x < sum(ax.Position([1 3]))
                    h = 22;
                    for i = 1:numel(obj.ui.TsAxes)
                        ax = obj.ui.TsAxes(i);
                        if i == numel(obj.ui.TsAxes)
                            h = h + 22;
                        end
                        if y > ax.Position(2) - h && y < ax.Position(2)
                            % obj.ui.TsAxes(i) x axis bottom
                            obj.popupContextMenu(obj.xaxisContextMenu(ax));
                        end
                    end
                    
                end
            elseif fig.SelectionType == "extend" % middle
            elseif fig.SelectionType == "open" % double-click
            end
        end
        
        % Key press/release in figure
        function UNUSED_windowKeyPress(obj, fig, event)
        end
        function UNUSED_windowKeyRelease(obj, fig, event)
        end
        
        % Set uimode to brush, pan or zoom.
        % Also updates toggle button state for each mode in toolbar panel.
        function setBrushMode(obj, tf)
            brush(obj.ui.Figure, char(matlab.lang.OnOffSwitchState(tf)));
            obj.ui.BrushButton.Value = tf;
            obj.ui.PanButton.Value = false;
            obj.ui.ZoomButton.Value = false;
        end
        function setPanMode(obj, tf)
            if tf
                h = pan(obj.ui.Figure);
                h.Enable = char(matlab.lang.OnOffSwitchState(tf));
                h.ActionPostCallback = @obj.onPanned;
            else
                pan(obj.ui.Figure, matlab.lang.OnOffSwitchState(tf));
            end
            if isfield(obj.ui, 'brushButton')
                obj.ui.BrushButton.Value = false;
            end
            obj.ui.PanButton.Value = tf;
            obj.ui.ZoomButton.Value = false;
        end
        function setZoomMode(obj, tf)
            if tf
                h = zoom(obj.ui.Figure);
                h.Enable = char(matlab.lang.OnOffSwitchState(tf));
                h.ActionPostCallback = @obj.onZoomed;
            else
                zoom(obj.ui.Figure, matlab.lang.OnOffSwitchState(tf));
            end
            if isfield(obj.ui, 'brushButton')
                obj.ui.BrushButton.Value = false;
            end
            obj.ui.PanButton.Value = false;
            obj.ui.ZoomButton.Value = tf;
        end
        
        % Respond to zoom/pan
        function onZoomed(obj, varargin)
            obj.replotHistograms();
        end
        function onPanned(obj, varargin)
            obj.replotHistograms();
        end
        
        % Popup context menu under mouse
        function popupContextMenu(obj, menu)
            menu.Position = obj.ui.Figure.CurrentPoint(1,1:2);
            menu.Visible = true;
        end
        
        %------------------------------------------------------------------
        % Axes
        %------------------------------------------------------------------
        
        % Create a pair of axes for grouped time series and histograms
        function [tsAxes, histAxes] = createAxes(obj)
            tsAxes = axes(obj.ui.MainPanel, ...
                'Units', 'pixels', ...
                'Box', 1);
            tsAxes.Interactions = [];
            tsAxes.Toolbar.Visible = 0;
            tsAxes.ButtonDownFcn = @obj.tsAxesButtonDown;
            try
                prevTsAxes = obj.ui.TsAxes(end);
                tsAxes.XScale = prevTsAxes.XScale;
                tsAxes.YScale = prevTsAxes.YScale;
                tsAxes.XLim = prevTsAxes.XLim;
                tsAxes.YLim = prevTsAxes.YLim;
                tsAxes.XLimMode = prevTsAxes.XLimMode;
                tsAxes.YLimMode = prevTsAxes.YLimMode;
            catch
            end
%             h = addlistener(tsAxes.YRuler, 'MarkedClean', ...
%                 @obj.tsAxesYRulerMarkedClean);
%             setappdata(tsAxes, 'YRulerMarkedCleanListener', h);
            hold(tsAxes, 'on');
            
            histAxes = axes(obj.ui.MainPanel, ...
                'Units', 'pixels', ...
                'Box', 1, ...
                'XTick', [], ...
                'YTick', []);
            histAxes.Interactions = [];
            histAxes.Toolbar.Visible = 0;
            histAxes.ButtonDownFcn = @obj.histAxesButtonDown;
            try
                prevHistAxes = obj.ui.HistAxes(end);
                histAxes.Position(3) = prevHistAxes.Position(3);
                nbins = getappdata(prevHistAxes, 'NumBins');
                setappdata(histAxes, 'NumBins', nbins);
            catch
                histAxes.Position(3) = 100;
                setappdata(histAxes, 'NumBins', 50);
            end
            hold(histAxes, 'on');
            
            setappdata(tsAxes, 'HistAxes', histAxes);
            setappdata(histAxes, 'TsAxes', tsAxes);
        end
        
        % Button click in axes
        function tsAxesButtonDown(obj, ax, varargin)
            cp = get(ax, 'CurrentPoint');
            x = cp(1,1);
            y = cp(1,2);
            if x < ax.XLim(1) % left
                return
            elseif y < ax.YLim(1) % bottom
                return
            elseif x > ax.XLim(2) % right
                return
            elseif y > ax.YLim(2) % top
                return
            end
            if ismember(obj.ui.Figure.SelectionType, ["normal","extend"]) % left, middle (shift+left)
                if obj.hasXROIMgr() && obj.ui.XROIsButton.Value
                    % drag xrange ROI
                    if obj.ui.Figure.SelectionType == "normal"
                        delete(obj.XROIMgr.ROIs);
                        obj.XROIMgr.ROIs = gobjects(0);
                    end
                    numOriginalRois = obj.XROIMgr.numROIs();
                    obj.XROIMgr.addROI(x);
                    j = size(obj.XROIMgr.ROIs, 2);
                    if j > numOriginalRois 
                        for i = 1:size(obj.XROIMgr.ROIs, 1)
                            roi = obj.XROIMgr.ROIs(i,j);
                            roiAxes = obj.ui.TsAxes(i);
                            roi.UIContextMenu = obj.tsAxesContextMenu(roiAxes, roi);
                        end
                    end
                    obj.updateXROIsButton();
                end
            elseif obj.ui.Figure.SelectionType == "alt" % right
                obj.popupContextMenu(obj.tsAxesContextMenu(ax));
            elseif obj.ui.Figure.SelectionType == "extend" % middle (shift+left)
            elseif obj.ui.Figure.SelectionType == "open" % double-click
            end
        end
        function histAxesButtonDown(obj, ax, varargin)
            cp = get(ax, 'CurrentPoint');
            x = cp(1,1);
            y = cp(1,2);
            if x < ax.XLim(1) % left
                return
            elseif y < ax.YLim(1) % bottom
                return
            elseif x > ax.XLim(2) % right
                return
            elseif y > ax.YLim(2) % top
                return
            end
            if obj.ui.Figure.SelectionType == "normal" % left
            elseif obj.ui.Figure.SelectionType == "alt" % right
                obj.popupContextMenu(obj.histAxesContextMenu(ax));
            elseif obj.ui.Figure.SelectionType == "extend" % middle
            elseif obj.ui.Figure.SelectionType == "open" % double-click
            end
        end
        
        % Page tsAxes.XLim left or right
        function pageLeft(obj)
            dataXlim = [];
            for i = 1:numel(obj.ui.TsAxes)
                h = findobj(obj.ui.TsAxes(i), 'Type', 'line');
                h(~isvalid(h)) = [];
                for j = 1:numel(h)
                    dataXlim = [dataXlim; h(j).XData(1) h(j).XData(end)];
                end
            end
            dataXlim = [min(dataXlim(:,1)) max(dataXlim(:,2))];
            axXlim = obj.ui.TsAxes(1).XLim;
            if dataXlim(1) < axXlim(1)
                dx = min(axXlim(1) - dataXlim(1), diff(axXlim));
                [obj.ui.TsAxes.XLim] = deal(axXlim - dx);
            end
        end
        function pageRight(obj)
            dataXlim = [];
            for i = 1:numel(obj.ui.TsAxes)
                h = findobj(obj.ui.TsAxes(i), 'Type', 'line');
                h(~isvalid(h)) = [];
                for j = 1:numel(h)
                    dataXlim = [dataXlim; h(j).XData(1) h(j).XData(end)];
                end
            end
            dataXlim = [min(dataXlim(:,1)) max(dataXlim(:,2))];
            axXlim = obj.ui.TsAxes(1).XLim;
            if dataXlim(2) > axXlim(2)
                dx = min(dataXlim(2) - axXlim(2), diff(axXlim));
                [obj.ui.TsAxes.XLim] = deal(axXlim + dx);
            end
        end
        
        % Set linear/log scale of trace axes
        function setXScale(obj, scale)
            [obj.ui.TsAxes.XScale] = deal(scale);
        end
        function setYScale(obj, scale, ax)
            if ~exist('ax', 'var') || isempty(ax)
                ax = obj.ui.TsAxes;
            end
            [ax.YScale] = deal(scale);
            % update yscale in associated histograms
            for i = 1:numel(ax)
                j = find(obj.ui.TsAxes == ax, 1);
                obj.ui.HistAxes(j).YScale = scale;
            end
        end
        
        % Autoscale in x and/or y.
        % If no axes given, defaults to autoscale all axes.
        function autoscaleXY(obj, ax)
            if ~exist('ax', 'var') || isempty(ax)
                % default is to autoscale all trace axes
                ax = obj.ui.TsAxes;
            end
            obj.autoscaleX(ax);
            obj.autoscaleY(ax);
        end
        function autoscaleX(obj, ax)
            if ~exist('ax', 'var') || isempty(ax)
                % default is to autoscale all trace axes
                ax = obj.ui.TsAxes;
            end
            % Setting XLimMode to auto does not play well with linked axes.
            % Thus, we manually determine the correct XLim.
            xlim = obj.visibleDataGlobalXLim(ax);
            if ~isempty(xlim)
                [obj.ui.TsAxes.XLim] = deal(xlim);
            end
            obj.replotHistograms();
        end
        function autoscaleY(obj, ax)
            if ~exist('ax', 'var') || isempty(ax)
                % default is to autoscale all trace axes
                ax = obj.ui.TsAxes;
            end
            [ax.YLimMode] = deal('auto');
            obj.replotHistograms();
        end
        
        % Scale all tsAxes y axis to match input axes
        function scaleAllYToMatch(obj, ax)
            [obj.ui.TsAxes.YLim] = deal(ax.YLim);
            obj.replotHistograms();
        end
        
        % Get XLim that includes all visible lines in all axes
        function xlim = visibleDataGlobalXLim(obj, ax)
            if ~exist('ax', 'var') || isempty(ax)
                ax = obj.ui.TsAxes;
            end
            xlim = [];
            for i = 1:numel(ax)
                h = findobj(ax(i), 'Type', 'line');
                h(~isvalid(h)) = [];
                for j = 1:numel(h)
                    xlim = [xlim; h(j).XData(1) h(j).XData(end)];
                end
            end
            if ~isempty(xlim)
                xlim = [min(xlim(:,1)) max(xlim(:,2))];
            end
        end
        
        % menus
        function menu = xscaleMenu(obj, parent, ax)
            menu = uimenu(parent, 'Text', 'XScale');
            isax = exist('ax', 'var') && ~isempty(ax) && isvalid(ax);
            uimenu(menu, 'Text', 'linear', ...
                'Checked', isax && ax.XScale == "linear", ...
                'MenuSelectedFc', @(varargin) obj.setXScale('linear'));
            uimenu(menu, 'Text', 'log', ...
                'Checked', isax && ax.XScale == "log", ...
                'MenuSelectedFc', @(varargin) obj.setXScale('log'));
        end
        function menu = yscaleMenu(obj, parent, ax)
            menu = uimenu(parent, 'Text', 'YScale');
            isax = exist('ax', 'var') && ~isempty(ax) && isvalid(ax);
            if isax
                uimenu(menu, 'Text', 'linear', ...
                    'Checked', ax.YScale == "linear", ...
                    'MenuSelectedFc', @(varargin) obj.setYScale('linear', ax));
                uimenu(menu, 'Text', 'log', ...
                    'Checked', ax.YScale == "log", ...
                    'MenuSelectedFc', @(varargin) obj.setYScale('log', ax));
            end
            uimenu(menu, 'Text', 'all linear', ...
                'Separator', isax, ...
                'MenuSelectedFc', @(varargin) obj.setYScale('linear'));
            uimenu(menu, 'Text', 'all log', ...
                'MenuSelectedFc', @(varargin) obj.setYScale('log'));
        end
        function menu = autoscaleMenu(obj, parent, ax)
            menu = uimenu(parent, 'Text', 'Autoscale');
            isax = exist('ax', 'var') && ~isempty(ax) && isvalid(ax);
            if isax
                uimenu(menu, 'Text', 'X', ...
                    'MenuSelectedFc', @(varargin) obj.autoscaleX(ax));
                uimenu(menu, 'Text', 'Y', ...
                    'MenuSelectedFc', @(varargin) obj.autoscaleY(ax));
                uimenu(menu, 'Text', 'XY', ...
                    'MenuSelectedFc', @(varargin) obj.autoscaleXY(ax));
                uimenu(menu, 'Text', 'Y -> all Y', ...
                    'Separator', true, ...
                    'MenuSelectedFc', @(varargin) obj.scaleAllYToMatch(ax));
            end
            uimenu(menu, 'Text', 'all X', ...
                'Separator', isax, ...
                'MenuSelectedFc', @(varargin) obj.autoscaleX());
            uimenu(menu, 'Text', 'all Y', ...
                'MenuSelectedFc', @(varargin) obj.autoscaleY());
            uimenu(menu, 'Text', 'all XY', ...
                'MenuSelectedFc', @(varargin) obj.autoscaleXY());
        end
        
        % context menus
        function menu = xaxisContextMenu(obj, ax)
            menu = uicontextmenu(obj.ui.Figure);
            obj.autoscaleMenu(menu, ax);
            uimenu(menu, 'Text', 'Set XLabel', ...
                'MenuSelectedFc', @(varargin) obj.setXLabel());
            fileMenu = obj.fileMenu(menu);
            fileMenu.Separator = true;
            obj.editMenu(menu);
            obj.groupMenu(menu);
            obj.viewMenu(menu, ax);
        end
        function menu = yaxisContextMenu(obj, ax)
            groupid = obj.axesGroupIndex(ax);
            isgroup = ~isempty(groupid) && groupid;
            menu = uicontextmenu(obj.ui.Figure);
            obj.autoscaleMenu(menu, ax);
            if isgroup
                uimenu(menu, 'Text', 'Set Group YLabel', ...
                    'MenuSelectedFc', @(varargin) obj.setGroupLabel(groupid));
            end
            fileMenu = obj.fileMenu(menu);
            fileMenu.Separator = true;
            obj.editMenu(menu);
            obj.groupMenu(menu);
            obj.viewMenu(menu, ax);
        end
        function menu = tsAxesContextMenu(obj, ax, roi)
            menu = uicontextmenu(obj.ui.Figure);
            group = obj.axesGroupIndex(ax);
            isroi = exist('roi', 'var') && isvalid(roi);
            
            obj.measureMenu(menu, ax);
            if isroi
                obj.fitMenu(menu, ax, roi);
            else
                obj.fitMenu(menu, ax);
            end
            obj.operationMenu(menu, ax);
            
            % Idealization
            methods = {'DISC', 'MDL', 'ChangePoint'};
            idealMenu = uimenu(menu, 'Text', 'Idealize', ...
                'Separator', true);
            idealGroupMenu = uimenu(idealMenu, 'Text', 'Idealize Group');
            for k = 1:numel(methods)
                uimenu(idealGroupMenu, 'Text', methods{k}, ...
                    'MenuSelectedFc', @(varargin) obj.idealizeGroup(methods{k}, group));
            end
            for k = 1:numel(methods)
                uimenu(idealMenu, 'Text', methods{k}, ...
                    'Separator', k == 1, ...
                    'MenuSelectedFc', @(varargin) obj.idealizeInAxes(methods{k}, ax));
            end
            uimenu(menu, 'Text', 'Join Ideal', ...
                ...%'Accelerator', 'J', ...
                'MenuSelectedFc', @(varargin) obj.joinIdeal(ax));
            
            % xrangesMgr ROI
            if obj.numXROIs() && (obj.ui.XROIsButton.Value || isroi)
                roisMenu = uimenu(menu, 'Text', 'ROIs', ...
                    'Separator', true);
                if isroi
                    uimenu(roisMenu, 'Text', 'Delete This ROI', ...
                        'MenuSelectedFc', @(varargin) obj.XROIMgr.deleteROI(roi));
                end
                uimenu(roisMenu, 'Text', 'Clear All ROIs', ...
                    'Separator', isroi, ...
                    'MenuSelectedFc', @(varargin) obj.XROIMgr.clearROIs());
                uimenu(roisMenu, 'Text', 'Merge Overlapping ROIs', ...
                    'MenuSelectedFc', @(varargin) obj.mergeOverlappingXROIs());
            end
            
            % Revert to raw
            uimenu(menu, 'Text', 'Revert to Raw', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.revertToRawAllVisibleInAxes(ax));
        end
        function menu = histAxesContextMenu(obj, ax)
            menu = uicontextmenu(obj.ui.Figure);
            if numel(obj.ui.HistAxes) > 1
                uimenu(menu, 'Text', 'Set # Bins (This Group)', ...
                    'MenuSelectedFc', @(varargin) obj.setHistogramNumBins([], ax));
                uimenu(menu, 'Text', 'Set # Bins (All Groups)', ...
                    'MenuSelectedFc', @(varargin) obj.setHistogramNumBins([], obj.ui.HistAxes));
            else
                uimenu(menu, 'Text', 'Set # Bins', ...
                    'MenuSelectedFc', @(varargin) obj.setHistogramNumBins([], ax));
            end
            uimenu(menu, 'Text', 'Set Axes Width', ...
                'MenuSelectedFc', @(varargin) obj.setHistogramAxesWidth());
        end
        
        %------------------------------------------------------------------
        % Plots
        %------------------------------------------------------------------
        
        % Create new or update existing line plot or histogram
        function h = createPlot(obj, ax, x, y)
            n = length(y);
            if obj.ui.HasPlotBig && n >= obj.BigN
                try
                    h = plotBig(ax, x, y, ...
                        'HitTest', 'off', ...
                        'PickableParts', 'none');
                    % Need to force rerender probably because linkaxes
                    % confuses initial rendering.
                    bp = getappdata(h, 'big_plot__data_pointer');
                    bp.big_plot_ref.forceRerender();
                    return
                catch
                end
            end
            h = line(ax, x, y, ...
                'HitTest', 'off', ...
                'PickableParts', 'none');
        end
        function h = updatePlot(obj, h, x, y)
            ax = h.Parent;
            n = length(y);
            if obj.ui.HasPlotBig
                bp = getappdata(h, 'big_plot__data_pointer');
                if ~isempty(bp)
                    if n < obj.BigN
                        % convert back to normal plot
                        delete(h);
                        h = obj.createPlot(ax, x, y);
                        return
                    end
                    % update big plot
                    bp.big_plot_ref.data.x{1} = x;
                    bp.big_plot_ref.data.y{1} = y;
                    bp.big_plot_ref.forceRerender();
                    return
                elseif n >= obj.BigN
                    % convert to big plot
                    delete(h);
                    h = obj.createPlot(ax, x, y);
                    return
                end
            end
            h.XData = x;
            h.YData = y;
        end
        function h = createHistogram(obj, ax, y)
            y = y(~isnan(y));
            nbins = getappdata(ax, 'NumBins');
            tax = getappdata(ax, 'TsAxes');
            edges = linspace(tax.YLim(1), tax.YLim(2), nbins + 1);
            centers = (edges(1:end-1) + edges(2:end)) / 2;
            counts = histcounts(y, edges);
            h = barh(ax, centers, counts, ...
                'BarWidth', 1, ...
                'LineStyle', 'none', ...
                'FaceAlpha', 0.25, ...
                'HitTest', 'off', ...
                'PickableParts', 'none');
%             ax.XTick = [];
%             ax.YTick = [];
        end
        function h = updateHistogram(obj, h, y)
            ax = h.Parent;
            y = y(~isnan(y));
            nbins = getappdata(ax, 'NumBins');
            tax = getappdata(ax, 'TsAxes');
            edges = linspace(tax.YLim(1), tax.YLim(2), nbins + 1);
            centers = (edges(1:end-1) + edges(2:end)) / 2;
            counts = histcounts(y, edges);
            h.XData = centers;
            h.YData = counts;
        end
        
        % Get (x,y) data from plot object
        function [x,y] = getXYFromPlot(obj, h)
            if obj.ui.HasPlotBig
                bp = getappdata(h, 'big_plot__data_pointer');
                if ~isempty(bp)
                    x = bp.big_plot_ref.data.x{1};
                    y = bp.big_plot_ref.data.y{1};
                    if class(x) == "big_plot.time"
                        n = length(y);
                        x = reshape((0:n-1) .* x.dt, [], 1);
                    else
                        x = reshape(x, [], 1);
                    end
                    y = reshape(y, [], 1);
                    return
                end
            end
            x = reshape(h.XData, [], 1);
            y = reshape(h.YData, [], 1);
        end
        function [x,y] = getVisibleXYFromPlot(obj, h)
            ax = h.Parent;
            [x,y] = obj.getXYFromPlot(h);
            [x,y] = obj.visibleXY(x, y, ax.XLim, ax.YLim);
        end
        
        % Get visible data points within limits
        function [x,y] = visibleXY(obj, x, y, xlim, ylim)
            if exist('xlim', 'var') && ~isempty(xlim)
                first = find(x >= xlim(1), 1);
                if isempty(first)
                    x = [];
                    y = [];
                    return
                end
                afterLast = find(x(first:end) > xlim(2), 1);
                if isempty(afterLast)
                    idx = first:length(x);
                elseif afterLast == 1
                    x = [];
                    y = [];
                    return
                else
                    last = first - 1 + afterLast - 1;
                    idx = first:last;
                end
                x = x(idx);
                y = y(idx);
            end
            if exist('ylim', 'var') && ~isempty(ylim)
                idx = find((y >= ylim(1)) & (y <= ylim(2)));
                x = x(idx);
                y = y(idx);
            end
        end
        
        % Data <--> plot
        function updateDataFromPlot(obj, h)
            tsi = getappdata(h, 'TsIndex');
            name = getappdata(h, 'TsName');
            [x,y] = obj.getXYFromPlot(h);
            obj.Data(tsi).("y"+name) = y;
            if isfield(obj.Data, "x"+name)
                obj.Data(tsi).("x"+name) = x;
            end
        end
        function updatePlotFromData(obj, h)
            tsi = getappdata(h, 'TsIndex');
            name = getappdata(h, 'TsName');
            [x,y] = obj.getXY(obj.Data(tsi), name);
            ax = h.Parent;
            allh = getappdata(ax, 'TsPlots');
            idx = find(allh == h, 1);
            allh(idx) = obj.updatePlot(h, x, y);
            setappdata(ax, 'TsPlots', allh);
        end
        
        % Find plot handles
        function h = getPlot(obj, ax, tsi, name)
            hplots = getappdata(ax, 'TsPlots');
            h = gobjects(0);
            for i = 1:numel(hplots)
                itsi = getappdata(hplots(i), 'TsIndex');
                iname = getappdata(hplots(i), 'TsName');
                if itsi == tsi && iname == name
                    h = hplots(i);
                    return
                end
            end
        end
        
        %------------------------------------------------------------------
        % Histograms
        %------------------------------------------------------------------
        
        % Show/hide histograms
        function toggleHistogramVisiblity(obj)
            [obj.ui.HistAxes.Visible] = deal(~obj.ui.HistAxes(1).Visible);
            for i = 1:numel(obj.ui.HistAxes)
                h = obj.ui.HistAxes.Children;
                h(~isvalid(h)) = [];
                if ~isempty(h)
                    [h.Visible] = deal(obj.ui.HistAxes(i).Visible);
                end
            end
            obj.resize();
            obj.replot();
        end
        
        % Set histogram axes width
        function setHistogramAxesWidth(obj, w)
            if ~exist('w', 'var')
                w = obj.ui.HistAxes(1).Position(3);
                answer = inputdlg({'Histogram Axes Width '}, 'Hist Axes Width', 1, {num2str(w)});
                if isempty(answer)
                    return
                end
                w = str2num(answer{1});
            end
            for i = 1:numel(obj.ui.HistAxes)
                obj.ui.HistAxes(i).Position(3) = w;
            end
            obj.resize();
        end
        
        % Set number of histogram bins for axes.
        % If axes is not specified, apply to all histogram axes.
        function setHistogramNumBins(obj, nbins, ax)
            if ~exist('ax', 'var') || isempty(ax)
                ax = obj.ui.HistAxes;
            end
            if ~exist('nbins', 'var') || isempty(nbins)
                nbins = getappdata(ax(1), 'NumBins');
                answer = inputdlg({'Histogram Axes Width '}, 'Hist Axes Width', 1, {num2str(nbins)});
                if isempty(answer)
                    return
                end
                nbins = str2num(answer{1});
            end
            for i = 1:numel(ax)
                setappdata(ax(i), 'NumBins', nbins);
            end
            obj.replot();
        end
        
        %------------------------------------------------------------------
        % Menubar
        %------------------------------------------------------------------
        
        % Create menubar in containing figure
        function createMenubar(obj)
            obj.ui.fileMenu = obj.fileMenu(obj.ui.Figure);
            obj.ui.editMenu = obj.editMenu(obj.ui.Figure);
            obj.ui.groupMenu = obj.groupMenu(obj.ui.Figure);
            obj.ui.viewMenu = obj.viewMenu(obj.ui.Figure);
            obj.ui.Menubar = [obj.ui.fileMenu obj.ui.editMenu obj.ui.groupMenu obj.ui.viewMenu];
        end
        
        function menu = fileMenu(obj, parent)
            menu = uimenu(parent, 'Text', 'File');
            uimenu(menu, ...
                'Text', 'New', ...
                'Accelerator', 'N', ...
                'MenuSelectedFc', @(varargin) obj.newWindow());
            uimenu(menu, ...
                'Text', 'Open', ...
                'Accelerator', 'O', ...
                'MenuSelectedFc', @(varargin) obj.loadData());
            uimenu(menu, ...
                'Text', 'Import HEKA', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.importHEKA());
            uimenu(menu, ...
                'Text', 'Save', ...
                'Accelerator', 'S', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.saveData());
        end
        function menu = editMenu(obj, parent)
            menu = uimenu(parent, 'Text', 'Edit');
            uimenu(menu, ...
                'Text', 'Edit Data Table', ...
                'MenuSelectedFc', @(varargin) obj.editData());
            uimenu(menu, ...
                'Text', 'Set Sample Interval', ...
                'MenuSelectedFc', @(varargin) obj.setSampleInterval());
            uimenu(menu, 'Text', 'Edit Group YLabels', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.editGroupLabels());
            uimenu(menu, 'Text', 'Edit XLabel', ...
                'MenuSelectedFc', @(varargin) obj.setXLabel());
            uimenu(menu, 'Text', 'Revert All to Raw', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.revertToRaw());
        end
        function menu = groupMenu(obj, parent)
            menu = uimenu(parent, 'Text', 'Group');
            uimenu(menu, 'Text', 'Show/Hide Groups', ...
                'MenuSelectedFc', @(varargin) obj.editVisibleGroups());
            uimenu(menu, 'Text', 'Edit Group YLabels', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.editGroupLabels());
            uimenu(menu, 'Text', 'Group by YLabel', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.groupTracesByYLabel());
            uimenu(menu, 'Text', 'Merge All', ...
                'MenuSelectedFc', @(varargin) obj.mergeAll());
            uimenu(menu, 'Text', 'Separate All', ...
                'MenuSelectedFc', @(varargin) obj.separateAll());
            uimenu(menu, 'Text', 'Group Every N', ...
                'MenuSelectedFc', @(varargin) obj.groupEveryN());
            uimenu(menu, 'Text', 'Group Blocks of N', ...
                'MenuSelectedFc', @(varargin) obj.groupBlocksOfN());
        end
        function menu = viewMenu(obj, parent, ax)
            menu = uimenu(parent, 'Text', 'View');
            isax = exist('ax', 'var') && ~isempty(ax) && isvalid(ax);
            if ~isax
                ax = gobjects(0);
            end
            uimenu(menu, 'Text', 'Show/Hide Data', ...
                'MenuSelectedFc', @(varargin) obj.editVisibleNames());
            submenu = obj.xscaleMenu(menu, ax);
            submenu.Separator = true;
            obj.yscaleMenu(menu, ax);
            obj.autoscaleMenu(menu, ax);
            uimenu(menu, 'Text', 'Plot Line Style/Width', ...
                'MenuSelectedFc', @(varargin) obj.showLineStyleMessageBox());
            histMenu = uimenu(menu, 'Text', 'Histograms');
            uimenu(histMenu, 'Text', 'Show/Hide Histograms', ...
                'MenuSelectedFc', @(varargin) obj.toggleHistogramVisiblity());
            uimenu(histMenu, 'Text', 'Set # Bins (All Groups)', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.setHistogramNumBins([], obj.ui.HistAxes));
            uimenu(histMenu, 'Text', 'Set Axes Width', ...
                'MenuSelectedFc', @(varargin) obj.setHistogramAxesWidth());
            uimenu(menu, 'Text', 'Refresh', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.refresh());
        end
        
        %------------------------------------------------------------------
        % UI
        %------------------------------------------------------------------
        
        % Create the toolbar panel that sits above the plots
        % (but still within the main panel)
        function panel = createToolbarPanel(obj)
            panel = uipanel(obj.ui.MainPanel, ...
                'BorderType', 'none', ...
                'Units', 'pixels', ...
                'Position', [0 0 2 24], ...
                'AutoResizeChildren', 0);
            x = 1;
            
%             [iptIconsDir, matlabIconsDir] = ipticondir();
            
            obj.ui.VisibleSweepsEdit = uicontrol(panel, ...
                'Style', 'edit', ...
                'String', '0 of 0', ...
                'Tooltip', 'Visible Sweep(s)', ...
                'Units', 'pixels', ...
                'Position', [x 1 100 22], ...
                'Callback', @(varargin) obj.setVisibleSweeps());
            x = x + 100;
            obj.ui.PrevSweepButton = uicontrol(panel, ...
                'Style', 'pushbutton', ...
                'String', char(hex2dec('25bd')), ...
                'Tooltip', 'Previous Sweep', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 11], ...
                'Callback', @(varargin) obj.prevSweep());
            obj.ui.NextSweepButton = uicontrol(panel, ...
                'Style', 'pushbutton', ...
                'String', char(hex2dec('25b3')), ...
                'Tooltip', 'Next Sweep', ...
                'Units', 'pixels', ...
                'Position', [x 12 22 11], ...
                'Callback', @(varargin) obj.nextSweep());
            x = x + 32;
            
            obj.ui.PageLeftButton = uicontrol(panel, ...
                'Style', 'pushbutton', ...
                'String', '<', ...
                'Tooltip', 'Page Left', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(varargin) obj.pageLeft());
            x = x + 23;
            obj.ui.PageRightButton = uicontrol(panel, ...
                'Style', 'pushbutton', ...
                'String', '>', ...
                'Tooltip', 'Page Right', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(varargin) obj.pageRight());
            x = x + 33;
            
            obj.ui.XROIsButton = uicontrol(panel, ...
                'Style', 'togglebutton', ...
                'String', '[ ]', ...
                'Tooltip', 'Draw/Edit X-Axis ROIs', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(varargin) obj.XROIsButtonChanged());
            if obj.hasXROIMgr()
                h = addlistener(obj.XROIMgr, 'NumberOfROIsChanged', ...
                    @(varargin) obj.updateXROIsButton());
                setappdata(obj.ui.XROIsButton, 'NumberOfROIsChangedListener', h);
            end
            x = x + 32;
            
            warning('off', 'MATLAB:structOnObject');
            [tb, btns] = axtoolbar(obj.ui.HistAxes(1), {'brush', 'pan', 'zoomin', 'zoomout', 'restoreview'});
            
            try
                brushBtn = btns(end);
                I = brushBtn.Image;
                I(I == 64) = nan;
                I = I ./ 64;
                rgb = cat(3, I, I, I);
            catch
                rgb = [];
            end
            obj.ui.BrushButton = uicontrol(panel, ...
                'Style', 'togglebutton', ...
                'CData', rgb, ...
                'Tooltip', 'Brush', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(s,e) obj.setBrushMode(s.Value));
            if isempty(rgb)
                obj.ui.BrushButton.String = 'B';
            end
            x = x + 23;
            
            try
                panBtn = findall(obj.ui.TsAxes(1).Toolbar, 'Tooltip', 'Pan');
                I = panBtn.Image;
                I(I == 64) = nan;
                I = I ./ 64;
                rgb = cat(3, I, I, I);
            catch
                rgb = [];
            end
            obj.ui.PanButton = uicontrol(panel, ...
                'Style', 'togglebutton', ...
                'CData', rgb, ...
                'Tooltip', 'Pan', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(s,e) obj.setPanMode(s.Value));
            if isempty(rgb)
                obj.ui.PanButton.String = 'P';
            end
            x = x + 23;
            
            try
                zoomInBtn = findall(obj.ui.TsAxes(1).Toolbar, 'Tooltip', 'Zoom In');
                I = zoomInBtn.Image;
                I(I == 64) = nan;
                I = I ./ 64;
                rgb = cat(3, I, I, I);
            catch
                rgb = [];
            end
            obj.ui.ZoomButton = uicontrol(panel, ...
                'Style', 'togglebutton', ...
                'CData', rgb, ...
                'Tooltip', 'Zoom', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(s,e) obj.setZoomMode(s.Value));
            if isempty(rgb)
                obj.ui.ZoomButton.String = 'Z';
            end
            x = x + 23;
            
            try
                restoreViewBtn = findall(obj.ui.TsAxes(1).Toolbar, 'Tooltip', 'Restore View');
                I = restoreViewBtn.Image;
                I(I == 64) = nan;
                I = I ./ 64;
                rgb = cat(3, I, I, I);
            catch
                rgb = [];
            end
            obj.ui.RestoreViewButton = uicontrol(panel, ...
                'Style', 'pushbutton', ...
                'CData', rgb, ...
                'Tooltip', 'Restore View', ...
                'Units', 'pixels', ...
                'Position', [x 1 22 22], ...
                'Callback', @(varargin) obj.autoscaleXY());
            if isempty(rgb)
                obj.ui.RestoreViewButton.String = 'A';
            end
            x = x + 23;
            
            panel.Position(3) = x;
        end
        function refreshToolbarPanel(obj)
            x = 1;
            
            tf = obj.numGroups() && max(obj.numTsPerGroup()) > 1;
            obj.ui.VisibleSweepsEdit.Visible = tf;
            obj.ui.PrevSweepButton.Visible = tf;
            obj.ui.NextSweepButton.Visible = tf;
            if tf
                obj.ui.VisibleSweepsEdit.Position(1) = x;
                x = x + obj.ui.VisibleSweepsEdit.Position(3);
                obj.ui.PrevSweepButton.Position(1) = x;
                obj.ui.NextSweepButton.Position(1) = x;
                x = x + obj.ui.NextSweepButton.Position(3) + 10;
            end
            
            obj.ui.PageLeftButton.Position(1) = x;
            x = x + obj.ui.PageLeftButton.Position(3) + 1;
            obj.ui.PageRightButton.Position(1) = x;
            x = x + obj.ui.PageRightButton.Position(3) + 10;
            
            if obj.hasXROIMgr()
                obj.ui.XROIsButton.Position(1) = x;
                x = x + obj.ui.XROIsButton.Position(3) + 10;
                obj.ui.XROIsButton.Visible = true;
            else
                obj.ui.XROIsButton.Visible = false;
            end
            
            obj.ui.BrushButton.Visible = false;
%             obj.ui.BrushButton.Position(1) = x;
%             x = x + obj.ui.BrushButton.Position(3) + 1;
            obj.ui.PanButton.Position(1) = x;
            x = x + obj.ui.PanButton.Position(3) + 1;
            obj.ui.ZoomButton.Position(1) = x;
            x = x + obj.ui.ZoomButton.Position(3) + 1;
            obj.ui.RestoreViewButton.Position(1) = x;
            x = x + obj.ui.RestoreViewButton.Position(3) + 1;
            
            obj.ui.ToolbarPanel.Position(3) = x;
        end
        
        % Open a new app window
        function newWindow(obj)
            TimeSeriesAnalyzer();
        end
        
        % Reparent the main panel which contains all graphics objects.
        % Also update figure listeners to new outer figure container.
        function set.Parent(obj, parent)
            obj.ui.MainPanel.Parent = parent;
            obj.ui.Figure = ancestor(parent, {'Figure', 'matlab.ui.Figure'}, 'toplevel');
            h = addlistener(obj.ui.Figure, 'WindowMousePress', @obj.windowMousePress);
            setappdata(obj.ui.MainPanel, 'WindowMousePressListener', h);
%             h = addlistener(obj.ui.Figure, 'WindowKeyPress', @obj.windowKeyPress);
%             setappdata(obj.ui.MainPanel, 'WindowKeyPressListener', h);
%             h = addlistener(obj.ui.Figure, 'WindowKeyRelease', @obj.windowKeyRelease);
%             setappdata(obj.ui.MainPanel, 'WindowKeyReleaseListener', h);
            obj.resize();
        end
        function parent = get.Parent(obj)
            parent = obj.ui.MainPanel.Parent;
        end
        
        % Update all graphics objects.
        % Creates/destroys plot axes as needed.
        % Call this whenever you reset the data or change the groups.
        function refresh(obj)
            % add/remove plot axes to match visible groups
            n = obj.numVisibleGroups();
            while numel(obj.ui.TsAxes) < n
                [obj.ui.TsAxes(end+1), obj.ui.HistAxes(end+1)] = obj.createAxes();
            end
            i = max(n+1, 2); % don't delete the first axes
            delete(obj.ui.TsAxes(i:end));
            delete(obj.ui.HistAxes(i:end));
            obj.ui.TsAxes(i:end) = [];
            obj.ui.HistAxes(i:end) = [];
            % link axes
            for i = 1:numel(obj.ui.TsAxes)
                linkaxes([obj.ui.TsAxes(i) obj.ui.HistAxes(i)], 'y');
            end
            if numel(obj.ui.TsAxes) > 1
                linkaxes(obj.ui.TsAxes, 'x');
            end
            % update x-axis ROIs
            if obj.hasXROIMgr()
                obj.XROIMgr.Axes = obj.ui.TsAxes; % will notify to reset ROI context menus
                if ~isempty(obj.XROIMgr.ROIs)
                    [obj.XROIMgr.ROIs.Visible] = deal(obj.ui.XROIsButton.Value);
                end
            end
            % refresh toolbar controls
            obj.refreshToolbarPanel();
            % reposition components
            obj.resize();
            % update visible traces and replot
            obj.setVisibleSweeps(); % calls replot()
        end
        
        % Reposition all graphics objects to fit in the main panel.
        function resize(obj)
            % reposition all graphics objects within main panel
            bbox = getpixelposition(obj.ui.MainPanel);
            margin = 5;
            left = margin;
            right = bbox(3) - margin;
            bottom = margin;
            top = bbox(4) - margin;
            
            % top control panel
            if isfield(obj.ui, 'ControlPanel') && isvalid(obj.ui.ControlPanel) && obj.ui.ControlPanel.Visible
                y = top - obj.ui.ControlPanel.Position(4);
                obj.ui.ControlPanel.Position(1:3) = [left, y, right-left];
                top = y - margin;
            end
            
            % top toolbar
            if isfield(obj.ui, 'ToolbarPanel') && isvalid(obj.ui.ToolbarPanel) && obj.ui.ToolbarPanel.Visible
                x = (left + right) / 2 - obj.ui.ToolbarPanel.Position(3) / 2;
                y = top - obj.ui.ToolbarPanel.Position(4);
                obj.ui.ToolbarPanel.Position(1:2) = [x, y];
                top = y - margin;
            end
            
            % plot axes
            tx = left + 72;
            y = bottom + 22;
            if any(horzcat(obj.ui.HistAxes.Visible))
                hx = right - obj.ui.HistAxes(1).Position(3);
                tw = hx - tx - 3;
                hw = obj.ui.HistAxes(1).Position(3);
            else
                hx = right;
                hw = obj.ui.HistAxes(1).Position(3);
                tw = right - tx - 5;
            end
            nax = numel(obj.ui.TsAxes);
            h = (top - y) / nax;
            for i = 1:nax
                y = top - i * h;
                obj.ui.TsAxes(i).Position = [tx, y+22, tw, h-22];
                obj.ui.HistAxes(i).Position = [hx, y+22, hw, h-22];
            end
        end
        
        % Replot all visible plots and remove nonvisible plots.
        % Applies consitent colormap to plots based on name.
        function replot(obj)
            % visible groups/ts/names
            visGroups = obj.visibleGroups();
            visTsPerGroup = obj.visibleTsPerGroup();
            names = obj.names();
            visNames = obj.visibleNames();
            % colormap (index will be based on plot name)
            cmap = lines();
            if ~isempty(names) && names(1) == "raw"
                % raw data shown in gray
                cmap = [0.5 0.5 0.5; cmap];
            end
            % for each ts axes...
            for i = 1:numel(obj.ui.TsAxes)
                ax = obj.ui.TsAxes(i);
                hold(ax, 'on');
                % grab plot handles for axes
                hplots = getappdata(ax, 'TsPlots');
                if isempty(hplots)
                    hplots = gobjects(0);
                else
                    % remove any invalid plot handles
                    hplots(~isvalid(hplots)) = [];
                end
                try
                    group = visGroups(i);
                    tsi = visTsPerGroup{group};
                    nts = numel(tsi);
                catch
                    group = nan;
                    tsi =[];
                    nts = 0;
                end
                hidx = 1;
                fitind = [];
                % time series
                for j = 1:nts
                    ts = obj.Data(tsi(j));
                    for k = 1:numel(visNames)
                    	name = visNames(k);
                        y = ts.("y" + name);
                        if isempty(y)
                            continue
                        end
                        if ~isnumeric(y) || obj.isideal(y)
                            ls = '-';
                            lw = 1.5;
                        else
                            ls = get(groot, 'defaultLineLineStyle');
                            lw = get(groot, 'defaultLineLineWidth');
                        end
                        [x,y] = obj.getXY(ts, name);
                        if hidx > numel(hplots)
                            h = obj.createPlot(ax, x, y);
                        else
                            h = obj.updatePlot(hplots(hidx), x, y);
                        end
                        if name == "fit"
                            fitind = [fitind hidx];
                            h.Color = 'r';
                            h.LineStyle = '-';
                            h.LineWidth = 1.5;
                            h.ContextMenu = obj.fitContextMenu(h);
                            h.HitTest = 'on';
                            h.PickableParts = 'visible';
                        else
                            colorIndex = find(names == name, 1);
                            h.Color = cmap(colorIndex,:);
                            h.LineStyle = ls;
                            h.LineWidth = lw;
                            delete(h.ContextMenu);
                            h.HitTest = 'off';
                            h.PickableParts = 'none';
                        end
                        setappdata(h, 'TsIndex', tsi(j));
                        setappdata(h, 'TsName', name);
                        hplots(hidx) = h;
                        hidx = hidx + 1;
                    end
                end
                % fits on top
                if ~isempty(fitind)
                    uistack(hplots(fitind), 'top');
                end
                % get rid of extra plot handles
                delete(hplots(hidx:end));
                hplots(hidx:end) = [];
                % store plot handles in axes
                setappdata(ax, 'TsPlots', hplots);
                % axes labels from first visible trace in group
                if i ~= numel(obj.ui.TsAxes)
                    ax.XLabel.String = '';
                else
                    try
                        ax.XLabel.String = obj.Data(tsi(1)).xlabel;
                    catch
                        ax.XLabel.String = '';
                    end
                end
                try
                    ax.YLabel.String = obj.Data(tsi(1)).ylabel;
                    if isempty(ax.YLabel.String)
                        ax.YLabel.String = ['Group ' num2str(group)];
                    end
                catch
                    ax.YLabel.String = ['Group ' num2str(group)];
                end
            end
            % ROI context menus
            obj.refreshXROIsContextMenus();
            % histograms
            obj.replotHistograms();
        end
        function replotHistograms(obj, ax)
            if ~exist('ax', 'var') || isempty(ax)
                ax = obj.ui.HistAxes;
            end
            for i = 1:numel(ax)
                hax = ax(i);
                if ~hax.Visible
                    continue
                end
                hhists = getappdata(hax, 'TsPlots');
                if isempty(hhists)
                    hhists = gobjects(0);
                else
                    hhists(~isvalid(hhists)) = [];
                end
                haxIdx = find(obj.ui.HistAxes == hax, 1);
                tax = obj.ui.TsAxes(haxIdx);
                hplots = getappdata(tax, 'TsPlots');
                if isempty(hplots)
                    hplots = gobjects(0);
                else
                    hplots(~isvalid(hplots)) = [];
                end
                hidx = 1;
                for j = 1:numel(hplots)
                    htsj = hplots(j);
                    name = getappdata(htsj, 'TsName');
                    if ismember(name, ["raw" "data"])
                        [~,y] = obj.getVisibleXYFromPlot(htsj);
                        % plot histogram
                        if hidx > numel(hhists)
                            h = obj.createHistogram(hax, y);
                        else
                            h = obj.updateHistogram(hhists(hidx), y);
                        end
                        h.FaceColor = htsj.Color;
                        setappdata(h, 'TsIndex', getappdata(htsj, 'TsIndex'));
                        setappdata(h, 'TsName', name);
                        hhists(hidx) = h;
                        hidx = hidx + 1;
                    end
                end
                % get rid of extra plot handles
                delete(hhists(hidx:end));
                hhists(hidx:end) = [];
                % store plot handles in axes
                setappdata(hax, 'TsPlots', hhists);
                % same y-axis as trace
                hax.YLim = tax.YLim;
            end
        end
        
        % Show message box describing how to adjust default line style
        function showLineStyleMessageBox(obj)
            msgbox([ ...
                "Plots use the default line style/width." ...
                "" ...
                "e.g. to change these:" ...
                ">> set( groot, 'defaultLineLineStyle', '--' );" ...
                ">> set( groot, 'defaultLineLineWidth', 2.0 );" ...
                "" ...
                "Then refresh the plots: 'View->Refresh'" ...
                ], 'Line Style/Width', 'help');
        end
        
        %------------------------------------------------------------------
        % X-Axis Ranges ROI Manager
        %------------------------------------------------------------------
        
        % Query validity of xrangesMgr
        function tf = hasXROIMgr(obj)
            tf = ~isempty(obj.XROIMgr) && isvalid(obj.XROIMgr);
        end
        
        % Get current number of xrange ROIs per plot.
        % Will return 0 for an invalid xrangesMgr.
        function n = numXROIs(obj)
            if obj.hasXROIMgr()
                n = obj.XROIMgr.numROIs();
            else
                n = 0;
            end
        end
        
        % Update button to indicate the current number of ROIs.
        function updateXROIsButton(obj)
            n = obj.numXROIs();
            if n
                obj.ui.XROIsButton.String = ['[' num2str(n) ']'];
            else
                obj.ui.XROIsButton.String = '[ ]';
            end
        end
        
        % Callback for xranges button toggled.
        function XROIsButtonChanged(obj)
            if obj.numXROIs()
                [obj.XROIMgr.ROIs.Visible] = deal(obj.ui.XROIsButton.Value);
            end
            obj.updateXROIsButton();
        end
        
        % Refresh the ROI context menus.
        function refreshXROIsContextMenus(obj)
            if ~obj.hasXROIMgr()
                return
            end
            for i = 1:size(obj.XROIMgr.ROIs, 1)
                ax = obj.ui.TsAxes(i);
                for j = 1:size(obj.XROIMgr.ROIs, 2)
                    roi = obj.XROIMgr.ROIs(i,j);
                    roi.UIContextMenu = obj.tsAxesContextMenu(ax, roi);
                end
            end
        end
        
        % Merge overlapping ROIs.
        function mergeOverlappingXROIs(obj)
            obj.XROIMgr.mergeOverlappingROIs();
            obj.refreshXROIsContextMenus();
        end
        
        %------------------------------------------------------------------
        % Measurement
        %------------------------------------------------------------------
        
        % Measurement for plot objects in axes
        function measure(obj, ax, method)
            hplots = gobjects(0);
            for i = 1:numel(ax)
                h = getappdata(ax(i), 'TsPlots');
                hplots = [hplots h];
            end
            hplots(~isvalid(hplots)) = [];
            if isempty(hplots)
                return
            end
            nplots = numel(hplots);
            nrois = max(1, obj.numXROIs());
            if ~obj.ui.XROIsButton.Value
                nrois = 1;
            end
            xy = nan(nrois, 2, nplots);
            htmp = gobjects(0);
            for i = 1:nplots
                % get plot data
                [xdata,ydata] = obj.getXYFromPlot(hplots(i));
                
                % get data indices included in measurement
                if obj.numXROIs() && obj.ui.XROIsButton.Value
                    idxPerROI = obj.XROIMgr.xindices(xdata);
                else
                    idxPerROI{1} = 1:length(ydata);
                end
                for j = 1:nrois
                    idx = idxPerROI{j};
                    idx(isnan(ydata(idx))) = [];
                    if isempty(idx)
                        continue
                    end
                    x = xdata(idx);
                    y = ydata(idx);
                    if method == "mean"
                        xy(j,:,i) = [mean(x) mean(y)];
                    elseif method == "std"
                        xy(j,:,i) = [mean(x) std(y)];
                    elseif method == "var"
                        xy(j,:,i) = [mean(x) var(y)];
                    elseif method == "min"
                        [ymin, imin] = min(y);
                        xmin = x(imin);
                        xy(j,:,i) = [xmin ymin];
                    elseif method == "max"
                        [ymax, imax] = max(y);
                        xmax = x(imax);
                        xy(j,:,i) = [xmax ymax];
                    elseif method == "absmax"
                        [yabsmax, iabsmax] = max(abs(y));
                        xabsmax = x(iabsmax);
                        xy(j,:,i) = [xabsmax yabsmax];
                    elseif method == "peak"
                        [~, ipeak] = max(abs(y));
                        ypeak = y(ipeak);
                        xpeak = x(ipeak);
                        xy(j,:,i) = [xpeak ypeak];
                    end
                end
                % brief graphical presentation of measurement
                if any(~isnan(xy(:,:,i)))
                    ax = gca;
                    hold(ax, 'on');
                    htmp(end+1) = plot(ax, xy(:,1,i), xy(:,2,i), 'ro', ...
                        'MarkerSize', 8, 'markerFaceColor', 'r');
                end
            end
            % store measurement organized as table
            tsk = [];
            ylabels = string.empty;
            groups = [];
            sweeps = [];
            names = string.empty;
            rois = [];
            x = [];
            y = [];
            for i = 1:numel(hplots)
                tsi = getappdata(hplots(i), 'TsIndex');
                name = getappdata(hplots(i), 'TsName');
                ts = obj.Data(tsi);
                gtsi = obj.tsPerGroup();
                gtsi = gtsi{ts.group};
                sweepi = find(gtsi == tsi, 1);
                for j = 1:nrois
                    tsk = [tsk; tsi];
                    ylabels = [ylabels; string(ts.ylabel)];
                    groups = [groups; ts.group];
                    sweeps = [sweeps; sweepi];
                    names = [names; name];
                    rois = [rois; j];
                    x = [x; xy(j,1,i)];
                    y = [y; xy(j,2,i)];
                end
            end
            obj.Measurement = table;
            if numel(obj.Data) > 1
                obj.Measurement.Index = tsk;
            end
            if any(arrayfun(@(ylabel) strlength(ylabel) > 0, ylabels))
                obj.Measurement.YLabel = ylabels;
            end
            if obj.numGroups() > 1
                obj.Measurement.Group = groups;
                obj.Measurement.Sweep = sweeps;
            end
            obj.Measurement.Name = names;
            if nrois > 1
                obj.Measurement.ROI = rois;
            end
            xtitle = "mean(x)";
            if ismember(string(method), ["min", "max", "absmax", "peak"])
                xtitle = "x";
            end
            ytitle = string(method) + "(y)";
            obj.Measurement.(xtitle) = x;
            obj.Measurement.(ytitle) = y;
            % print measurement in command window
            disp(obj.Measurement);
            % brief graphical presentation of measurement
            drawnow;
            pause(2);
            delete(htmp);
        end
        
        % menu
        function menu = measureMenu(obj, parent, ax)
            menu = uimenu(parent, 'Text', 'Measure');
            methods = {'mean', 'std', 'var', 'min', 'max', 'absmax', 'peak'};
            for i = 1:numel(methods)
                uimenu(menu, 'Text', methods{i}, ...
                    'MenuSelectedFc', @(varargin) obj.measure(ax, methods{i}));
            end
        end
        
        %------------------------------------------------------------------
        % Curve Fitting
        %------------------------------------------------------------------
        
        % Curve fit plot objects
        function curveFit(obj, hplotsOrAxes, methodOrExpr, params, lowerbound, upperbound, xfitlim)
            if class(hplotsOrAxes) == "matlab.graphics.axis.Axes"
                ax = hplotsOrAxes;
                hplots = gobjects(0);
                for i = 1:numel(ax)
                    h = getappdata(ax(i), 'TsPlots');
                    hplots = [hplots h];
                end
                hplots(~isvalid(hplots)) = [];
                tf = arrayfun(@(h) getappdata(h, 'TsName') == "data", hplots);
                hplots = hplots(tf);
            else
                hplots = hplotsOrAxes;
                hplots(~isvalid(hplots)) = [];
            end
            if isempty(hplots)
                return
            end
            
            % input method or expression if not given
            if ~exist('methodOrExpr', 'var') || isempty(methodOrExpr) || methodOrExpr == "custom"
                answer = inputdlg({'Expression '}, 'Curve Fit', 1, {'a*exp(-x/b)+c'});
                if isempty(answer) || isempty(answer{1})
                    return
                end
                methodOrExpr = answer{1};
            end
            
            % fit each plot in h
            for i = 1:numel(hplots)
                % get plot data
                [xdata,ydata] = obj.getXYFromPlot(hplots(i));
                
                xfit = [];
                yfit = [];
                
                % get data indices included in fit
                if obj.numXROIs() && obj.ui.XROIsButton.Value
                    [~,idx] = obj.XROIMgr.xindices(xdata);
                else
                    idx = 1:length(ydata);
                end
                idx(isnan(ydata(idx))) = [];
                if isempty(idx)
                    continue
                end
                x = reshape(xdata(idx), [], 1);
                y = reshape(ydata(idx), [], 1);
                
                % limit fit range?
                if exist('xfitlim', 'var') && ~isempty(xfitlim)
                    xfit = xdata;
                    xfit(xfit < xfitlim(1)) = [];
                    xfit(xfit > xfitlim(2)) = [];
                end
                
                % fit based on method
                if methodOrExpr == "mean"
                    yfit = cfit(fittype('poly1'), 0, mean(y));
                elseif methodOrExpr == "line"
                    p = polyfit(x, y, 1);
                    yfit = cfit(fittype('poly1'), p(1), p(2));
                elseif methodOrExpr == "polynomial"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'Polynomial Order '}, 'Polynomial Fit', 1);
                        if isempty(answer)
                            return
                        end
                        params = str2num(answer{1});
                    end
                    order = params;
                    p = polyfit(x, y, order);
                    lm = {};
                    np = numel(p);
                    for j = 1:np
                        if j == np
                            lm{j} = '1';
                        else
                            lm{j} = ['x.^' num2str(np - j)];
                        end
                    end
                    p = mat2cell(p, ones(1,size(p,1)), ones(1,size(p,2)));
    %                 yfit = polyval(p, x);
                    yfit = cfit(fittype(lm), p{:});
                elseif methodOrExpr == "spline"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'# of Spline Segments '}, 'Spline Fit', 1);
                        if isempty(answer)
                            return
                        end
                        params = str2num(answer{1});
                    end
                    try
                        numSegments = params;
                        pp = splinefit(x, y, numSegments);
    %                     yfit = ppval(pp, x);
                        yfit = pp;
                    catch err
                        disp(err);
                        msgbox("!!! Requires package 'splinefit'. See Add-On Explorer.", ...
                            'splinefit');
                        return
                    end
                else % e.g. 'a * exp(-x/b)'
                    try
                        ft = fittype(methodOrExpr);
                        options = fitoptions(ft);
                        if ~exist('params', 'var') || isempty(params)
                            coeffs = coeffnames(ft);
                            answer = inputdlg(coeffs, 'Starting Params', 1);
                            if isempty(answer)
                                return
                            end
                            params = zeros(1, numel(coeffs));
                            for j = 1:numel(params)
                                val = str2num(answer{j});
                                if ~isempty(val)
                                    params(j) = val;
                                end
                            end
                        end
                        options.StartPoint = params;
                        if exist('lowerbound', 'var') && ~isempty(lowerbound)
                            options.Lower = lowerbound;
                        end
                        if exist('upperbound', 'var') && ~isempty(upperbound)
                            options.Upper = upperbound;
                        end
                        yfit = fit(x, y, ft, options);
                    catch err
                        disp(err);
                        return
                    end
                end
                
                % add fit to time series data
                if ~isempty(yfit)
                    if ~isfield(obj.Data, 'yfit')
                        [obj.Data.yfit] = deal([]);
                    end
                    tsi = getappdata(hplots(i), 'TsIndex');
                    obj.Data(tsi).xfit = xfit;
                    obj.Data(tsi).yfit = yfit;
                    if ~isnumeric(yfit)
                        disp(yfit); % print fit result
                    end
                end
            end
                
            % show fits
            visNames = obj.visibleNames();
            if ~ismember("fit", visNames)
                visNames(end+1) = "fit";
                obj.setVisibleNames(visNames);
            end
            obj.replot();
        end
        
        % Delete fit for time series indices tsind
        function deleteFit(obj, hfit)
            tsi = getappdata(hfit, 'TsIndex');
            if isfield(obj.Data, 'xfit')
                obj.Data(tsi).xfit = [];
            end
            if isfield(obj.Data, 'yfit')
                obj.Data(tsi).yfit = [];
            end
            obj.replot();
        end
        
        % Delete fit for time series indices tsind
        function printFitParams(obj, hfit)
            tsi = getappdata(hfit, 'TsIndex');
            if isfield(obj.Data, 'yfit') && ~isnumeric(obj.Data(tsi).yfit)
                obj.Data(tsi).yfit
            end
        end
        
        % Subtract yfit from ydata
        function subtractFit(obj, hfit)
            % find plot handle for ydata associated with hfit
            ax = hfit.Parent;
            tsi = getappdata(hfit, 'TsIndex');
            hdata = obj.getPlot(ax, tsi, "data");
            if isempty(hdata) || ~isvalid(hdata)
                return
            end
            % subtract hfit from hydata
            [xdata,ydata] = obj.getXYFromPlot(hdata);
            [xfit,yfit] = obj.getXYFromPlot(hfit);
            if isequal(xdata, xfit)
                ydata = ydata - yfit;
            else
                try
                    first = find(xdata >= xfit(1), 1);
                    last = length(xdata) + 1 - find(xdata(end:-1:1) <= xfit(end), 1);
                    idx = first:last;
                    ydata(idx) = ydata(idx) - interp1(xfit, yfit, xdata(idx));
                catch
                    return
                end
            end
            ylim = ax.YLim;
            ylimmode = ax.YLimMode;
            ax.YLim = [-max(yfit - ylim(1)), max(ylim(2) - yfit)];
            yfit(:) = 0;
            obj.updatePlot(hdata, xdata, ydata);
            obj.updatePlot(hfit, xfit, yfit);
            % update histogram
            i = find(obj.ui.TsAxes == ax, 1);
            hax = obj.ui.HistAxes(i);
            hhist = obj.getPlot(hax, tsi, "data");
            if ~isempty(hhist) && isvalid(hhist)
                obj.updateHistogram(hhist, ydata);
            end
            % ask if we should keep result?
            if questdlg('Subtract fit?', 'Subtract Fit') == "Yes"
                % update time series data
                obj.Data(tsi).ydata = ydata;
                obj.Data(tsi).yfit = [];
                if isfield(obj.Data, 'xfit')
                    obj.Data(tsi).xfit = [];
                end
            else
                ax.YLim = ylim;
                ax.YLimMode = ylimmode;
            end
            obj.replot();
        end
        
        % Normalize data to fit for time series indices tsind
        function normalizeToFit(obj, hfit)
            % find plot handle for ydata associated with hfit
            ax = hfit.Parent;
            tsi = getappdata(hfit, 'TsIndex');
            hdata = obj.getPlot(ax, tsi, "data");
            if isempty(hdata) || ~isvalid(hdata)
                return
            end
            % normalize hdata to hfit
            [xdata,ydata] = obj.getXYFromPlot(hdata);
            [xfit,yfit] = obj.getXYFromPlot(hfit);
            if isequal(xdata, xfit)
                ydata = ydata ./ yfit;
            else
                try
                    first = find(xdata >= xfit(1), 1);
                    last = length(xdata) + 1 - find(xdata(end:-1:1) <= xfit(end), 1);
                    idx = first:last;
                    ydata(idx) = ydata(idx) ./ interp1(xfit, yfit, xdata(idx));
                catch
                    return
                end
            end
            ylim = ax.YLim;
            ylimmode = ax.YLimMode;
            ax.YLim = [1 - max((yfit - ylim(1)) ./ yfit), 1 + max((ylim(2) - yfit) ./ yfit)];
            yfit(:) = 1;
            obj.updatePlot(hdata, xdata, ydata);
            obj.updatePlot(hfit, xfit, yfit);
            % update histogram
            i = find(obj.ui.TsAxes == ax, 1);
            hax = obj.ui.HistAxes(i);
            hhist = obj.getPlot(hax, tsi, "data");
            if ~isempty(hhist) && isvalid(hhist)
                obj.updateHistogram(hhist, ydata);
            end
            % ask if we should keep result?
            if questdlg('Normalize to fit?', 'Normalize to Fit') == "Yes"
                % update time series data
                obj.Data(tsi).ydata = ydata;
                obj.Data(tsi).yfit = [];
                if isfield(obj.Data, 'xfit')
                    obj.Data(tsi).xfit = [];
                end
            else
                ax.YLim = ylim;
                ax.YLimMode = ylimmode;
            end
            obj.replot();
        end
        
        % menu
        function menu = fitMenu(obj, parent, ax, roi)
            menu = uimenu(parent, 'Text', 'Fit');
            isroi = exist('roi', 'var') && isvalid(roi);
            methods = {'mean', 'line', 'polynomial', 'spline', 'custom'};
            if isroi
                xfitlim = cumsum(roi.Position([1 3]));
                submenu = uimenu(menu, 'Text', 'Fit in ROI');
                for i = 1:numel(methods)
                    uimenu(submenu, 'Text', methods{i}, ...
                        'MenuSelectedFc', @(varargin) obj.curveFit( ...
                        ax, methods{i}, [], [], [], xfitlim));
                end
            end
            for i = 1:numel(methods)
                uimenu(menu, 'Text', methods{i}, ...
                    'Separator', (i == 1) && isroi, ...
                    'MenuSelectedFc', @(varargin) obj.curveFit(ax, methods{i}));
            end
        end
        
        % Fit plot object context menu
        function menu = fitContextMenu(obj, hfit)
            menu = uicontextmenu(obj.ui.Figure);
            uimenu(menu, 'Text', 'Delete Fit', ...
                'MenuSelectedFc', @(varargin) obj.deleteFit(hfit));
            uimenu(menu, 'Text', 'Print Fit Parameters', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.printFitParams(hfit));
            uimenu(menu, 'Text', 'Subtract Fit', ...
                'Separator', true, ...
                'MenuSelectedFc', @(varargin) obj.subtractFit(hfit));
            uimenu(menu, 'Text', 'Normalize to Fit', ...
                'MenuSelectedFc', @(varargin) obj.normalizeToFit(hfit));
        end
        
        %------------------------------------------------------------------
        % Data Operations
        %------------------------------------------------------------------
        
        function operation(obj, hplotsOrAxes, method, params)
            if class(hplotsOrAxes) == "matlab.graphics.axis.Axes"
                ax = hplotsOrAxes;
                hplots = gobjects(0);
                for i = 1:numel(ax)
                    h = getappdata(ax(i), 'TsPlots');
                    hplots = [hplots h];
                end
                hplots(~isvalid(hplots)) = [];
                tf = arrayfun(@(h) getappdata(h, 'TsName') ~= "raw", hplots);
                hplots = hplots(tf);
            else
                hplots = hplotsOrAxes;
                hplots(~isvalid(hplots)) = [];
            end
            if isempty(hplots)
                return
            end
            
            % operation on each plot in h
            for i = 1:numel(hplots)
                % get plot data
                [xdata,ydata] = obj.getXYFromPlot(hplots(i));
                
                % get data indices included in operation
                if obj.numXROIs() && obj.ui.XROIsButton.Value
                    [~,idx] = obj.XROIMgr.xindices(xdata);
                else
                    idx = 1:length(ydata);
                end
                idx(isnan(ydata(idx))) = [];
                if isempty(idx)
                    continue
                end
                x = reshape(xdata(idx), [], 1);
                y = reshape(ydata(idx), [], 1);
                
                % operations
                if method == "mask"
                    y(:) = nan;
                elseif method == "zero"
                    y(:) = 0;
                elseif method == "interpolate"
                    notidx = setdiff(1:length(ydata), idx);
                    if exist('params', 'var')
                        try
                            y = interp1(xdata(notidx), ydata(notidx), x, params);
                        catch
                            y = interp1(xdata(notidx), ydata(notidx), x);
                        end
                    else
                        y = interp1(xdata(notidx), ydata(notidx), x);
                    end
                elseif method == "set"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'Value '}, 'Set ', 1);
                        if isempty(answer) || isempty(answer{1})
                            return
                        end
                        params = str2num(answer{1});
                    end
                    y(:) = params;
                elseif method == "add" || method == "+"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'Value '}, 'Add ', 1);
                        if isempty(answer) || isempty(answer{1})
                            return
                        end
                        params = str2num(answer{1});
                    end
                    y = y + params;
                elseif method == "subtract" || method == "-"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'Value '}, 'Subtract ', 1);
                        if isempty(answer) || isempty(answer{1})
                            return
                        end
                        params = str2num(answer{1});
                    end
                    y = y - params;
                elseif method == "multiply" || method == "*"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'Value '}, 'Multiply ', 1);
                        if isempty(answer) || isempty(answer{1})
                            return
                        end
                        params = str2num(answer{1});
                    end
                    y = y .* params;
                elseif method == "divide" || method == "/"
                    if ~exist('params', 'var') || isempty(params)
                        answer = inputdlg({'Value '}, 'Divide ', 1);
                        if isempty(answer) || isempty(answer{1})
                            return
                        end
                        params = str2num(answer{1});
                    end
                    y = y ./ params;
                end
                % show operation in plot
                ynew = ydata;
                ynew(idx) = y;
                hplots(i) = obj.updatePlot(hplots(i), xdata, ynew);
            end
            % update histograms
            ... % TODO
            % ask if we should keep the result?
            if questdlg('Apply operation?', 'Apply Operation?') == "Yes"
                for i = 1:numel(hplots)
                    obj.updateDataFromPlot(hplots(i));
                end
            end
            obj.replot();
        end
        
        % menu
        function operationMenu(obj, menu, ax)
            ops = {'Mask', 'Zero', 'Interpolate', 'Set'};
            for i = 1:numel(ops)
                uimenu(menu, 'Text', ops{i}, ...
                    'Separator', i == 1, ...
                    'MenuSelectedFc', @(varargin) obj.operation(ax, lower(ops{i}), []));
            end
            mathMenu = uimenu(menu, 'Text', 'Math');
            mathops = {'+', '-', '*', '/'};
            for i = 1:numel(mathops)
                uimenu(mathMenu, 'Text', mathops{i}, ...
                    'MenuSelectedFc', @(varargin) obj.operation(ax, mathops{i}, []));
            end
        end
        
        %------------------------------------------------------------------
        % Idealization for piecewise continuous data
        %------------------------------------------------------------------
        
        function tf = isideal(obj, y)
            % If y does not change for more than half the sample steps,
            % assume it's an idealization of some sort.
            y(isnan(y)) = [];
            tf = sum(diff(y) == 0) > length(y) / 2;
        end
        
        function idealize(obj, tsIndices, method, params)
            if ~exist('params', 'var')
                if method == "DISC"
                    params = initDISC();
                    labels = fieldnames(params);
                    defaults = {};
                    for i = 1:numel(labels)
                        defaults{i} = char(string(params.(labels{i})));
                    end
                    answer = inputdlg(labels, 'DISC', 1, defaults);
                    if isempty(answer)
                        return
                    end
                    for i = 1:numel(labels)
                        if isnumeric(params.(labels{i}))
                            params.(labels{i}) = str2num(answer{i});
                        else
                            params.(labels{i}) = answer{i};
                        end
                    end
                elseif method == "ChangePoint"
                    answer = inputdlg('Threshold ', 'ChangePoint', 1);
                    if isempty(answer)
                        return
                    end
                    params = str2num(answer{1});
%                 elseif method == "MinStepHeight"
%                     answer = inputdlg('Min step height ', 'Min Step Height', 1);
%                     if isempty(answer)
%                         return
%                     end
%                     params = str2num(answer{1});
                end
            end
            if ~isfield(obj.Data, 'yideal')
                [obj.Data.yideal] = deal([]);
            end
            wb = waitbar(0, 'Idealizing...');
            for i = 1:numel(tsIndices)
                j = tsIndices(i);
                % get data indices included in idealization
                if obj.numXROIs() && obj.ui.XROIsButton.Value
                    [~,idx] = obj.XROIMgr.xindices(obj.Data(j).xdata);
                else
                    idx = 1:length(obj.Data(j).ydata);
                end
                idx(isnan(obj.Data(j).ydata(idx))) = [];
                if isempty(idx)
                    continue
                end
                if method == "DISC"
                    try
                        obj.Data(j).yideal = nan(size(obj.Data(j).ydata));
                        obj.Data(j).yideal(idx) = runDISC(obj.Data(j).ydata(idx), params).ideal;
                    catch err
                        disp(err);
                        msgbox("!!! Requires DISC. See https://github.com/ChandaLab/DISC", ...
                            'DISC');
                        return
                    end
                elseif method == "MDL"
                    try
                        obj.Data(j).yideal = nan(size(obj.Data(j).ydata));
                        [breaks, ideal, steplength, stepvalue, jumps] = scan_for_breaks(obj.Data(j).ydata(idx), 1);
                        obj.Data(j).yideal(idx) = ideal;
                    catch err
                        disp(err);
                        msgbox("!!! Requires Scan_for_breaks by Jakob Dreyer. See Add-On Explorer.", ...
                            'DISC');
                        return
                    end
                elseif method == "ChangePoint"
                    N = length(obj.Data(j).ydata(idx));
                    ipt = findchangepts(obj.Data(j).ydata(idx), 'MinThreshold', params);
                    ipt = [1; ipt; N+1];
                    obj.Data(j).yideal = nan(size(obj.Data(j).ydata));
                    jdx = find(idx);
                    for k = 1:numel(ipt)-1
                        kdx = jdx(ipt(k):ipt(k+1)-1);
                        obj.Data(j).yideal(kdx) = mean(obj.Data(j).ydata(kdx));
                    end
%                 elseif method == "MinStepHeight"
%                     obj.Data(j).yideal(idx) = aggSmallSteps(obj.Data(j).ydata(idx), obj.Data(j).yideal(idx), params);
                end
                % only update waitbar at most 10 times as it is slow
                if mod(i, ceil(numel(tsIndices) / 10)) == 0
                    waitbar(double(i) / numel(tsIndices), wb);
                end
            end
            close(wb);
            visNames = obj.visibleNames();
            if ~ismember("ideal", visNames)
                visNames(end+1) = "ideal";
                obj.setVisibleNames(visNames);
            end
            obj.replot();
        end
        
        function idealizeInAxes(obj, method, ax)
            group = obj.axesGroupIndex(ax);
            visTraces = obj.visibleTsPerGroup();
            obj.idealize(visTraces{group}, method);
        end
        
        function idealizeGroup(obj, method, group)
            traces = obj.tsPerGroup();
            obj.idealize(traces{group}, method);
        end
        
        % Set selected portion of idealization to uniform value equal to
        % the first selected idealization point.
        function joinIdeal(obj, ax)
            if obj.numXROIs() == 0 || ~obj.ui.XROIsButton.Value || ~isfield(obj.Data, 'yideal')
                return
            end
            tsi = obj.visibleTsInAxes(ax);
            for i = 1:numel(tsi)
                j = tsi(i);
                [x,y] = obj.getXY(obj.Data(j), "ideal");
                if isempty(y) || numel(x) ~= numel(y)
                    continue
                end
                tf = ~isnan(y);
                if ~any(tf)
                    continue
                end
                x = x(tf);
                y = y(tf);
                indPerROI = obj.XROIMgr.xindices(x);
                for k = 1:numel(indPerROI)
                    indk = indPerROI{k};
                    y(indk) = y(indk(1));
                end
                obj.Data(j).yideal(tf) = y;
            end
            obj.replot();
        end
    end
    
    methods (Static)
        function tsa = mytest()
            tsa = TimeSeriesAnalyzer();
            tmp = load('test_traces.mat');
            tsa.setData(tmp.data);
            for i = 1:numel(tsa.Data)
                tsa.Data(i).yraw = tsa.Data(i).yraw - 10;
                tsa.Data(i).yideal = repmat(mean(tsa.Data(i).ydata), size(tsa.Data(i).ydata));
            end
            tsa.groupEveryN(3);
        end
        function tsa = test()
            tsa = TimeSeriesAnalyzer();
            data = {};
            for i = 1:3
                n = 1e5;
                data{i} = sin(2 * pi * 5/n * (0:n-1))' + rand(n, 1);
            end
            for i = 4:6
                n = 1e7;
                data{i} = sin(2 * pi * 7/n * (0:n-1))' + rand(n, 1);
            end
            tsa.setData(data);
            tsa.groupEveryN(3);
            tsa.setGroupLabels(["Current, pA", "Voltage, mV", "Temperature, C"]);
            tsa.setXLabel("Time, s");
%             for i = 1:12
%                 ui.ts(i).yraw = rand(100 * i,1);
%                 ui.ts(i).yideal = rand(100 * i,1);
%             end
            tsa.replot();
            tsa.autoscaleXY();
        end
    end
end

