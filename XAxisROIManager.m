classdef XAxisROIManager < handle
    % Manage, display and interact with a set of x-axis range ROIs.
    %
    % An ROI is a range in the x-axis, e.g. [xmin xmax]
    % ROIs are depicted graphically as interactive Rectangles that shade
    % the selected ranges in the plot axes.
    %
    % Add ROIs interactively:
    % .addROI() --> click and drag to draw a range in the current axes
    % .addROI(x0) --> drag to draw a range in the current axes with one edge at x0
    %
    % Get/Set ROIs programatically via .XLims:
    % .XLims --> an Mx2 matrix of [xmin xmax] for each of M ranges
    % e.g. xlims = obj.XLims; % grab the ROI lims
    % e.g. obj.XLims = [...]; % set the ROIs
    %
    % Standard Rectangle interactions for dragging/resizing ROIs
    % interactively.
    %
    % ROIs provide a context menu for interactions.
    % You can customize the context menu for your needs.
    % .ROIs(...).UIContextMenu = ...
    
    properties
        % List of plot axes in which ROIs should be displayed and managed.
        % For multiple axes, all axes share the same ROIs. Thus, changing an
        % ROI in one axes changes it for all.
        Axes = gobjects(0);
        
        % NxM matrix of Rectangles for visualizing the ROIs.
        % ROis(i,j) -> ith axes, jth roi
        ROIs = images.roi.Rectangle.empty;
        
        % ROI attributes
        % Setting these will set them for all ROIs
        ROIProperties = struct('Color', [1 0 0], 'FaceAlpha', 0.05, 'Visible', true);
    end
    
    properties (Dependent)
        % This makes querying and setting the ROIs very simple.
        % XLims(j,:) = [xmin xmax] for ROIs(:,j)
        XLims
    end
    
    events
        NumberOfROIsChanged
        ROIChanged
    end
    
    methods
        function obj = XAxisROIManager(axes)
            obj.Axes = unique(axes);
        end
        function delete(obj)
            % Removing all the axes will clean everything up
            obj.Axes = gobjects(0);
        end
        
        function xlims = get.XLims(obj)
            allpos = vertcat(obj.ROIs(1,:).Position);
            xlims = cumsum(allpos(:,[1 3]), 2);
        end
        function set.XLims(obj, xlims)
            delete(obj.ROIs);
            obj.ROIs = gobjects(0);
            N = size(xlims, 1);
            for j = 1:N
                obj.addROI(xlims(j,:));
            end
        end
        
        % Set managed axes
        function set.Axes(obj, axes)
            % Current ranges
            try
                xlims = obj.XLims;
            catch
                xlims = [];
            end
            
            % Delete current ROIs
            delete(obj.ROIs);
            obj.ROIs = gobjects(0);
            
            % Remove current axes
            for i = 1:numel(obj.Axes)
                ax = obj.Axes(i);
                % remove axes data
                if isvalid(ax)
                    rmappdata(ax, 'XAxisROIManagerData');
                end
            end
            obj.Axes = gobjects(0);
            
            % Set new axes
            obj.Axes = axes(isvalid(axes));
            for i = 1:numel(obj.Axes)
                ax = obj.Axes(i);
                % store listeners and data in axes
                axdata.obj = obj;
                axdata.YRulerMarkedCleanListener = ...
                    addlistener(ax.YRuler, 'MarkedClean', @obj.onAxesYRulerChanged);
                setappdata(ax, 'XAxisROIManagerData', axdata);
            end
            
            % Reset ROIs
            obj.XLims = xlims;
        end
        
        % Union of figures for all managed axes
        function figs = figures(obj)
            figs = gobjects(0);
            if ~isempty(obj.Axes)
                for i = 1:numel(obj.Axes)
                    ax = obj.Axes(i);
                    if isvalid(ax)
                        figs = [figs ancestor(ax, 'Figure', 'toplevel')];
                    end
                end
            end
        end
        
        % Set ROI attributes
        function set.ROIProperties(obj, props)
            names = fieldnames(props);
            for i = 1:numel(names)
                [obj.ROIs(:).(names{i})] = deal(props.(names{i}));
            end
            obj.ROIProperties = props;
        end
        
        % Number of ROIs
        function n = numROIs(obj)
            n = size(obj.ROIs, 2);
        end
        
        % ROI indices with respect to input x
        function [indPerROI, indAllROI] = xindices(obj, x)
            indPerROI = {};
            indAllROI = [];
            xlims = obj.XLims;
            for j = 1:size(xlims,1)
                xmin = xlims(j,1);
                first = find(x >= xmin, 1);
                if isempty(first)
                    indPerROI{j} = [];
                    continue
                end
                xmax = xlims(j,2);
                afterLast = find(x(first:end) > xmax, 1);
                if isempty(afterLast)
                    ind = first:length(x);
                elseif afterLast == 1
                    indPerROI{j} = [];
                    continue
                else
                    last = first - 1 + afterLast - 1;
                    ind = first:last;
                end
                indPerROI{j} = ind;
                indAllROI = union(indAllROI, ind);
            end
            indAllROI = unique(indAllROI);
        end
        
        % Add ROIs
        % addROI() --> click and drag to draw ROI
        % addROI(x0) --> drag to draw ROI with one edge at x0
        % addROI(xrange) --> add ROI with input range (noninteractive)
        function addROI(obj, xlim)
            if isempty(obj.Axes)
                obj.Axes = gca;
            end
            % add ROI to all axes
            j = size(obj.ROIs, 2) + 1;
            obj.ROIs(1,j) = images.roi.Rectangle(obj.Axes(1), ...
                'Color', [1,0,0], ...
                'FaceAlpha', 0.05, ...
                'SelectedColor', 'none', ...
                'LineWidth', 1e-9, ...
                'LabelVisible', 'off', ...
                'Visible', true, ...
                'InteractionsAllowed', 'all');
            names = fieldnames(obj.ROIProperties);
            for i = 1:numel(names)
                [obj.ROIs(1,j).(names{i})] = deal(obj.ROIProperties.(names{i}));
            end
            for i = 1:numel(obj.Axes)
                if i > 1
                    obj.ROIs(i,j) = copyobj(obj.ROIs(1,j), obj.Axes(i));
                end
                roi = obj.ROIs(i,j);
                setappdata(roi, 'MovingROIListener', ...
                    addlistener(roi, 'MovingROI', @obj.onROIChanged));
                setappdata(roi, 'ROIMovedListener', ...
                    addlistener(roi, 'ROIMoved', @obj.onROIChanged));
                roi.UIContextMenu = obj.ROIContextMenu(roi);
            end
            % draw or set ROI range in current axes
            ica = find(obj.Axes == gca);
            if isempty(ica)
                ica = 1;
            end
            ax = obj.Axes(ica);
            roi = obj.ROIs(ica,j);
            if ~exist('xlim', 'var') || isempty(xlim)
                draw(roi);
                if ~isvalid(roi)
                    delete(obj.ROIs(:,j));
                    obj.ROIs(:,j) = [];
                    return
                end
                xlim = cumsum(roi.Position([1 3]));
            elseif numel(xlim) == 1
                beginDrawingFromPoint(roi, [xlim ax.YLim(1)]);
                if ~isvalid(roi)
                    delete(obj.ROIs(:,j));
                    obj.ROIs(:,j) = [];
                    return
                end
                xlim = cumsum(roi.Position([1 3]));
            end
            if diff(xlim) == 0
                % delete added ROIs if they have zero width
                delete(obj.ROIs(:,j));
                obj.ROIs(:,j) = [];
                return
            else
                % update width of all added ROIs in all axes
                for i = 1:numel(obj.Axes)
                    ax = obj.Axes(i);
                    obj.ROIs(i,j).Position = [xlim(1) ax.YLim(1) diff(xlim) diff(ax.YLim)];
                end
                notify(obj, 'NumberOfROIsChanged');
            end
        end
        
        % Merge overlapping ROIs
        function xlims = mergedXLims(obj)
            xlims = obj.XLims;
            xmin = xlims(:,1);
            xmax = xlims(:,2);
            [xmin, idx] = sort(xmin);
            xmax = xmax(idx);
            i = 1;
            while i+1 <= length(xmin)
                if xmax(i) >= xmin(i+1)
                    % merge i and i+1
                    xmax(i) = max(xmax(i:i+1));
                    xmin(i+1) = [];
                    xmax(i+1) = [];
                else
                    i = i+1;
                end
            end
            xlims = [xmin xmax];
        end
        function mergeOverlappingROIs(obj)
            obj.XLims = obj.mergedXLims();
        end
        
        % Remove ROIs
        function deleteROI(obj, roi)
            if ~exist('roi', 'var')
                roi = gco;
            end
            idx = find(obj.ROIs(:) == roi);
            if ~isempty(idx)
                [i,j] = ind2sub(size(obj.ROIs), idx);
                delete(obj.ROIs(:,j));
                obj.ROIs(:,j) = [];
                notify(obj, 'NumberOfROIsChanged');
            end
        end
        function clearROIs(obj)
            if questdlg('Clear all ROIs?', 'Clear ROIs?') == "Yes"
                delete(obj.ROIs);
                obj.ROIs = gobjects(0);
                notify(obj, 'NumberOfROIsChanged');
            end
        end
        
        % Update ROI in all axes when moved/resized in any axes
        function onROIChanged(obj, roi, varargin)
            ax = roi.Parent;
            roi.Position([2 4]) = [ax.YLim(1) diff(ax.YLim)];
            xw = roi.Position([1 3]);
            idx = find(obj.ROIs(:) == roi);
            if ~isempty(idx)
                [i0,j0] = ind2sub(size(obj.ROIs), idx);
                for i = 1:size(obj.ROIs, 1)
                    if i ~= i0
                        ax = obj.Axes(i);
                        obj.ROIs(i,j0).Position = [xw(1) ax.YLim(1) xw(2) diff(ax.YLim)];
                    end
                end
            end
            notify(obj, 'ROIChanged');
        end
        
        % Update ROIs to fill vertical height of axes
        function onAxesYRulerChanged(obj, yruler, varargin)
            ax = yruler.Parent;
            i = find(obj.Axes == ax);
            if isempty(i)
                return
            end
            yh = [ax.YLim(1) diff(ax.YLim)];
            for j = 1:size(obj.ROIs, 2)
                obj.ROIs(i,j).Position([2 4]) = yh;
            end
        end
        
        % ROI context menu
        function menu = ROIContextMenu(obj, roi, varargin)
            if ~exist('roi', 'var')
                roi = gco;
            end
            fig = ancestor(roi, 'Figure', 'toplevel');
            menu = uicontextmenu(fig);
            
            uimenu(menu, 'Label', 'Delete This ROI', ...
                'Callback', @(varargin) obj.deleteROI(roi));
            uimenu(menu, 'Label', 'Clear All ROIs', ...
                'Callback', @(varargin) obj.clearROIs());
            uimenu(menu, 'Label', 'Merge Overlapping ROIs', ...
                'Callback', @(varargin) obj.mergeOverlappingROIs());
        end
    end
    
    methods (Static)
        function mgr = test()
            ax1 = subplot(2,1,1);
            ax2 = subplot(2,1,2);
            mgr = XAxisROIManager([ax1 ax2]);
            mgr.addROI();
            mgr.addROI();
        end
    end
end

