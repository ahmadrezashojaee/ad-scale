classdef ResultHandler < handle
    % Class for storing and retrieving simulation results, either in memory or stored to disk
    %
    % SYNOPSIS:
    %   handler = ResultHandler()
    %
    %   handler = ResultHandler('dataPrefix', 'mydata', 'writeToDisk', true);
    % 
    % DESCRIPTION:
    %   This class can be used to store and retrieve simulation results. It 
    %   is somewhat similar to a cell array in use, although more limited.
    %  
    %   Take a look at the class declaration to get more information.
    %
    % EXAMPLE:
    %   % Setup handler
    %   handler = ResultHandler('dataprefix', 'mydata', 'writeToDisk', true);
    %   % Write result
    %   handler{1} = {'hello'};
    %   % Read result from disk and print
    %   disp(handler{1})
    %
    % RETURNS:
    %   Class instance that in some limited aspects acts like a cell array

    properties
        % Boolean indicating if results should be written to disk
        writeToDisk
        % Boolean indicating if the class should store results in memory
        storeInMemory
        
        % Directory where data in general is stored.
        dataDirectory
        % The folder under directory where we will store results. Will be
        % created if it does not exist.
        dataFolder
        % Data will be stored in the format <dataPrefix><index> so that
        % calling handler{51} will store a file as "state51" if dataprefix
        % is 'state'
        dataPrefix
        % Flags passed on to MATLAB builtin 'save'. Consider '-v7' if
        % results are huge.
        saveflags
        % Clear directory on startup
        cleardir
        
        % Internal data storage
        data
        % Flag indicating if verbose output is on. Will output storage and
        % reading notifications if enabled.
        verbose
    end
    
    methods
        
        function handler = ResultHandler(varargin)
            handler.writeToDisk = true;
            handler.storeInMemory = false;

            handler.dataDirectory = fullfile(mrstOutputDirectory(), 'tmp');
            handler.dataPrefix = 'state';
            handler.dataFolder = 'cache';
            handler.saveflags = '';
            handler.cleardir = false;
            
            handler.verbose = mrstVerbose();
            
            handler = merge_options(handler, varargin{:});
            
            if ~(handler.writeToDisk || handler.storeInMemory)
                warning('ResultHandler:Noop',...
                    ['Configured without any storage targets!'...
                    ' All data will be discarded!']);
            end
            
            % Storage stuff
            handler.data = {};
            if handler.writeToDisk
                p = handler.getDataPath();
                if ~(exist(p, 'dir') == 7)
                    ok = mkdir(p);
                    if ~ok
                        error(['Unable to create output directory! ',...
                        'Ensure that you have write permissions to ''',...
                        handler.dataDirectory, ...
                        '''', ...
                        ])
                    end
                end
                try
                    d = ls(fullfile(p, [handler.dataPrefix, '*.mat']));
                    if ~isempty(d)
                        if ~handler.cleardir
                            dispif(handler.verbose > 1, ...
                                ['Input directory not clean, consider calling', ...
                                ' ''resetData'' to remove pre-existing results.\n']);
                        else
                            handler.resetData();
                        end
                    end
                catch e
                    if ~strcmp(e.identifier, 'MATLAB:ls:OSError')
                        % ls returns error for empty dir on some systems
                        rethrow(e)
                    end
                end
            end
            
        end
        
        function n = numelData(handler)
            if handler.writeToDisk
                n = numel(handler.getValidIds);
            elseif handler.storeInMemory
                n = numel(handler.data); 
            else
                n = 0;
            end
        end
        
        function varargout = sizeData(handler, dim)
            N = [numelData(handler), 1];
            if nargin == 1
                varargout{1} = N;
            else
                varargout{1} = N(dim);
            end
        end
        
        function varargout = subsref(handler, s)
            switch s(1).type
                case '.'
                    % Methods and so on can be dispatched to matlab
                    [varargout{1:nargout}] = builtin('subsref', handler, s);
                case {'()', '{}'}
                    if handler.storeInMemory
                        [varargout{1:nargout}] = builtin('subsref', handler.data, s);
                        return
                    end
                    
                    if handler.writeToDisk
                        sub = s(1).subs{1};
                        if ischar(sub) && strcmp(sub, ':')
                            sub = 1:handler.numelData();
                        end
                        tmp = handler.readFromFile(sub);
                        s(1).subs{1} = 1:numel(sub);
                        [varargout{1:nargout}] = builtin('subsref', tmp, s);
                        return
                    end
            end
            
        end
        
        function handler = subsasgn(handler, s, v)
            switch s.type
                case '.'
                    handler = builtin('subsasgn', handler, s, v);
                case {'()', '{}'}
                    if iscell(s.subs)
                        assert(numel(s.subs) == 1);
                        s.subs = s.subs{1};
                    end
                    
                    if ~iscell(v)
                        if numel(v) == 1
                            [tmp{1:numel(s.subs)}] = deal(v);
                            v = tmp;
                        else
                            v = arrayfun(@(x) x, v, 'UniformOutput', false);
                        end
                    end
                    
                    if numel(v) ~= numel(s.subs)
                        error('ResultHandler:dimMismatch',...
                            'Left hand side and right hand side must have equal dimensions')
                    end
                    
                    if handler.writeToDisk
                        for i = 1:numel(s.subs)
                            handler.writeToFile(v{i}, s.subs(i));
                        end
                    end
                    
                    if handler.storeInMemory
                        for i = 1:numel(s.subs)
                            handler.data{s.subs(i)} = v{i};
                        end
                    end
                otherwise
                
            end
        end
        
        function data = readFromFile(handler, ids)
            n = numel(ids);
            data = cell(1, n);
            for i = 1:n
                p = handler.getDataPath(ids(i));
                tmp = load(p, 'data');
                data{i} = tmp.data;
            end
        end
        
        function ids = getValidIds(handler)
            ids = [];
            
            if handler.writeToDisk
                p = handler.getDataPath();
                
                l = ls(p);
                if size(l, 1) > 1
                    % Windows behavior
                    % Pad with a space at the end and reshape into a single
                    % long string.
                    pad = repmat(' ', size(l, 1), 1);
                    l = reshape([l, pad]', 1, []);
                end
                l = regexp(l,'\s+','split');
                for i = 1:numel(l)
                    line = l{i};
                    [s, e] = regexp(line, [handler.dataPrefix, '\d+']);
                    if isempty(s); continue; end
                    ids = [ids, str2double(line((s+numel(handler.dataPrefix)):e))];
                end
                ids = sort(ids);
                return
            end
            
            if handler.storeInMemory
                ids = find(~cellfun(@isempty, handler.data));
                return
            end
        end
        
        function handler = writeToFile(handler, data, id) %#ok
            p = handler.getDataPath(id);
            save(p, 'data', handler.saveflags);
            dispif(handler.verbose, 'Writing data to %s\n', p);
        end
        
        function handler = deleteFile(handler, id)
            p = handler.getDataPath(id);
            delete(p, 'data');
            dispif(handler.verbose, 'Deleting data at %s\n', p);
        end
        
        function p = getDataPath(handler, i)
            
            p = fullfile(handler.dataDirectory, handler.dataFolder);
            if nargin > 1
                assert(numel(i) == 1 && isnumeric(i));
                p = fullfile(p, [handler.dataPrefix, num2str(i), '.mat']);
            end
        end
        
        function resetData(handler, subs)
            if nargin == 1
                handler.data = {};
                if handler.writeToDisk
                    p = handler.getDataPath();
                    fp = fullfile(p, [handler.dataPrefix, '*.mat']);
                    delete(fp);
                end
            else
                if handler.storeInMemory
                    [handler.data{subs}] = deal([]);
                end
                if handler.writeToDisk
                    if islogical(subs)
                        subs = find(subs);
                    end
                    p = handler.getDataPath();
                    for i = 1:numel(subs)
                        fp = fullfile(p, [handler.dataPrefix, num2str(subs(i)), '.mat']);
                        delete(fp);
                    end
                end
            end
        end
        
        function [state0, restartStep] = getRestart(handler)
            nd = handler.numelData();
            if nd > 0
                s = substruct('{}', {nd});
                state0 = handler.subsref(s);
                restartStep = nd + 1;
            else
                state0 = [];
                restartStep = 1;
            end
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


