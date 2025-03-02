% file slog.m
% This is a function used to print log.

% Log Type
% E, errors, this type will break the code execution.
% W, warnings, this type won't break the code execution, it is used to warn the behavior that may create problems
% D, debug notice, this type will display information for debugging;
% B, basic log.

function slog(message, varargin)

slogdisplaylevel = {'ERR', 'WRN', 'DBG', 'BSCIF'};
textColors = {[1 0 0], [1, 0.5, 0], [0 0.7 1],[0 0 0]}; % Black, orange, red, cyan

persistent LogToFile
persistent LogFile
persistent LogLevel
persistent LogInBlack

if isempty(LogToFile)
	LogToFile = getpref('slog', 'LogToFile', false);
	LogFile = getpref('slog', 'LogFile', 'slog.txt');
	LogLevel = getpref('slog', 'LogLevel', length(slogdisplaylevel));
	LogInBlack = getpref('slog', 'LogInBlack', 0);
end

if LogLevel < 1
	error('slog log level must >= 1');
end

Judge = 0;
if nargin>1 && ischar(varargin{1})
	[Judge, LogTypeIdx] = ismember(varargin{1}, slogdisplaylevel);
	if Judge
		LogType = varargin{1};
		ArgIdx = 2;
	end
end

if nargin == 1 || ~Judge
	LogType = 'BASICINFO';
	LogTypeIdx = find(strcmp(slogdisplaylevel, LogType));
	ArgIdx = 1;
end

if LogTypeIdx > LogLevel
	return
end

currentTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');

db = dbstack(1);
if ~isempty(db)
	callerName = strsplit(db(1).name, '.');
	callerName = callerName(1);
	callerName = callerName{:};
else
	callerName = 'main script';
end

if strcmpi(callerName,'unit')
   callerName = evalin('caller','class(obj)');
end

if strcmp(LogType, 'ERROR')
	message = [message '\n\nTo get help' callerName];
end

message = strrep(message, '%%', '%%%%');
message = strrep(message, '\\', '\\\\');

if nargin > 1
    %What is passed as string argument to robolog (%s) should be treated as string
    %but since we are passing it into sprintf twice we need to escape % and \\.
    %There cannot be escape sequences into the string arguments.
    idx = cellfun(@isstr, varargin);
    varargin(idx) = strrep(varargin(idx), '%', '%%');
    varargin(idx) = strrep(varargin(idx), '\', '\\');
    message=sprintf(['(S %s)\tIn %s: ', message], LogType, callerName, varargin{ArgIdx:end});
else
    message=sprintf(['(S %s)\tIn %s: ', message], LogType, callerName);
end

prefixLength = 16; % length('(Robo NFO) In : ')
spacedNewLine = ['\n' repmat(' ', 1, prefixLength + length(callerName))];
LogMsg = strrep(message, sprintf('\n'), spacedNewLine);
LogMsg = [LogMsg '\n'];

if LogToFile
    fileId = fopen(LogFile, 'a');
    fprintf(fileId, [currentTime ':\t' LogMsg]);
    fclose(fileId);
end

if strcmp(LogType, 'ERR')
    me = MException('slog:genericError', LogMsg);
    throwAsCaller(me);
else

if LogToFile == 2 || LogToFile == 0
    if strcmp(LogType, 'BASICINFO') || LogInBlack
        fprintf(1, LogMsg);
    else
        cprintf(textColors{LogTypeIdx}, LogMsg);
    end
end
end