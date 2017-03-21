% Performs automatic note detection using the 'Notes' algorithm
% from the pYin VAMP plugin, by Matthias Mauch.
% The VAMP plugin is invoked directly using Sonic Annotator.
% The output is a midi file with the same name as the input wav file.
function pYinNotes(wavFileFolder, wavFileName, outputFilePath, pYinPath, sonicAnnotatorPath, transformFilePath)

if isunix(), sep = '/'; else sep = '\'; end

if nargin < 4
    pYinPath = ['..' sep 'resources' sep 'pYin'];
    sonicAnnotatorPath = ['..' sep 'resources' sep 'sonic-annotator'];
    transformFilePath = [sonicAnnotatorPath sep 'transform.rdf'];
end

if isunix
    % set VAMP_PATH to include pYinPath
    exportCommand = ['export VAMP_PATH=$(pwd)/' pYinPath '; '];
    moveCommand = 'mv ';
    if ismac
        platFolder = 'macosx';
    else
        platFolder = 'debian64';
    end
else
    exportCommand = ['set VAMP_PATH=%cd%\' pYinPath ' & '];
    moveCommand = 'move ';
    platFolder = 'win32';
end
system([exportCommand sonicAnnotatorPath sep platFolder sep 'sonic-annotator -t ' ...
     transformFilePath ' -w midi ' wavFileFolder sep wavFileName]);

if nargin >= 3
    system([moveCommand wavFileFolder sep 'performance.mid ' ...
        outputFilePath]);
end