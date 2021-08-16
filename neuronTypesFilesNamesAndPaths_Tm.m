% Author: Yu-Chuan Chen IOC AS Tawain
% Last update July 19, 2019

matFilePath = '/Yu-Chuan/MATLAB/Working_Tm/' ;

%%  Give paths to *.swc.txt files for all neurons. These files have the node ID data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadingPath = '/Yu-Chuan/MATLAB/Rigid_Registered_swc/';
directory_out = '/Yu-Chuan/MATLAB/Rigid_Registered_swc/';
file_end_in = '*.swc';

date = '0729'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Define folder name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for classified by Yu-Chuan (classified = 1 )
% for classified by original file name (classified =2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%% #3
lab_name = 'Lee'
folderName = {'Tm2_'; 'Tm20_'};
classified = 2 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% #4
lab_name = 'LeeCH'
folderName = {'Tm1'; 'Tm2'; 'Tm9'; 'Tm20'};
classified = 1
%}
%%%%%%%%
directory_in  = strcat(loadingPath, lab_name, '/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
numberOfClasses = size( folderName, 1 ); 
numberOfNeuronsInClass = zeros( numberOfClasses, 1 ) ;
totalNumberOfNeurons = 0 ;

neuron_names = folderName';

if classified == 1 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   for classified by Yu-Chuan (classified = 1 )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for cx = 1 : numberOfClasses
    
        files = dir(strcat(directory_in,'/', neuron_names{cx},'/' ,file_end_in));
    
        %%%Get number of neurons measured for each class
        numberOfNeuronsInClass( cx ) = size( files, 1 );
    
        totalNumberOfNeurons = totalNumberOfNeurons + size( files, 1 ) ;
 
    end
   
else   
    
    for cx = 1 : numberOfClasses
    
         files = dir(strcat(directory_in,'*', neuron_names{cx} ,file_end_in));
    
        %%%Get number of neurons measured for each class
        numberOfNeuronsInClass( cx ) = size( files, 1 );
    
        totalNumberOfNeurons = totalNumberOfNeurons + size( files, 1 ) ;
 
    end

end

maxNumberOfNeuronsInAnyClass = max( numberOfNeuronsInClass ) ;

%%
matFilePathOuput =  strcat(matFilePath, 'mat_files/');

save(strcat(matFilePathOuput, 'workingDateAndLab', '.mat'), 'lab_name', 'date')

save( strcat( matFilePathOuput , 'neuronTypesFileNamesAndPaths_Pyramidal_', lab_name, '_', date,'.mat' ), 'directory_in', 'directory_out', ...
'file_end_in', 'neuron_names', 'numberOfClasses','numberOfNeuronsInClass', 'maxNumberOfNeuronsInAnyClass', 'totalNumberOfNeurons', 'lab_name',...
'classified')


