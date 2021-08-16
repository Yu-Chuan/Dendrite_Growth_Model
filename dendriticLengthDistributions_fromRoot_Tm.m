
% Author: Yu-Chuan Chen IOC Academia Sinca Tawain
% Last update: July 19 2019

%%%%%%%%%%
%% This script is used to extract the node-to-node lengths along the dendrites.
%% Rigid data used.
%%%%%%%%%%

%% Path to file with the path data
matFilePath = '/Yu-Chuan/MATLAB/Working_Tm/mat_files/';

load( strcat( matFilePath, 'workingDateAndLab', '.mat'))

load( strcat( matFilePath, 'neuronTypesFileNamesAndPaths_Pyramidal_', lab_name, '_', date,'.mat' ) )

data = {};
trees = {};
Tm = {};
newTrees = {};

start_trees;

 
%%%%%
%% Find largest branch length (which will be used later in defining the size of ceretain arrays).
%% Also, find the largest order of branching
%%%%%
maxDendriticBranchLength = -10000000.0 ;
maximumBranchingOrderForClass = zeros( numberOfClasses, 1 ) ;

%%%%%
%% check the  neuron type (Yu)
%%%%%
numOfSoma = zeros(numberOfClasses, maxNumberOfNeuronsInAnyClass);
numOfAxion = zeros(numberOfClasses, maxNumberOfNeuronsInAnyClass);
numOfBasal_dendrite =  zeros(numberOfClasses, maxNumberOfNeuronsInAnyClass);
numOfApical_dendrite =  zeros(numberOfClasses, maxNumberOfNeuronsInAnyClass);
numOfOther = zeros(numberOfClasses, maxNumberOfNeuronsInAnyClass);

numberOfSegmentsForNeuron = {};
numberOfSegmentsInNeuron = zeros( numberOfClasses, maxNumberOfNeuronsInAnyClass ) ;

structureFlag = zeros( numberOfClasses, maxNumberOfNeuronsInAnyClass );

%%%time Start
timer1 = tic;

for cx = 1:numberOfClasses
   
    loadingPath = strcat(directory_in);
        
    for nx = 1: numberOfNeuronsInClass(cx)
        
        if classified == 1
        files = dir(strcat(directory_in, '/', neuron_names{cx},'/' ,file_end_in));
        daa = importdata( strcat( loadingPath, neuron_names{cx},'/' ,files(nx).name ) ) ;
        trees{cx, nx} = load_tree(  strcat(loadingPath , neuron_names{cx},'/' ,files(nx).name ) );
        
        else
        files = dir(strcat(loadingPath,'*', neuron_names{cx}, file_end_in));
        daa = importdata( strcat( loadingPath ,files(nx).name ) ) ;
        trees{cx, nx} = load_tree(  strcat(loadingPath ,files(nx).name ) );

        end % end of if

        data{cx, nx} = daa; 

        % check the structure type  number

        soma = data{cx, nx}.data(data{cx, nx}.data(:,2) == 1);
        numOfSoma(cx, nx) = sum(data{cx, nx}.data(:,2) == 1);
        
        axon = data{cx, nx}.data(data{cx, nx}.data(:,2) == 2);
        numOfAxion(cx, nx)  = sum(data{cx, nx}.data(:,2) == 2);
        
        dendrite = data{cx, nx}.data(data{cx, nx}.data(:,2) == 3);
        numOfBasal_dendrite(cx, nx) = sum(data{cx, nx}.data(:,2) == 3);
        
        apical_dendrite = data{cx, nx}.data(data{cx, nx}.data(:,2) == 4);
        numOfApical_dendrite(cx, nx)  = sum(data{cx, nx}.data(:,2) == 4);
        
        Other = data{cx, nx}.data(data{cx, nx}.data(:,2) >= 5);
        numOfOther(cx, nx) = sum(data{cx, nx}.data(:,2) >= 5);
        
        
        %%% check the lost %%%
        all = sum([numOfSoma(cx, nx) , numOfAxion(cx, nx) , numOfBasal_dendrite(cx, nx) , ...
            numOfApical_dendrite(cx, nx), numOfOther(cx,nx)]);
        test = all -size(data{cx, nx}.data(:,2),1);
        if test > 0
            disp("lose some segment in", data{cx, nx})
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%segment selection%%%
        
        lenFromRoot = Pvec_tree(trees{cx, nx}); 
        
        if numOfSoma(cx, nx) == size(data{cx, nx}.data, 1)      
            %%% neuron with all dendrite
            lengthBr =  lenFromRoot(B_tree(trees{cx, nx}));
            lengthTe =  lenFromRoot(T_tree(trees{cx, nx}));
            structureFlag(cx, nx) = 1;
           
            
        elseif (numOfAxion(cx, nx) + numOfSoma(cx, nx) ==  all )
            %%% dendrite and soma   
            lengthBr = lenFromRoot(data{cx, nx}.data(:,2) == 2 & B_tree(trees{cx, nx}));
            lengthTe = lenFromRoot(data{cx, nx}.data(:,2) == 2 & T_tree(trees{cx, nx}));
            structureFlag(cx, nx)  = 2;

        
        elseif (numOfApical_dendrite(cx, nx) == 0 && numOfBasal_dendrite(cx, nx) ~= 0)
            %%% Basal dendrite and soma
            lengthBr =  lenFromRoot(data{cx, nx}.data(:,2) == 3 & B_tree(trees{cx, nx}));
            lengthTe = lenFromRoot(data{cx, nx}.data(:,2) == 3 & T_tree(trees{cx, nx}));
            structureFlag(cx, nx)  = 3;
             
        
        elseif (numOfApical_dendrite(cx, nx) ~= 0 && numOfBasal_dendrite(cx, nx) == 0)
            %%% Apical dendrite and soma     
            lengthBr =  lenFromRoot(data{cx, nx}.data(:,2) == 4 & B_tree(trees{cx, nx}));
            lengthTe = lenFromRoot(data{cx, nx}.data(:,2) == 4 & T_tree(trees{cx, nx}));
            structureFlag(cx, nx)  = 4;
            
        
        else
            %%% Apical dendrite        
            lengthBr_A =  lenFromRoot(data{cx, nx}.data(:,2) == 4 & B_tree(trees{cx, nx}));
            lengthTe_A = lenFromRoot(data{cx, nx}.data(:,2) == 4 & T_tree(trees{cx, nx}));
            %%% Basal dendrite 
            lengthBr_B =  lenFromRoot(data{cx, nx}.data(:,2) == 3 & B_tree(trees{cx, nx}));
            lengthTe_B = lenFromRoot(data{cx, nx}.data(:,2) == 3 & T_tree(trees{cx, nx}));
            structureFlag(cx, nx)  = 5;
                        
        end
    
       %%% For different dentrite type%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       if structureFlag(cx, nx) == 5 
        %%% Apical and basal dendrite 
           
        Br_structureID_A = zeros(length(lengthBr_A), 1) + 4;
        Br_structureID_B = zeros(length(lengthBr_B), 1) + 3;
        Br_structureID = [Br_structureID_A; Br_structureID_B];
        lengthBr = [lengthBr_A; lengthBr_B];
        
        Te_structureID_A = zeros(length(lengthTe_A), 1) + 4;
        Te_structureID_B = zeros(length(lengthTe_B), 1) + 3;
        Te_structureID = [Te_structureID_A; Te_structureID_B];
        lengthTe = [lengthTe_A; lengthTe_B]; 
        %%%
        numberOfSegmentsInNeuronB( cx, nx ) = length([lengthBr_B; lengthTe_B]);
        
        numberOfSegmentsInNeuronA( cx, nx ) = length([lengthBr_A; lengthTe_A]);
        
       
        %%%
        flagOfBr = zeros(length(lengthBr), 1) + 2;
        Br = [lengthBr, flagOfBr ];
        
        flagOfTe = zeros(length(lengthTe), 1);
        Te = [lengthTe, flagOfTe ];
        
        gene{cx, nx}= [Br; Te];
        
        structureID = [Br_structureID; Te_structureID];  
           
       else
        
        %%%% only one type of dendrite
   
        flagOfBr = zeros(length(lengthBr), 1) + 2;
        Br = [lengthBr, flagOfBr ];
        
        flagOfTe = zeros(length(lengthTe), 1);
        Te = [lengthTe, flagOfTe ];
        
        gene{cx, nx}= [Br; Te];
        
        structureID = zeros(size(gene{cx, nx}, 1), 1) + 2;
        
        %%%
        numberOfSegmentsInNeuron( cx, nx ) = length(gene{cx, nx}) ;
        
        %%%
        
       end
        
        Tm{cx, nx} = struct('length', gene{cx, nx}(:, 1), 'termFlag', gene{cx, nx}(:, 2), 'structureID', structureID);
        
        %numberOfSegmentsForNeuronT( cx, nx ) = size( Tm{cx, nx}.length, 1 ) ;
                
        for ix = 1:  size( Tm{cx, nx}.length, 1 ) 
            
            if maxDendriticBranchLength < Tm{cx, nx}.length( ix )
                
                maxDendriticBranchLength = Tm{cx, nx}.length( ix );
                
            end
           
        end  
        
         disp(files(nx).name); 
        
    end % end of nx
    
end % end of cx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cx = 1: numberOfClasses

    for nx = 1: numberOfNeuronsInClass(cx)
    
        if structureFlag(cx, nx) == 5
        %%%combine two data set
        %%% for basal
            numberOfSegmentsForNeuron{1} = numberOfSegmentsInNeuronB;
        %%% for apical
            numberOfSegmentsForNeuron{2} = numberOfSegmentsInNeuronA;

        else
    
            numberOfSegmentsForNeuron{1} = numberOfSegmentsInNeuron;
    
        end
        
    end % end of nx
    
end % end of cx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Overall time of part1= %f sec\n', toc(timer1));

%%%%%%
%% Set up some auxilliary variables and arrays
%%%%%%

minLen = 0 
maxLen = maxDendriticBranchLength 
deltaLen = 2.0  % 0.2 for fruit fly
numberOfLengthValues = 1 + floor( ( maxLen - minLen )/deltaLen );
numberOfLengthValues = numberOfLengthValues + 1 
%% To account for initial time of t = 0 


lenValues = zeros( numberOfLengthValues, 1 ) ;
lenValues( 1 ) = 0 ;
for ix = 2:numberOfLengthValues 
    
    lenValues( ix ) = lenValues( ix - 1 ) + deltaLen ;

end
    
lowerBoundSHatCurve = -100 ;

%%%%%
%% Gather branching and terminating segment data from the neurons 
%%%%%

numBrNeu = zeros( numberOfClasses, maxNumberOfNeuronsInAnyClass, numberOfLengthValues ) ; 
numTeNeu = zeros( numberOfClasses, maxNumberOfNeuronsInAnyClass, numberOfLengthValues ) ;

numBrWC = zeros( numberOfClasses, numberOfLengthValues ) ;  
numTeWC = zeros( numberOfClasses, numberOfLengthValues ) ;

numBrWC_AB = zeros( 3, numberOfClasses, numberOfLengthValues ) ;  
numTeWC_AB = zeros( 3, numberOfClasses, numberOfLengthValues ) ;

%%% time start
timer2 = tic;

for cx = 1:numberOfClasses
    
    for nx = 1:numberOfNeuronsInClass(cx)
        
        
        for ix = 1: size( Tm{cx, nx}.length, 1 )  
            
            
            lld = Tm{cx, nx}.length( ix ) ;
            
            inx = floor( ( lld - minLen )/deltaLen ) + 1 ;
            inx = inx + 1 ;
            
            if inx < 2
                
                inx = 2 ;
                
            else
                
                if inx > numberOfLengthValues
                    
                    inx = numberOfLengthValues ;
                    
                end
                
            end
            
            if Tm{cx, nx}.termFlag( ix ) == 0 
                
                numTeNeu( cx, nx, inx ) = numTeNeu( cx, nx, inx ) + 1 ;
                
                numTeWC( cx, inx ) = numTeWC( cx, inx ) + 1 ;
                
                if Tm{cx, nx}.structureID(ix) == 2
                    
                   numTeWC_AB(1, cx, inx) =  numTeWC_AB(1, cx, inx) +1;
                    
                elseif Tm{cx, nx}.structureID(ix) == 3
                    
                   numTeWC_AB(2, cx, inx) =  numTeWC_AB(2, cx, inx) +1;
                    
                elseif Tm{cx, nx}.structureID(ix) == 4
                    
                   numTeWC_AB(3, cx, inx) =  numTeWC_AB(3, cx, inx) +1;
                    
                end     
                                     
            else
                
                numBrNeu( cx, nx, inx ) = numBrNeu( cx, nx, inx ) + 1 ;
                
                numBrWC( cx, inx ) = numBrWC( cx, inx ) + 1 ;
                
                if Tm{cx, nx}.structureID(ix) == 2
                    
                   numBrWC_AB(1, cx, inx) =  numBrWC_AB(1, cx, inx) +1;
                    
                elseif Tm{cx, nx}.structureID(ix) == 3
                    
                   numBrWC_AB(2, cx, inx) =  numBrWC_AB(2, cx, inx) +1;
                    
                elseif Tm{cx, nx}.structureID(ix) == 4
                    
                   numBrWC_AB(3, cx, inx) =  numBrWC_AB(3, cx, inx) +1;
                    
                end     
                               
            end
         
        end
        
        disp(files(cx).name);
        
    end

end


fprintf('Overall time of part2= %f sec\n', toc(timer2));

%% define the segment of different region
t_section = 1.0
lenRange = 0: t_section: 200;
numOfBr = {};
numOfTe = {};
groupOfKB = {};
groupOfKT = {};
groupOfAtRisk = {};

AtRisk = {};

numOfBrOfType ={};
numOfTeOfType = {};
numOfAtriskOfType = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate kb and kt (Yu)
%%% time Start
timer3 = tic;
%%%%%%%%%%%%%%%%%%%

for cx = 1:numberOfClasses
    
     atrisk = zeros(numberOfNeuronsInClass(cx),  size(lenRange, 2) -1 );
     Br = {};
     Te = {};
     numOfBr = {};
     numOfTe = {};
     AtRisk = {};
    
    for nx = 1:numberOfNeuronsInClass(cx)
        
        if structureFlag(cx, nx) == 5 
    
            numberOfType = 2;
    
        else  
            numberOfType = 1;

        end
        
        
         for tx = 1: numberOfType
        
            if structureFlag(cx, nx) == 1   
            %%% neuron with all dendrite
                atrisk(nx, 1) = sum(data{cx, nx}.data(: , 7) == 1 & data{cx, nx}.data(:,2) == 1);  

            elseif structureFlag(cx, nx) == 2
            %%% dendrite and soma   
                AA = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 1, 1);
                BB = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 2, 7);
                % find the same element of two sequence
                C = intersect(AA, BB);
                
               
                %{
                AA = data{4, 8}.data(data{4, 8}.data(: , 2) == 1, 1);
                BB = data{4, 8}.data(data{4, 8}.data(: , 2) == 2, 7);
                if BB(1) == -1
                   disp('root is dendrite');
                end
                
                % find the same element of two sequence
                C = intersect(AA, BB);
                
                AA = data{4, 11}.data(data{4, 11}.data(: , 2) == 1, 1);
                BB = data{4, 11}.data(data{4, 11}.data(: , 2) == 2, 7);
                if BB(1) == -1
                   disp('root is dendrite');
                end
                
                % find the same element of two sequence
                C = intersect(AA, BB);
                D = data{4, 11}.data(data{4, 11}.data(: , 2) == 1, 1);
                atrisk = 0;
                for cx = 1: length(C)
                    same = nnz(D ==C(cx));
                    atrisk = atrisk + same;
                end
                %}
                
                
            elseif structureFlag(cx, nx) == 3
                %%% basal dendrite    
                AA = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 1, 7);
                BB = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 3, 7);
                % find the same element of two sequence
                C = intersect(AA, BB );
                
                atrisk(nx, 1) = length(C );
                
            elseif structureFlag(cx, nx) == 4
                %%% apical dendrite   
                AA = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 1, 7);
                BB = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 4, 7);
                % find the same element of two sequence
                C = intersect(AA, BB );
                
                atrisk(nx, 1) = length(C );

                 
            else
                %%% apical and basal dendrite 
                AA = data{cx, nx}.data(data{cx, nx}.data(: , 2) == 1, 7);
                BB = data{cx, nx}.data(data{cx, nx}.data(: , 2) == (tx + 2), 7);
                % find the same element of two sequence
                C = intersect(AA, BB );
                
                atrisk(nx, 1) = length(C );
           
            end % end of if
            
            if structureFlag(cx, nx) == 5
                
                Br{tx}= Tm{cx, nx}.length(Tm{cx, nx}.termFlag == 2 & Tm{cx, nx}.structureID == (tx + 2));
                Te{tx} = Tm{cx, nx}.length(Tm{cx, nx}.termFlag == 0 & Tm{cx, nx}.structureID == (tx + 2));
            
            else
                
                Br{tx}= Tm{cx, nx}.length(Tm{cx, nx}.termFlag == 2 & Tm{cx, nx}.structureID == 2);
                Te{tx} = Tm{cx, nx}.length(Tm{cx, nx}.termFlag == 0 & Tm{cx, nx}.structureID == 2);
            
            end
        
            for R = 1: size(lenRange , 2) - 1
                
                %%% For Branching
                
                inRangeBr = (Br{tx} > lenRange(R) & Br{tx} <= lenRange(R+1));
        
                numOfBr{nx, tx}.data(R) = sum( inRangeBr); 
                
                  
                %%% For Terminal 
           
                inRangeTe = (Te{tx} > lenRange(R) & Te{tx} <= lenRange(R+1));
        
                numOfTe{nx, tx}.data(R) = sum( inRangeTe); 
            
                %%% For atrisk 
                
                atrisk(nx, R+1) = atrisk(nx, R) - sum( inRangeTe) + sum( inRangeBr); 
            
                AtRisk{nx, tx}.data(R) = atrisk(nx, R);
            

                %%% Warning Message %%%
                
                if atrisk(nx, R+1) <= 0
                
                    atrisk(nx, R+1) = 0;
            
                end
                
                %%%%%%%%%%%%%%
            
             end % end of R
 
            
         end %end of tx
         
    end % end of nx   
      
    numOfBrOfType{cx} = numOfBr;
    numOfTeOfType{cx} = numOfTe;
    numOfAtriskOfType{cx} = AtRisk;
    
end %end of cx

 %%% time stop
fprintf('Overall time of part3= %f sec\n', toc(timer3));

%% Calculate the KB and KT 
%%% time Start
timer4 = tic;
groupOfBr = {};
groupOfTe = {};
groupOfAtRisk = {};
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cx = 1: numberOfClasses

    BrMatrix = zeros(numberOfNeuronsInClass(cx), size(lenRange, 2) -1);
    TeMatrix = zeros(numberOfNeuronsInClass(cx), size(lenRange, 2) -1);
    AtRiskMatrix = zeros(numberOfNeuronsInClass(cx), size(lenRange, 2) -1);
    
     dendriteBr = numOfBrOfType{cx};
     dendriteTe =  numOfTeOfType{cx};
     dendriteAtRisk = numOfAtriskOfType{cx};
    
    
    for nx = 1:  numberOfNeuronsInClass(cx)
   
         if structureFlag(cx, nx) == 5 
    
            numberOfType = 2;
    
        else  
            numberOfType = 1;

         end
        
        for tx = 1: numberOfType
             %{        
            BrMatrix(nx, :) = dendriteBr{nx, tx}.data;
            TeMatrix(nx, :) =  dendriteTe{nx, tx}.data;
            AtRiskMatrix(nx, :) =  dendriteAtRisk{nx, tx}.data; 
            %}
            groupOfBr{cx, tx}(nx, :) =dendriteBr{nx, tx}.data;
            groupOfTe{cx, tx}(nx, :) = dendriteTe{nx, tx}.data;
            groupOfAtRisk{cx, tx}(nx, :) = dendriteAtRisk{nx, tx}.data; 
            
        end % end of tx
        
    end % end of nx
    
end % end of cx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cx =1 : numberOfClasses
    
    for tx =1: size(groupOfAtRisk, 2)
        
        groupOfKB{cx, tx} =  groupOfBr{cx, tx}(:, 2:end) ./  (groupOfAtRisk{cx, tx}(:, 1: end-1) .* t_section);
        groupOfKT{cx, tx} = groupOfTe{cx, tx}(:, 2:end) ./ (groupOfAtRisk{cx, tx}(:, 1: end-1) .* t_section);

    end
end

 %%% time stop
fprintf('Overall time of part4= %f sec\n', toc(timer4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

bins = t_section*( 1: size(lenRange , 2)-2) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save mat
save( strcat( matFilePath, 'dendriticLengthDistributionData_TmRoot_', lab_name, '_', date,'.mat' ), 'maxDendriticBranchLength', ...
    'deltaLen', 'numberOfLengthValues', 'lenValues', 'numBrNeu', 'numTeNeu', 'numBrWC', 'numTeWC', 'numberOfSegmentsForNeuron', 't_section', 'bins', ...
    'lenRange', 'date', 'lab_name', 'numOfBr', 'numOfTe', 'AtRisk', 'groupOfKB', 'groupOfKT', 'numberOfClasses','numberOfNeuronsInClass', ...
    'groupOfAtRisk', 'numOfBrOfType', 'numOfTeOfType', 'numOfAtriskOfType', 'structureFlag', 'neuron_names', 'numberOfType', 'Tm', 'groupOfBr', 'groupOfTe')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Coding ends here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
