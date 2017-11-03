% A model of the distribution of multi-temperature food items
% Programmed by: Daniel Botz and Abhinaya Rajendran
% The goal is to develop a route for a reefer trailer to deliver food items
% stored at different temperatures while minimizing the time the food 
path(path,'j:/eos/lockers/research/ie/kay/matlog');
disp('This procedure constructs the short route (time-wise)')
disp('with the US Highway Network data.');
Tcoord = readtable('FoodLions.xlsx','Range','P2:Q12',...
    'ReadVariableNames', true);
Tdrive = table2array(Tcoord);
expansionAroundXY = 0.15;
[flXY,IJD,isXY,isIJD] = subgraph(usrdnode('XY'),...
   isinrect(usrdnode('XY'),boundrect(Tdrive,expansionAroundXY)),...
   usrdlink('IJD'));
s = usrdlink(isIJD);
isI = s.Type == 'I';         % Interstate highways
isIR = isI & s.Urban == ' '; % Rural Interstate highways
isIU = isI & ~isIR;          % Urban Interstate highways
isR = s.Urban == ' ' & ~isI; % Rural non-Interstate roads
isU = ~isI & ~isR;           % Urban non-Interstate roads
% Add connector roads from cities to road network
[IJD11,IJD12,IJD22] = addconnector(Tdrive,flXY,IJD);
% Convert road distances to travel times (needs to be after ADDCONNECTOR)
IJT = IJD;
IJT(isIR,3) = IJD(isIR,3)/75;
IJT(isIU,3) = IJD(isIU,3)/65;
IJT(isR,3) = IJD(isR,3)/50;
IJT(isU,3) = IJD(isU,3)/25;
IJT22 = IJD22;               % road to road
IJT22(:,3) = IJT(:,3);
IJT12 = IJD12;               % facility to road
IJT12(:,3) = IJD12(:,3)/20;  % (IJD11 facility to facility arcs ignored)
% ***
n = size(Tdrive,1); % Shortest routes based on time taken to travel.
[T,P] = dijk(list2adj([IJT12; IJT22]),1:n,1:n);
% ***
% The size of each Food Lion order is uniformly distributed between 50 and 250.
orderFL = 50+randi(200, [9,1]);
tempRange = [0.7 0.15 0.15]; % The portion of the orders stored at 
    % refrigerated, frozen, and super-frozen temperatures, respectfully.
itemMatx = round(orderFL * tempRange);
palletNo = zeros(size(itemMatx));
unloadTime = zeros(size(orderFL));
for ord = 1:size(orderFL,1)
    [unloadTime(ord),palletNo(ord,:)] = unlFood(itemMatx(ord,1),...
        itemMatx(ord,2), itemMatx(ord,3), 90);
end
shT = vec2struct('b',1, 'e',2:height(Tcoord), 'q',orderFL * 1/100, 's',10,...
    'tU',unloadTime/60); % structure of the shipments to each store.
reef = vec2struct('b',1, 'e',1, 'Kcu',(53*8*13), 'Kwt',25);
Tminute = readtable('FoodLions.xlsx','Range','D15:N25',...
    'ReadVariableNames', true);
Tarray = table2array(Tminute(:,2:11)); % later algorithms can't accept data
    % from tables as inputs.
% Question: Is the road data from Oak Ridge National Highway Network
% accuracy enough to determine shortest distance?
rTDh = @(rte) rteTC(rte,shT,T,reef);
IJS = pairwisesavings(rTDh,shT);
[r,~] = twoopt(savings(rTDh,shT,IJS,true),rTDh,true);
[TC,Xflg,out] = rTDh(r);
out;
sprintf('Total time taken to complete all deliveries (US Highway Network): %4.1f_minutes',...
    60 * (sum(out.Cost)+sum(out.LU)))
%% Here the shortest route is found with the tspnneighbor function.
disp('The shortest route is determined using the nearest neighbor')
disp('algorithm for the Traveling Salesperson Problem.')
Tmiles = readtable('FoodLions.xlsx','Range','D2:N12',...
    'ReadVariableNames', true);
Tminute = readtable('FoodLions.xlsx','Range','D15:N25',...
    'ReadVariableNames', true);
FoodLion = table2cell(Tminute(:,1));
Tarray = table2array(Tminute(:,2:11)); % later algorithms can't accept data
    % from tables as inputs.
% The size of each Food Lion order is uniformly distributed between 50 and 250.
orderFL = 50+randi(200, [9,1]);
% using list2adj function, then locTC to determine the shortest time route.
tempRange = [0.7 0.15 0.15]; % The portion of the orders stored at 
    % refrigerated, frozen, and super-frozen temperatures, respectfully.
itemMatx = round(orderFL * tempRange);
palletNo = zeros(size(itemMatx));
unloadTime = zeros(size(orderFL));
for ord = 1:size(orderFL,1)
    [unloadTime(ord),palletNo(ord,:)] = unlFood(itemMatx(ord,1),...
        itemMatx(ord,2), itemMatx(ord,3), 90);
    % We assume here that 90 boxes can be stacked on 1 pallet.
end
% Dimension of reefer trailer: (Length × Width × Height): (53’ × 8’ × 13').
% We also assume for now that the weight capacity of the trailer is 25 tons.
reef.Kcu = (53*8*13); reef.Kwt = 25;
% We consider each box to be 20lb and 2 cubic feet.
boxes = [2, 20];
[maxPL, truckC] = maxpayld(orderFL, reef, boxes);
% Next step would be to find the shortest path with respect to time 
% traveled. This is done with the nearest neighbor algorithm for the 
% Traveling Salesperson Problem. 
[loc,TCost] = tspnneighbor(Tarray,1);
locSeq = cell(length(loc)-1,2);% Addresses and distances between stores.
for pt = 1:(length(loc)-1)
    locSeq{pt,1} = FoodLion{loc(pt)};
    locSeq{pt,2} = Tarray(loc(pt),loc(pt+1));
end
sprintf('Total time taken to complete all deliveries (TSPnneighbor): %3d_minutes',...
    TCost+sum(unloadTime)) 