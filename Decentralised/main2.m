clear
clc
%%%HOW TO ADD BASELOAD?
%%%CHARGING EARLIER = BETTER?

Smin = 0.1;
Smax = 1;
capBat = 70; % in kWh (Tesla Model Y)
Emax = 15; % in kW
cap = 200; % Grid (Line) capacity

offPeak = 0.37067; % Prices per kWh
onPeak = 0.48695;

EVnumber = 10;
EVs = zeros(EVnumber,5);
for i = 1:EVnumber
    %TODO: Stochastic battery capacity? 
    EVs(i,1) = unifrnd(40,80); % Battery capacity of the car
    EVs(i,2) = unifrnd(Smin,Smax-0.3); % 1st column: initial SoC
    EVs(i,3) = unifrnd(max(EVs(i,2),0.6),Smax); % 2nd column: Sobj
    EVs(i,4) = EVs(i,1)*(EVs(i,3)-EVs(i,2)); % 3rd column: energy to get to Sobj
    EVs(i,5) = max(offPeak,normrnd(.44,0.1)); % price for which EV will defer charging
    EVs(i,6) = round(normrnd(17,3)); % arrival time (16th index equals 17h)
    EVs(i,7) = min(round(normrnd(8,2)),max(EVs(:,6))); % departure time
end

% EVs = [52.3747696969100	0.518423838753594	0.730119596997771	11.0875165838515	0.585274196098984	21	8;
% 50.5185613492063	0.644452961843560	0.727527645661536	4.19681351102432	0.427498322405135	15	8;
% 79.6681488384782	0.704176321342264	0.994218081658574	23.1070901302540	0.370670000000000	18	6;
% 70.2750188393849	0.581583813020408	0.942867483648884	25.3892167597782	0.549922948349201	19	7;
% 40.0136584887256	0.532702658210435	0.683092248798424	6.01763771804705	0.370670000000000	16	5;
% 69.9013464695321	0.698807165914634	0.862445061922510	11.4385092643919	0.402159224149845	17	11;
% 75.7133327727201	0.385202893701604	0.818560757196755	32.8109681284834	0.370670000000000	18	9;
% 69.8349960616260	0.200428985083856	0.928957602703608	50.8767931422572	0.370670000000000	16	10;
% 71.2549600110385	0.393828732105095	0.897947142496669	35.9209371732799	0.603449867890992	15	5;
% 49.0027152851165	0.380011135200613	0.714833845258471	16.4072219319563
% 0.497668116986182	11	8]; DOES NOT WORK

beginPlot = min(EVs(:,6));
t1 = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
% baseload = [92.684 82.153 79.172 73.166 74.093 79.684 109.067 142.416 158.256 156.271 152.508 162.324 152.579 139.126 136.246 142.246 155.712 230.033 216.728 204.236 222.338 186.18 161.535 130.007];
rawdata = readmatrix('Data/fulldata.xlsx');
% fixedLoad5min = sum(rawdata(:,1:10),2);
fixedLoad = zeros(EVnumber,24);
for i=1:24
    for j=1:EVnumber
        avg = reshape(rawdata(:,j),288,8); % Reshape vector to get the profile for one day in a column, for eight days. 
        avg = mean(avg,2); % Take the average of the eigth days. 
        fixedLoad(j,i) = sum(avg(1+12*(i-1):12*i)); % Rework 5 min resolution to 1 hour resolution.
    end
end


residentBehaviour = mean(fixedLoad([1 10],:),1); % first and last row correspond to a regular resident
totalFixed = residentBehaviour*EVnumber;

%% Dumb Charging
demandDumb = zeros(EVnumber,size(t1,2)); % For each EV, we generate a dumb charging load profile (time window is 16h-6h)

for i=1:EVnumber
    energy = EVs(i,4);
    demandCurr = zeros(1,24);
    tarr = EVs(i,6);
    tdep = EVs(i,7);
    t = 1;
    while t ~= (tdep + 24 - tarr)
        if energy > Emax
            demandCurr(t) = Emax;
            energy = energy - Emax;
        else 
            demandCurr(t) = energy;
            break
        end
        t = t+1;
    end
    demandDumb(i,:) = circshift(demandCurr,tarr);
end

totalDumb = sum(demandDumb, 1) + totalFixed;

dumbShifted = circshift(totalDumb,-beginPlot);
% timeShifted = circshift(t1,-beginPlot);
% timeSeriesDumb = timeseries(dumbShifted, timeShifted);
xLabels = [beginPlot mod(beginPlot+5,24) mod(beginPlot+10,24) mod(beginPlot+15,24) mod(beginPlot+20,24)];
figure(1)
plot(t1,dumbShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
title('Dumb charging profile')
xlabel('Time (h)')
ylabel('Load (kW)')
ylim([0 cap+30]);
hold on 
ax = gca; 
ax.FontSize = 12;
yline(cap, '-r')
hold off
 
%% DSO request

requestDSO = totalDumb;
requestDSO(requestDSO > cap) = cap; % requestDSO(requestDSO<=cap) = 0;
% requestDSO(requestDSO>0) = cap;
congestionHours = find(requestDSO==max(requestDSO));
% requestDSO = requestDSO(congestionHours);
% requestDSO = cap*ones(1,24);

%%% INITIAL PRICE SETTING

alpha = 1; % learning rate 
% sumEVs = totalDumb;
% sumEVs(sumEVs<=cap) = 0;
% sumEVs(congestionHours) = totalDumb(congestionHours);

normcst = 1/max(totalDumb-requestDSO);
currPrice = 10^-7*sum(fixedLoad,1).^3; %Price proportional to fixed load profile
newPrice = currPrice + alpha*normcst*max(totalDumb-requestDSO,0); 
epsilon = 0.1;
forecastedSchedule = zeros(EVnumber,24);
%%% ITERATIONS
while norm(newPrice-currPrice) > epsilon
    currPrice = newPrice;
%     forecastedSchedule = zeros(EVnumber,24);
    
    % EV MINIMIZING THEIR COSTS
    for i=1:EVnumber 
        energy = EVs(i,4);
        tarr = EVs(i,6);
        tdep = EVs(i,7);
        
        if tarr > 24
            chargingTime = tdep-tarr;
        else 
            chargingTime = 24 - tarr + tdep;
        end
        
        price = circshift(currPrice,-tarr); price = price(1:chargingTime);
        request = circshift(requestDSO/EVnumber,-tarr); request = request(1:chargingTime);

        objective = @(x) price*transpose(x-request);
        startOpt = circshift(forecastedSchedule(i,:),-tarr); startOpt = startOpt(1:chargingTime);
%         equality = ones(1,24); % sum of energy charged throughout chargingTime is equal to energy needed to get to objective SoC

        profile = fmincon(objective,startOpt, [], [], ones(1,chargingTime), energy, zeros(1,chargingTime), Emax*(ones(1,chargingTime))); % Linear optimization to find EV charging profile
        profile = [profile zeros(1, 24-length(profile))];
        
        forecastedSchedule(i,:) = circshift(profile,tarr);
    end
    
    %PRICE SETTING BY THIRD PARTY (DSO/AGGREGATOR)
    total = sum(forecastedSchedule, 1);
    priceUpdate = max(sum(forecastedSchedule,1)+totalFixed-cap*ones(1,24),0);
    

%     newPrice = currPrice + alpha*normcst*max(sum(forecastedSchedule,1)-requestDSO,0); 
    newPrice = currPrice + alpha*normcst*priceUpdate;
end

totalGame = sum(forecastedSchedule,1)+totalFixed;
GameShifted = circshift(totalGame,-beginPlot);

figure(2)
plot(t1,GameShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
% xlim([0 23])
title('Game Theory charging profile')
xlabel('Time (h)')
ylabel('Load (kW)')
ax = gca; 
ax.FontSize = 12;
ylim([0 cap+30]);
hold on 

yline(cap, '-r')
hold off

