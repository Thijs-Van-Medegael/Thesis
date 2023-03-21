%% Parameters
clear
close all
Smin = 0.1;
Smax = 1;
capBat = 60; % in kWh (Tesla Model Y)
Emax = 10; % in kW
cap = 30; % Grid (Transformer) capacity


offPeak = 0.37067; % Prices per kWh
onPeak = 0.48695;

% offPeak = 0.1; % Prices per kWh
% onPeak = 0.5;

EVnumber = 10;

EVs = zeros(EVnumber,5);

% falpha = [/39 10/39 6/39 10/39 3/39 2/39 0/39 0/39 0/39 0/39 0/39 0/39 0/39 0/39];
t1 = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
% t2 = [16 17 18 19 20 21 22 23 0 1 2 4 5 6];
% n = EVnumber;
% 
% c = cumsum([0,falpha(:).']);
% c = c/c(end); % make sur the cumulative is 1 
% [~,i] = histc(rand(1,n),c);
% r = t1(i); % map to v values
% 
% EVs(:,5) = r'; % When the charging task starts


for i = 1:EVnumber
    %TODO: Stochastic battery capacity? 
    EVs(i,1) = unifrnd(Smin,Smax-0.1); % 1st column: initial SoC
    EVs(i,2) = unifrnd(max(EVs(i,1),0.6),Smax); % 2nd column: Sobj
    EVs(i,3) = capBat*(EVs(i,2)-EVs(i,1)); % 3rd column: energy to get to Sobj
    EVs(i,4) = max(offPeak,normrnd(.44,0.1)); % price for which EV will defer charging
    EVs(i,5) = round(normrnd(17,3)); % arrival time (16th index equals 17h)
    EVs(i,6) = min(round(normrnd(8,2)),max(EVs(:,5))); % departure time
end


% EVs = [0.751778949114543	0.976615575619027	13.4901975902690	0.370670000000000	20	9;
% 0.178032323999528	0.711399287546819	32.0020178128375	0.474262446653865	28	14;
% 0.226090465342039	0.988237112704246	45.7287988417325	0.512540422494611	17	9;
% 0.213509070901772	0.768704513050510	33.3117265289243	0.588969760778547	21	11;
% 0.624592559325270	0.637998989181005	0.804385791344109	0.511723865132884	22	9;
% 0.706192104462667	0.924530291027913	13.1002911939148	0.409655907521398	18	6;
% 0.664836870415687	0.675506066831120	0.640151784925978	0.370670000000000	15	2;
% 0.758766262661834	0.926382368191822	10.0569663317993	0.370670000000000	21	5;
% 0.450995487725119	0.752623382837203	18.0976737067251	0.471920673916550	18	6;
% 0.491811516630585	0.778234480284360	17.1853778192265	0.502770728752873	20	10];

% EVs = [0.602012037521139	0.814792820673305	12.7668469891299	0.370670000000000	22	6;
% 0.420023772495711	0.807089405757887	23.2239379957306	0.370670000000000	14	5;
% 0.179024028698803	0.956110182218791	46.6251692111993	0.370670000000000	20	8;
% 0.190973152835646	0.796178262784808	36.3123065969497	0.480474457080366	14	10;
% 0.622415018062455	0.796531583736583	10.4469939404477	0.501434740533673	13	9;
% 0.363180244051477	0.976520212459366	36.8003981044734	0.370670000000000	23	8;
% 0.551663887446005	0.997486452640385	26.7493539116628	0.550295034306817	14	6;
% 0.135259033455307	0.925177587869021	47.3951132648228	0.426430191232172	16	8;
% 0.813719173513560	0.889334405932335	4.53691394512651	0.470928218227206	10	7;
% 0.215873470769043	0.697952211897265	28.9247244676933	0.370670000000000	16	7];

beginPlot = min(EVs(:,5));

%% Dumb Charging

demandDumb = zeros(EVnumber,size(t1,2)); % For each EV, we generate a dumb charging load profile (time window is 16h-6h)
% for i=1:EVnumber
%     energy = EVs(i,3);
% %     t = mod(EVs(i,5),24)+1;
%     t = mod(EVs(i,5),24);
%     while t ~= mod(EVs(i,5),24)-1
% %     while t ~= EVs(i,5)
%         if t == 24
%             t = 0;
%         end
%         if energy > Emax
%             demandDumb(i,t+1) = Emax;
%             energy = energy - Emax;
%         else 
%             demandDumb(i,t+1) = energy;
%             break
%         end
%         t = t+1;
%     end
% end
% 
% totalDumb = sum(demandDumb, 1);

for i=1:EVnumber
    energy = EVs(i,3);
    demandCurr = zeros(1,24);
    tarr = EVs(i,5);
    tdep = EVs(i,6);
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

totalDumb = sum(demandDumb, 1);


dumbShifted = circshift(totalDumb,-beginPlot);
% timeShifted = circshift(t1,-beginPlot);
% timeSeriesDumb = timeseries(dumbShifted, timeShifted);
xLabels = [beginPlot mod(beginPlot+5,24) mod(beginPlot+10,24) mod(beginPlot+15,24) mod(beginPlot+20,24)];
figure(1)
plot(t1,dumbShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
title('Dumb charging profile')
xlabel('Time (h)')
ylabel('Load (kWh)')
ylim([0 cap+30]);
hold on 
ax = gca; 
ax.FontSize = 12;
yline(cap, '-r')
hold off
% plotResult(totalDumb,cap)

%% ToU
priceElec = zeros(1,24);
priceElec(8:22) = onPeak; % Day tarrif is from 7am until 22pm
priceElec(priceElec == 0) = offPeak; % Night tariff is from 22pm until 7am

demandToU = zeros(EVnumber,size(t1,2));
for i=1:EVnumber
    energy = EVs(i,3);
    tarr = EVs(i,5);
    tdep = EVs(i,6);
    priceElecCurr = circshift(priceElec,-tarr);
    demandToUCurr = zeros(1,24);
    t = 1;
    while t ~= (tdep + 24 - tarr)
        if priceElecCurr(t) <= EVs(i,4)
            if energy >= Emax
                demandToUCurr(t) = Emax;
                energy = energy - Emax;
            else 
                demandToUCurr(t) = energy;
                break
            end
        end
        t=t+1;
    end
    demandToU(i,:) = circshift(demandToUCurr,tarr);
end

totalToU = sum(demandToU, 1);
ToUShifted = circshift(totalToU,-beginPlot);

figure(2)
plot(t1,ToUShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
% xlim([0 23])
title('ToU charging profile')
xlabel('Time (h)')
ylabel('Load (kW)')
ax = gca; 
ax.FontSize = 12;
ylim([0 cap+30]);
hold on 

yline(cap, '-r')
hold off

% plotResult(totalToU,cap)


%% Congestion Pricing (Iterative)

price = zeros(1,24); % index 1 corresponds to 0h-1h etc. 
totalIteration = totalDumb;
w = 0.1; % weighting factor for congestion pricing 


while max(totalIteration) > cap
    demandIteration = zeros(EVnumber,size(t1,2));
    congestionPrice = arrayfun(@(x)w*max(x-cap,0)/cap,totalIteration);
    price = price + congestionPrice;

    for i=1:EVnumber
        energy = EVs(i,3);
        tarr = EVs(i,5);
        tdep = EVs(i,6);
        priceCurrIter = circshift(price,-tarr);
        demandIterCurr = zeros(1,24);
        t = 1;

        while t ~= (tdep + 24 - tarr)
%             if t == 24
%                 t = 0;
%             end
            if priceCurrIter(t) <= EVs(i,4)
                if energy >= Emax
                    demandIterCurr(t) = Emax;
                    energy = energy - Emax;
                else 
                    demandIterCurr(t) = energy;
                    break
                end
            end
            t=t+1;
        end
        demandIteration(i,:) = circshift(demandIterCurr,tarr);
    end
    totalIteration = sum(demandIteration, 1);
end


totalFinal = totalIteration;
CongShifted = circshift(totalFinal,-beginPlot);

figure(3)
plot(t1,CongShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
% xlim([0 23])
title('Congestion Pricing charging profile')
xlabel('Time (h)')
ylabel('Load (kW)')
ax = gca; 
ax.FontSize = 12;
ylim([0 cap+30]);
hold on 

yline(cap, '-r')
hold off
% plotResult(totalFinal,cap)

% if sum(totalFinal) - sum(totalDumb) > 1^(-10)
%        disp('ITERATIVE FAILED');
% end



%% Linear programming
% Linprog takes last timeslots if costs are the same. More accurate cost
% function needed including discomfort

price2 = zeros(1,24); % index 1 corresponds to 0h-1h etc. 
totalIteration2 = totalDumb; 
w = 0.1; % weighting factor for congestion pricing 


while max(totalIteration2) > cap
    demandIteration2 = zeros(EVnumber,size(t1,2));
    congestionPrice2 = arrayfun(@(x)w*max(x-cap,0)/cap,totalIteration2); 
    price2 = price2 + congestionPrice2;  

    for i=1:EVnumber
        % Charging each hour is < Emax and total charging is equal to (SoCobj-SoCinit)*capBat
%         demandCurrent = fmincon((@(x) price2*transpose(x)+sum(demandDumb(i,:)-x)),demandDumb(i,:),eye(24),Emax*ones(24,1),ones(1,24),EVs(i,3));

        tdep = EVs(i,6);
        tarr = EVs(i,5);
        energy = EVs(i,3);
        if tarr >= 24
            weight = [price2(mod(tarr,24)+1:tdep)];
            schedule = transpose(linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2)))); % works as it should
            schedule = [zeros(1,mod(tarr,24)) schedule zeros(1,24-tdep)];
        else
            weight = [price2(tarr+1:end) price2(1:tdep)];
            schedule = transpose(linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2)))); % works as it should
            schedule = [schedule(24-tarr+1:end) zeros(1,24-size(weight,2)) schedule(1:24-tarr)];
        end
        
%         schedule = transpose(linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2)))); % works as it should
%         schedule = [schedule(24-tarr+1:end) zeros(1,24-size(weight,2)) schedule(1:24-tarr)];
        demandIteration2(i,:) = schedule;
        
    end
    totalIteration2 = sum(demandIteration2, 1);
end


totalFinal2 = totalIteration;

LinearShifted = circshift(totalFinal2,-beginPlot);

figure(4)
plot(t1,LinearShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
% xlim([0 23])
title('Linear Programming charging profile')
xlabel('Time (h)')
ylabel('Load (kW)')
ax = gca; 
ax.FontSize = 12;
ylim([0 cap+30]);
hold on 

yline(cap, '-r')
hold off
% plotResult(totalFinal,cap)

% if sum(totalFinal2) - sum(totalDumb) > 1^(-10)
%        disp('LINEAR PROGRAMMING FAILED');
% end

