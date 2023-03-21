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

% EVs = [0.815562772330580	0.886705460945767	4.26856131691120	0.375580787466621	14	4;
% 0.404244850617017	0.837916056899460	26.0202723769465	0.370670000000000	18	10;
% 0.183177006818200	0.891696122998194	42.5111469707996	0.496548467187182	17	12;
% 0.177154943473722	0.839642767346171	39.7492694323469	0.370670000000000	15	9;
% 0.773758521872488	0.899751719630036	7.55959186545284	0.570484977907667	14	9;
% 0.239216709837526	0.715815927511889	28.5959530604618	0.370670000000000	19	12;
% 0.699207468036466	0.926825415525194	13.6570768493237	0.445898042631383	13	6;
% 0.329054574259535	0.985253764092006	39.3719513899483	0.391963856392989	17	5;
% 0.447233035641351	0.656884788987650	12.5791052007780	0.370670000000000	19	11;
% 0.823032725167332	0.918661546228361	5.73772926366172	0.547326114248551	19
% 3]; DOEST WORK BECAUSE TO MUCH EVS HAVE SIMILAR THRESHOLD PRICE
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
    t = 1;
    while t ~= 25
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
    t = mod(EVs(i,5),24);
    while t ~= EVs(i,5)-1
        if t == 24
            t = 0;
        end
        if priceElec(t+1) <= EVs(i,4)
            if energy >= Emax
                demandToU(i,t+1) = Emax;
                energy = energy - Emax;
            else 
                demandToU(i,t+1) = energy;
                break
            end
        end
        t=t+1;
    end
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
        t = EVs(i,5); % arrival time 
        while t ~= EVs(i,5)-1 
            if t == 24
                t = 0;
            end
            if price(t+1) <= EVs(i,4)
                if energy >= Emax
                    demandIteration(i,t+1) = Emax;
                    energy = energy - Emax;
                else 
                    demandIteration(i,t+1) = energy;
                    break
                end
            end
            t=t+1;
        end
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

if logical(sum(totalFinal) ~= sum(totalDumb))
       disp('ITERATIVE FAILED');
end



%% Discomfort

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
        weight = [price2(tarr:end) price2(1:tdep-1)];
     
        
        schedule = linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2))); % works as it should
        schedule = [schedule(24-tarr:end) (zeros(1,24-size(weight,2))) (schedule(1:24-tarr))];
        demandIteration2(i,:) = schedule;
        
    end
    totalIteration2 = sum(demandIteration2, 1);
end


totalFinal2 = totalIteration;


