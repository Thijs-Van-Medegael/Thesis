%% Dataset

maxCap = 50;
alphaPeak = 0.3;
alpha = 0.2;
beta = 0.8;

EVnumber = 10;
% EV = zeros(EVnumber, 7);
% for i = 1:EVnumber
%     %TODO: Stochastic battery capacity? 
%     EV(i,1) = unifrnd(40,80); % Battery capacity of the car
%     EV(i,2) = unifrnd(Smin,Smax-0.1); % 1st column: initial SoC
%     EV(i,3) = unifrnd(max(EVs(i,2),0.6),Smax); % 2nd column: Sobj
%     EV(i,4) = EVs(i,1)*(EVs(i,3)-EVs(i,2)); % 3rd column: energy to get to Sobj
%     EV(i,5) = max(offPeak,normrnd(.44,0.1)); % price for which EV will defer charging
%     EV(i,6) = round(normrnd(17,3)); % arrival time (16th index equals 17h)
%     EV(i,7) = min(round(normrnd(8,2)),max(EVs(:,6))); % departure time
% end

initCharge = [7.5; 6; 7; 5.6; 6.7; 7.5; 8.1; 9; 7.2; 7.5];
finalCharge = [67.5; 72; 67.5; 63; 58.6; 67.5; 58.5; 72; 63; 67.5];
minCharge = 5*ones(EVnumber,1);
maxCharge = [75; 80; 75; 70; 65; 75; 65; 80; 70; 75];
rate = [15; 15; 10; 8; 8; 8; 8; 10; 8; 15];
tarr = [17; 18; 19; 17; 18; 12; 11; 12; 13; 14];
tdep = [22; 23; 6; 7; 8; 18; 23; 6; 7; 20];

EV = [initCharge finalCharge minCharge maxCharge rate tarr tdep];
beginPlot = min(EV(:,6));

%% Dumb Charging

demandDumb = zeros(EVnumber,24);

for i=1:EVnumber
    energy = EV(i,2) - EV(i,1);
    demandCurr = zeros(1,24);
    tarr = EV(i,6);
    tdep = EV(i,7);
    t = 1;
    rateCharging = EV(i,5);
    while t ~= (tdep + 24 - tarr)
        if energy > rateCharging
            demandCurr(t) = rateCharging;
            energy = energy - rateCharging;
        else 
            demandCurr(t) = energy;
            break
        end
        t = t+1;
    end
    demandDumb(i,:) = circshift(demandCurr,tarr);
end

%% Plotting Dumb Charging
totalDumb = sum(demandDumb, 1);
dumbShifted = circshift(totalDumb,-beginPlot);

t1 = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
xLabels = [beginPlot mod(beginPlot+5,24) mod(beginPlot+10,24) mod(beginPlot+15,24) mod(beginPlot+20,24)];
figure(1)
plot(t1,dumbShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
title('Dumb charging profile')
xlabel('Time (h)')
ylabel('Load (kWh)')
ylim([0 maxCap+30]);
hold on 
ax = gca; 
ax.FontSize = 12;
yline(maxCap, '-r')
hold off

figure(2)
for i=1:EVnumber
    plot(t1,circshift(demandDumb(i,:),-beginPlot))
    hold on 
end
plot(t1,dumbShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
title('Dumb charging profile')
xlabel('Time (h)')
ylabel('Load (kWh)')
ylim([0 maxCap+30]);
hold on 
ax = gca; 
ax.FontSize = 12;
yline(maxCap, '-r')
legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9','EV10','Total Demand','Line capacity')
hold off

%% Game theory implementation
ToUtariff = [alpha*ones(1,7) alphaPeak*ones(1,10) alpha*ones(1,7)];
rho = 0.1;
eta = 0.01;
crit = 0.01;

% schedule = zeros(EVnumber,24);
previousIteration = demandDumb;
previousTotal = totalDumb;
% price = ToUtariff.*totalDumb + beta;
while norm(previousIteration-schedule) > crit
   
    previousIteration = schedule;
    schedule = zeros(EVnumber, 24);


    for i = 1:EVnumber
        tarr = EV(i,6);
        tdep = EV(i,7);
%         if tdep > tarr 
%             duration = tdep-tarr;
%         else 
%             duration = 24 - tarr + tdep;
%         end
        energy = EV(i,2) - EV(i,1);
        OtherEVs = previousTotal - previousIteration(i,:); % Strategies of other EVs
        Zero = [];
        if tdep > tarr 
            for j=1:24
                if j <= tarr || j > tdep
                    Zero(end+1) = j;
                end
            end
        else 
            for j=1:24
                if j > tdep && j <= tarr
                    Zero(end+1) = j;
                end
            end
        end
        
        equalityMatrix = zeros(1,24);
        equalityMatrix(Zero) = 1;
        equalityMatrix = [ones(1,24);
                          diag(equalityMatrix)];
        equalityVector = [energy;
                          zeros(24,1)];
        cost = @(x)(energy/sum(previousTotal))*sum((ToUtariff.*(OtherEVs+x) + beta).*(OtherEVs + x)) + (rho/2)*norm(previousIteration(i,:)-x)^2;
        x = fmincon(cost, previousIteration(i,:), [], [], equalityMatrix, equalityVector, zeros(1,24), EV(i,5)*ones(1,24));
        x = round(x*100)/100;
        schedule(i,:) = x;
    end
end

%% Plotting Game Theory

totalSchedule = sum(schedule, 1);
gameShifted = circshift(totalSchedule,-beginPlot);

t1 = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
xLabels = [beginPlot mod(beginPlot+5,24) mod(beginPlot+10,24) mod(beginPlot+15,24) mod(beginPlot+20,24)];
figure(3)
plot(t1,gameShifted,'-om','MarkerFaceColor','b')
xticklabels(xLabels)
title('Game Theory profile')
xlabel('Time (h)')
ylabel('Load (kWh)')
ylim([0 maxCap+30]);
hold on 
ax = gca; 
ax.FontSize = 12;
yline(maxCap, '-r')
hold off

% figure(4)
% for k=1:EVnumber
%     plot(t1,circshift(gameShifted(k,:),-beginPlot))
%     hold on 
% end
% plot(t1,gameShifted,'-om','MarkerFaceColor','b')
% xticklabels(xLabels)
% title('Game Theory charging profile')
% xlabel('Time (h)')
% ylabel('Load (kWh)')
% ylim([0 maxCap+30]);
% hold on 
% ax = gca; 
% ax.FontSize = 12;
% yline(maxCap, '-r')
% legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9','EV10','Total Demand','Line capacity')
% hold off
