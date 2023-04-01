Smin = 0.1;
Smax = 1;
capBat = 60; % in kWh (Tesla Model Y)
Emax = 10; % in kW
cap = 240; % Grid (Transformer) capacity

offPeak = 0.37067; % Prices per kWh
onPeak = 0.48695;

EVnumber = 10;
EVs = zeros(EVnumber,5);
for i = 1:EVnumber
    %TODO: Stochastic battery capacity? 
    EVs(i,1) = unifrnd(40,80); % Battery capacity of the car
    EVs(i,2) = unifrnd(Smin,Smax-0.1); % 1st column: initial SoC
    EVs(i,3) = unifrnd(max(EVs(i,2),0.6),Smax); % 2nd column: Sobj
    EVs(i,4) = EVs(i,1)*(EVs(i,3)-EVs(i,2)); % 3rd column: energy to get to Sobj
    EVs(i,5) = max(offPeak,normrnd(.44,0.1)); % price for which EV will defer charging
    EVs(i,6) = round(normrnd(17,3)); % arrival time (16th index equals 17h)
    EVs(i,7) = min(round(normrnd(8,2)),max(EVs(:,6))); % departure time
end

beginPlot = min(EVs(:,6));
t1 = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
baseload = [92.684 82.153 79.172 73.166 74.093 79.684 109.067 142.416 158.256 156.271 152.508 162.324 152.579 139.126 136.246 142.246 155.712 230.033 216.728 204.236 222.338 186.18 161.535 130.007];

%%
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

totalDumb = sum(demandDumb, 1);

%%
currPrice = zeros(1,24);
totalSchedule = totalDumb + baseload;

iterationCount = 0;
k = 10^-7;

while max(totalSchedule) > cap
    demandEVs = zeros(EVnumber,size(t1,2));

    for i=1:EVnumber
        % Charging each hour is < Emax and total charging is equal to (SoCobj-SoCinit)*capBat
%         demandCurrent = fmincon((@(x) price2*transpose(x)+sum(demandDumb(i,:)-x)),demandDumb(i,:),eye(24),Emax*ones(24,1),ones(1,24),EVs(i,3));
        tdep = EVs(i,7);
        tarr = EVs(i,6);
        energy = EVs(i,4);
        if tarr >= 24
            weight = [currPrice(mod(tarr,24)+1:tdep)];
            schedule = transpose(linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2)))); % works as it should
            schedule = [zeros(1,mod(tarr,24)) schedule zeros(1,24-tdep)];
        else
            weight = [currPrice(tarr+1:end) currPrice(1:tdep)];
            schedule = transpose(linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2)))); % works as it should
            schedule = [schedule(24-tarr+1:end) zeros(1,24-size(weight,2)) schedule(1:24-tarr)];
        end
        
%         schedule = transpose(linprog(weight, [], [], ones(1,tdep + (24-tarr)), energy, zeros(1,size(weight,2)),Emax*ones(1,size(weight,2)))); % works as it should
%         schedule = [schedule(24-tarr+1:end) zeros(1,24-size(weight,2)) schedule(1:24-tarr)];
        demandEVs(i,:) = schedule;
        
    end

    iterationCount = iterationCount +1;
    totalSchedule = sum(demandEVs, 1) + baseload;

    % Price setting
    currPrice = k*(totalSchedule).^3; 
end


GameShifted = circshift(totalSchedule,-beginPlot);
xLabels = [beginPlot mod(beginPlot+5,24) mod(beginPlot+10,24) mod(beginPlot+15,24) mod(beginPlot+20,24)];

figure(5)
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
