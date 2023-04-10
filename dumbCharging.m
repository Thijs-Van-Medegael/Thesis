Smin = 0.1;
Smax = 1;
% capBat = 60; % in kWh (Tesla Model Y)
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
demandDumb = zeros(EVnumber,size(t1,2)); % For each EV, we generate a dumb charging load profile

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


dumbShifted = circshift(totalDumb+baseload,-beginPlot);
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

