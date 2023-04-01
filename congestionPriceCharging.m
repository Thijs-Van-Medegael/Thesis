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
price = zeros(1,24);
w = 0.4; % weighting factor for congestion pricing 

totalIteration = cap+1;
while max(totalIteration) > cap
    demandIteration = zeros(EVnumber,size(t1,2));
    congestionPrice = arrayfun(@(x)w*max(x-cap,0)/cap,totalIteration);
    price = price + congestionPrice;

    for i=1:EVnumber
        energy = EVs(i,4);
        tarr = EVs(i,6);
        tdep = EVs(i,7);
        priceCurrIter = circshift(price,-tarr);
        demandIterCurr = zeros(1,24);
        t = 1;

        while t ~= (tdep + 24 - tarr)
%             if t == 24
%                 t = 0;
%             end
            if priceCurrIter(t) <= EVs(i,5)
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
    totalIteration = sum(demandIteration, 1) + baseload;
end


totalFinal = totalIteration;
CongShifted = circshift(totalFinal,-beginPlot);
xLabels = [beginPlot mod(beginPlot+5,24) mod(beginPlot+10,24) mod(beginPlot+15,24) mod(beginPlot+20,24)];

figure(1)
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