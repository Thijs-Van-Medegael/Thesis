function output= DSOrequest(dumbCharging,cap)
%     for i=1:len(dumbCharging)
%         if dumbCharging(i) > cap
%             dumbCharging(i) = cap;
%         else
%             dumbCharging(i) = 0;
%         end
%     end

    dumbCharging(dumbCharging>cap) = 0;
    dumbCharging(dumbCharging>0) = cap;
end

