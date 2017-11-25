function showpos(count,framenumber,positions)
    for i = 1:framenumber,
        fprintf('%d %d \n',positions(count,i,1),positions(count,i,2));
    end
end