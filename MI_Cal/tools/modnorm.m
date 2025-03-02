function scale = modnorm(constellation, pwrType, pwr)
    if(strcmp(pwrType, 'peakpwr'))
        constpwr = max(abs(constellation).^2);
    end
    if (strcmp(pwrType, 'avgpwr'))
        constpwr = mean(abs(constellation).^2);
    end
    scale = sqrt(pwr/constpwr);
end