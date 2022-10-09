function wellData = WellDataProcess(wellSols,nstep,dt)
wellData = struct();
for t=1:nstep
    [qw,qw_inj]  = wellSols{t,1}.qWs;
    [qo,~]  = wellSols{t,1}.qWs;
    [SO4,~] = wellSols{t,1}.qW_SO4;
    [Ca,~]  = wellSols{t,1}.qW_Ca;
    [Sr,~]  = wellSols{t,1}.qW_Sr;
    [Ba,~]  = wellSols{t,1}.qW_Ba;
    [CO3,~] = wellSols{t,1}.qW_CO3;
    [Na,~] = wellSols{t,1}.qW_Na;
    [Cl,~] = wellSols{t,1}.qW_Cl;
    [Mg,~] = wellSols{t,1}.qW_Mg;
    wellData.qOs(t,1) = convertTo(qo,stb/day);
    wellData.qWs(t,1) = convertTo(qw,stb/day);
    wellData.qWs_inj(t,1) = convertTo(qw_inj,stb/day);
    if t==1
        [~,a] = wellSols{t,1}.qWs;
        [b,~] = wellSols{t,1}.qOs;
        a = a*dt;
        b = -b*dt;
        wellData.qw_inj_Total(t,1) = convertTo(a,stb);
        wellData.qo_prod_Total(t,1) = convertTo(b,stb);
    else
        [~,a] = wellSols{t,1}.qWs;
        [b,~] = wellSols{t,1}.qOs;
        a = a*dt;
        a = wellData.qw_inj_Total(t-1,1) + convertTo(a,stb);
        b = -b*dt;
        b = wellData.qo_prod_Total(t-1,1) + convertTo(b,stb);
        wellData.qw_inj_Total(t,1) = a;
        wellData.qo_prod_Total(t,1) = b;
    end
    x_SO4 = SO4/qw;
    x_Ca  = Ca/qw;
    x_Sr  = Sr/qw;
    x_Ba  = Ba/qw;
    x_Na  = Na/qw;
    x_Mg  = Mg/qw;
    x_CO3 = CO3/qw;
    x_Cl  = Cl/qw;
    sigma = x_SO4 + x_Ca + x_Sr + x_Ba + x_Na + x_Mg + x_CO3 + x_Cl;
    wellData.SO4(t,1) = x_SO4/(sigma+1)*1e6;
    wellData.Ca(t,1)  = x_Ca/(sigma+1)*1e6;
    wellData.Sr(t,1)  = x_Sr/(sigma+1)*1e6;
    wellData.Ba(t,1)  = x_Ba/(sigma+1)*1e6;
    wellData.Na(t,1)  = x_Na/(sigma+1)*1e6;
    wellData.Mg(t,1)  = x_Mg/(sigma+1)*1e6;
    wellData.CO3(t,1) = x_CO3/(sigma+1)*1e6;
    wellData.Cl(t,1)  = x_Cl/(sigma+1)*1e6;
end
    