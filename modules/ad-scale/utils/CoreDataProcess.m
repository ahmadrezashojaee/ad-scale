function [CoreData] = CoreDataProcess(states, nstep,NX ,dx, rockProps)
CoreData = struct();
CoreData.Pressure = deal(cell(nstep,1));
CoreData.Sat      = deal(cell(nstep,1));
CoreData.qw       = deal(zeros(nstep,1));
CoreData.qo       = deal(zeros(nstep,1));
CoreData.SO4      = deal(cell(nstep,1));
CoreData.Ca       = deal(cell(nstep,1));
CoreData.Mg       = deal(cell(nstep,1));
CoreData.Na       = deal(cell(nstep,1));
CoreData.Sr       = deal(cell(nstep,1));
CoreData.Ba       = deal(cell(nstep,1));
CoreData.CO3      = deal(cell(nstep,1));
CoreData.Cl       = deal(cell(nstep,1));
CoreData.SO4_out      = deal(zeros(nstep,1));
CoreData.Ca_out       = deal(zeros(nstep,1));
CoreData.Mg_out       = deal(zeros(nstep,1));
CoreData.Na_out       = deal(zeros(nstep,1));
CoreData.Sr_out       = deal(zeros(nstep,1));
CoreData.Ba_out       = deal(zeros(nstep,1));
CoreData.CO3_out      = deal(zeros(nstep,1));
CoreData.Cl_out       = deal(zeros(nstep,1));
CoreData.k_ave       = deal(zeros(nstep,2));
for t=1:nstep
    p=states{t,1};
    CoreData.Pressure{t,1} = p.pressure;
    CoreData.Sat{t,1}     = p.s;
    CoreData.qw(t) = convertTo(p.flux(NX+1,1),centi^3/hour);
    CoreData.qo(t) = convertTo(p.flux(NX+1,2),centi^3/hour);
    x_SO4 = p.c_SO4;
    x_Ca  = p.c_Ca;
    x_Sr  = p.c_Sr;
    x_Ba  = p.c_Ba;
    x_Na  = p.c_Na;
    x_Mg  = p.c_Mg;
    x_CO3 = p.c_CO3;
    x_Cl  = p.c_Cl;
    sigma = x_SO4 + x_Ca + x_Ba + x_Sr + x_Na + x_Cl + x_Mg + x_CO3;
    CoreData.SO4{t,1}     = x_SO4./(1+sigma).*1e6;
    a = CoreData.SO4{t,1};
    CoreData.SO4_out(t) = a(end);
    CoreData.Ca{t,1}      = x_Ca./(1+sigma).*1e6;
    a = CoreData.Ca{t,1};
    CoreData.Ca_out(t) = a(end);
    CoreData.Sr{t,1}      = x_Sr./(1+sigma).*1e6;
    a = CoreData.Sr{t,1};
    CoreData.Sr_out(t) = a(end);
    CoreData.Ba{t,1}      = x_Ba./(1+sigma).*1e6;
    a = CoreData.Ba{t,1};
    CoreData.Ba_out(t) = a(end);
    CoreData.Na{t,1}      = x_Na./(1+sigma).*1e6;
    a = CoreData.Na{t,1};
    CoreData.Na_out(t) = a(end);
    CoreData.Mg{t,1}      = x_Mg./(1+sigma).*1e6;
    a = CoreData.Mg{t,1};
    CoreData.Mg_out(t) = a(end);
    CoreData.CO3{t,1}     = x_CO3./(1+sigma).*1e6;
    a = CoreData.CO3{t,1};
    CoreData.CO3_out(t) = a(end);
    CoreData.Cl{t,1}      = x_Cl./(1+sigma).*1e6;
    a = CoreData.Cl{t,1};
    CoreData.Cl_out(t) = a(end);
    k1=1;
    k2=0;
    for i=1:NX
        perm = rockProps.perm{t,1};
        k1=k1*perm(i);
        k2=k2+dx/NX/perm(i);
    end
    k1=(k1)^(1/NX);
    CoreData.k_ave(t,1) = convertTo(k1,milli*darcy); %Geometric Average
    CoreData.k_ave(t,2) = convertTo(dx/k2,milli*darcy); %Harmonic Average
end