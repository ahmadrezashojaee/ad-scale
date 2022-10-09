Exp_SO4 = [0.046652268	0.015550756	0.075161987	0.098488121	0.158099352	0.287688985	0.339524838	0.531317495	0.782721382	0.715334773	0.666090713	0.785313175	0.826781857	0.94600432]';
Ch_SO4 = [0.002591793	0.002591793	0.002591793	0.002591793	0.002591793	0.002591793	0.002591793	0.002591793	0.005183585	0.007775378	0.010367171	0.015550756	0.025917927	0.03887689	0.054427646	0.064794816	0.075161987	0.095896328	0.11663067	0.142548596	0.178833693	0.212526998	0.24362851	0.269546436	0.295464363	0.318790497	0.355075594	0.399136069	0.430237581	0.453563715	0.474298056	0.500215983	0.539092873	0.570194384	0.609071274	0.647948164	0.694600432	0.730885529	0.756803456	0.79049676	0.834557235	0.855291577	0.881209503	0.90712743	0.925269978	0.935637149]';
Exp_Mg = [0.186609071	0.025917927	0.316198704	0.756803456	0.90712743	0.90712743	0.881209503	0.964146868	1.010799136	1.003023758	1.028941685	0.982289417	1.026349892]';
Ch_Mg = [0	0	0	0.002591793	0.007775378	0.028509719	0.106263499	0.212526998	0.373218143	0.474298056	0.611663067	0.733477322	0.847516199	0.98488121	1.067818575	1.093736501	1.11187905	1.117062635	1.117062635]';
Exp_Cl = [0.137365011	0.041468683	0.090712743	0.326565875	0.694600432	0.800863931	0.850107991	0.852699784	0.904535637	0.966738661	0.977105832	0.992656587	0.982289417	0.997840173]';
Ch_Cl = [0	0	0	0	0	0	0	0	0	0.002591793	0.007775378	0.018142549	0.033693305	0.04924406	0.067386609	0.093304536	0.108855292	0.132181425	0.160691145	0.184017279	0.21511879	0.248812095	0.274730022	0.29287257	0.349892009	0.388768898	0.432829374	0.487257019	0.513174946	0.554643629	0.593520518	0.62462203	0.666090713	0.679049676	0.704967603	0.733477322	0.782721382	0.816414687	0.873434125	0.904535637	0.930453564	0.951187905	0.969330454	0.987473002	0.997840173	1.003023758	1.010799136	1.015982721	1.018574514	1.021166307	1.021166307	1.021166307]';
Exp_Ca = [0.034754721 0.327812387 0.605654858 0.708900018 0.784182913 0.856933655 0.894260754 0.924069044 0.93913566	0.952417265	0.959706327	0.966079551]';
Ch_Ca = [0	0	0	0	0	0	0	0	0	0	0.005183585	0.010367171	0.018142549	0.023326134	0.03887689	0.090712743	0.155507559	0.202159827	0.272138229	0.37062635	0.443196544	0.559827214	0.653131749	0.725701944	0.821598272	0.878617711	0.912311015	0.943412527	0.953779698	0.958963283	0.961555076	0.961555076	0.964146868]';
Exp_Na = [0.03887689	0.090712743	0.318790497	0.692008639	0.798272138	0.850107991	0.852699784	0.901943844	0.966738661	0.974514039	0.987473002	0.99524838]';
Ch_Na = [0	0	0	0	0	0	0	0	0.002591793	0.007775378	0.018142549	0.033693305	0.051835853	0.085529158	0.160691145	0.277321814	0.37062635	0.438012959	0.507991361	0.572786177	0.619438445	0.668682505	0.694600432	0.733477322	0.803455724	0.876025918	0.940820734	0.977105832	0.990064795	0.99524838	0.992656587]';
figure
plot(so4_exp,Exp_SO4,'bo',so4_model,Ch_SO4,so4_model,y_modelSO4)
legend('Experimental','Chandrasekhar et al. model','ad-scale')
grid
title('Sulfate Normalized Concentration')
xlabel('Injected PV')
ylabel('Normalized Concentration')
figure
plot(Ca_exp,Exp_Ca,'bo',Ca_model,Ch_Ca,Ca_model,y_modelCa)
legend('Experimental','Chandrasekhar et al. model','ad-scale')
grid
title('Calcium Normalized Concentration')
xlabel('Injected PV')
ylabel('Normalized Concentration')
figure
plot(Na_exp,Exp_Na,'bo',Na_model,Ch_Na,Na_model,y_modelNa)
legend('Experimental','Chandrasekhar et al. model','ad-scale')
grid
title('sodium Normalized Concentration')
xlabel('Injected PV')
ylabel('Normalized Concentration')
figure
plot(Cl_exp,Exp_Cl,'bo',Cl_model,Ch_Cl,Cl_model,y_modelCl)
legend('Experimental','Chandrasekhar et al. model','ad-scale')
grid
title('Chloride Normalized Concentration')
xlabel('Injected PV')
ylabel('Normalized Concentration')
figure
plot(Mg_exp,Exp_Mg,'bo',Mg_model,Ch_Mg,Mg_model,y_modelMg)
legend('Experimental','Chandrasekhar et al. model','ad-scale')
grid
title('Magnesium Normalized Concentration')
xlabel('Injected PV')
ylabel('Normalized Concentration')