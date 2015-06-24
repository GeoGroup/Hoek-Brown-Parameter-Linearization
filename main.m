m=2.0121;s=3.87E-03;a=0.506;sgm0=10;sgmc=80;sgmi=0;
 [c,f]=elastic_plastic(m,s,a,sgmc,sgm0,sgmi);
 m=2.0121;s=3.87E-03;a=0.506;sgmc=80;
 mr=1.1095;sr=1.27E-03;ar=0.511;sgmcr=80;
 sgm0=10;sgmi=0;
 [cp,fp,cr,fr]=elastic_brittle_plastic(m,s,a,sgmc,mr,sr,ar,sgmcr,sgm0,sgmi);