#ifndef RegionOfInterest_h
#define RegionOfInterest_h

double ROI_func(int region, int nth, double eget){
  double p[4];

if( region == 1 &&  nth == 0 ){
p[0] = 0.000201154;
p[1] = -1.20559;
p[2] = -1.09273;
p[3] = -0.75382;
}


if( region == 1 &&  nth == 1 ){
p[0] = 0.00137851;
p[1] = -0.784711;
p[2] = -0.98718;
p[3] = -7.06446;
}


if( region == 1 &&  nth == 2 ){
p[0] = 0.000368831;
p[1] = -1.31085;
p[2] = -1.11991;
p[3] = -0.57669;
}


if( region == 1 &&  nth == 3 ){
p[0] = 0.00167382;
p[1] = -0.714873;
p[2] = -0.977424;
p[3] = -37.9691;
}

if( region == 2 &&  nth == 0 ){
p[0] = 0.000764339;
p[1] = -0.78772;
p[2] = -0.992126;
p[3] = -15.8588;
}


if( region == 2 &&  nth == 1 ){
p[0] = 0.000741569;
p[1] = -0.765971;
p[2] = -0.994122;
p[3] = -15.5259;
}


if( region == 2 &&  nth == 2 ){
p[0] = 0.00108463;
p[1] = -0.720655;
p[2] = -0.980618;
p[3] = -118.347;
}


if( region == 2 &&  nth == 3 ){
p[0] = -0.0006604;
p[1] = -1.21195;
p[2] = -1.13784;
p[3] = -0.653238;
}


if( region == 3 &&  nth == 0 ){
p[0] = -4.78824e-05;
p[1] = -0.657584;
p[2] = -0.964997;
p[3] = -0.685739;
}


if( region == 3 &&  nth == 1 ){
p[0] = 0.000276935;
p[1] = -0.920694;
p[2] = -0.992173;
p[3] = -0.305891;
}


if( region == 3 &&  nth == 2 ){
p[0] = -0.000111358;
p[1] = -0.630247;
p[2] = -0.97236;
p[3] = -0.688998;
}


if( region == 3 &&  nth == 3 ){
p[0] = -0.00020583;
p[1] = -0.5859;
p[2] = -0.964334;
p[3] = -0.694892;
}

if( region == 4 &&  nth == 0 ){
p[0] = -0.000546603;
p[1] = -0.518232;
p[2] = -1.00133;
p[3] = -0.680051;
}


if( region == 4 &&  nth == 1 ){
p[0] = -0.000773789;
p[1] = -0.537211;
p[2] = -1.02968;
p[3] = -0.67315;
}


if( region == 4 &&  nth == 2 ){
p[0] = -0.00333379;
p[1] = -0.283478;
p[2] = 0.0958272;
p[3] = 0.368868;
}


if( region == 4 &&  nth == 3 ){
p[0] = -0.000600583;
p[1] = -0.591073;
p[2] = -1.04961;
p[3] = -0.492056;
}


if( region == 5 &&  nth == 0 ){
p[0] = -4.83035e-05;
p[1] = -0.310813;
p[2] = -0.981197;
p[3] = -0.69223;
}


if( region == 5 &&  nth == 1 ){
p[0] = 0.000439714;
p[1] = -0.191644;
p[2] = -0.867853;
p[3] = -92.0521;
}


if( region == 5 &&  nth == 2 ){
p[0] = -0.00140326;
p[1] = -0.264927;
p[2] = -0.181628;
p[3] = 0.315235;
}


if( region == 5 &&  nth == 3 ){
p[0] = -0.00195718;
p[1] = -0.189465;
p[2] = 0.0589701;
p[3] = 0.370528;
}



if( region == 6 &&  nth == 0 ){
p[0] = -0.000216748;
p[1] = -0.440745;
p[2] = -1.07862;
p[3] = -0.130544;
}


if( region == 6 &&  nth == 1 ){
p[0] = -0.00190377;
p[1] = -0.122485;
p[2] = 0.183526;
p[3] = 0.408082;
}


if( region == 6 &&  nth == 2 ){
p[0] = -0.000904923;
p[1] = -0.335985;
p[2] = -1.2196;
p[3] = -0.600006;
}


if( region == 6 &&  nth == 3 ){
p[0] = -0.00152215;
p[1] = -0.185752;
p[2] = -0.105493;
p[3] = 0.372795;
}



  return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3]));

}

#endif