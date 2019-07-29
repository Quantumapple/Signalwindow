#ifndef trkIso_PtFitFunction_h
#define trkIso_PtFitFunction_h

float pT_fit(int region, int combi, float dPhi)
{
    Float_t p[2] = {};
    Float_t x = fabs(dPhi);

    if( region == 1 )
    {
        if( combi == 1 ) {
            p[0] = 0.066493; p[1] = -0.608042;
        }
        if( combi == 2 || combi == 4 ) {
            p[0] = 0.0879698; p[1] = -0.475425;
        }
        if( combi == 3 ) {
            p[0] = 0.0876691; p[1] = -0.3938;
        }
    }
      
    if( region == 2 )
    {
        if( combi == 1 ) {
            p[0] = 0.066493; p[1] = -0.608042;
        }
        if( combi == 2 ) {
            p[0] = 0.0717142; p[1] = -0.769743;
        }
        if( combi == 3 || combi == 4 ) {
            p[0] = 0.0861495; p[1] = -0.366929;
        }
    }
    
    if( region == 3 )
    {
        if( combi == 1 ) {
            p[0] = 0.0717142; p[1] = -0.769743;
        }
        if( combi == 2 || combi == 4 ) {
            p[0] = 0.0750887; p[1] = -0.172686;
        }
        if( combi == 3 ) {
            p[0] = 0.0491693; p[1] = -0.705809;
        }
    }
    
    if( region == 4 )
    {
        if( combi == 1 ) {
            p[0] = 0.0491693; p[1] = -0.705809;
        }
        if( combi == 2 || combi == 4 ) {
            p[0] = 0.0520133; p[1] = -0.1448;
        }
        if( combi == 3 ) {
            p[0] = 0.050818; p[1] = -0.583872;
        }
    }
   
    if( region == 5 ) {
        if( combi == 1 || combi == 3 ) {
            p[0] = 0.0520133; p[1] = -0.1448;
        }
        if( combi == 2 || combi == 4 ) { 
            p[0] = 0.0534344; p[1] = -0.0481642;
        }
    }

    if( region == 6 )
    {
        if( combi == 1 ) {
            p[0] = 0.0534344; p[1] = -0.0481642;
        }
        if( combi == 2 || combi == 3 ) {
            p[0] = 0.0586581; p[1] = -0.0167858;
        }
        if( combi == 4 ) {
            p[0] = 0.0568861; p[1] = 0.398326;
        }
    }
    
    Float_t result = p[0]*(1./x)+p[1];
    return result;
}
#endif
