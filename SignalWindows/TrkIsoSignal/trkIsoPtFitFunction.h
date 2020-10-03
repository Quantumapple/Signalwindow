#ifndef trkIso_PtFitFunction_h
#define trkIso_PtFitFunction_h

float pT_fit(int region, int combi, float dPhi)
{
	Float_t x = fabs(dPhi);

	// Function a/x + b
	// parameter 'a' 2-d array
	Float_t p1[6][4] = {
		{0.066493,  0.0879698, 0.0876691, 0.0879698}, // 1st region
		{0.066493,  0.0717142, 0.0861495, 0.0861495}, // 2nd region
		{0.0717142, 0.0750887, 0.0491693, 0.0750887}, // 3rd region
		{0.0491693, 0.0520133, 0.050818,  0.0520133}, // 4th region
		{0.0520133, 0.0534344, 0.0520133, 0.0534344}, // 5th region
		{0.0534344, 0.0586581, 0.0586581, 0.0568861}  // 6th region
	};
	// parameter 'b' 2-d array
	Float_t p2[6][4] = {
		{-0.608042,  -0.475425,  -0.3938,   -0.475425}, // 1st region
		{-0.608042,  -0.769743,  -0.366929, -0.366929}, // 2nd region
		{-0.769743,  -0.172686,  -0.705809, -0.172686}, // 3rd region
		{-0.705809,  -0.1448,    -0.583872, -0.1448  }, // 4th region
		{-0.1448,    -0.0481642, -0.1448,   -0.0481642}, // 5th region
		{-0.0481642, -0.0167858, -0.0167858, 0.398326 }  // 6th region
	};
	
	Float_t result = p1[region-1][combi-1]/x + p2[region-1][combi-1];
	return result;
}
#endif
