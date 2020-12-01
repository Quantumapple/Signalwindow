#ifndef trkIso_PtFitFunction_h
#define trkIso_PtFitFunction_h

float pT_fit(int region, int combi, float dPhi)
{
	Float_t x = fabs(dPhi);

	// Function a/x + b
	// parameter 'a' 2-d array
	Float_t p1[6][4] = {
		{0.05,  0.073, 0.073, 0.05}, // 1st region
		{0.05,  0.051, 0.070, 0.047}, // 2nd region
		{0.051, 0.059, 0.036, 0.036}, // 3rd region
		{0.036, 0.041, 0.041, 0.021}, // 4th region
		{0.021, 0.031, 0.03,  0.02}, // 5th region
		{0.02,  0.03,  0.03,  0.02}  // 6th region
	};
	
	Float_t result = p1[region-1][combi-1]/x;
	return result;
}
#endif
