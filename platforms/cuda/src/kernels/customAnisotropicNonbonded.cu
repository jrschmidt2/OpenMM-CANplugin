//CustomAnisotropicNonbonded

//Define useful functions.. first 4 in CustomHbond
//Vector simplification
inline __device__ real3 trim(real4 v) {
	return make_real3(v.x,v.y,v.z);
}
inline __device__ real3 trim(real3 v) {
	return v;
}
//Periodic and nonperiodic difference
inline __device__ real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
	real4 result = make_real4(vec1.x-vec2.x,vec1.y-vec2.y,vec1.z-vec2.z,0.0f);
	#ifdef USE_PERIODIC
		APPLY_PERIODIC_TO_DELTA(result)
	#endif
	result.w = result.x*result.x+result.y*result.y+result.z*result.z;
	return result
}
inline __device__ real4 computeCross(real4 vec1, real4 vec2) {
	real3 result = cross(vec1,vec2);
	return make_real4(result.x,result.y,result.z,result.x*result.x + result.y*result.y + result.z*result.z)
}

//Angle Calculations
inline __device__ real computePhi(real4 vec1, real4 vec2) {
	real angle;
	real dot = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
	double cosine = dot/RSQRT(vec1.w*vec2.w);
	if (cosine >= 1) angle = 0;
	else if (cosine <= -1) angle = M_PI;
	else angle = acos(cosine);
	return angle
}
inline __device__ real computeTheta(real4 vec1, real4 vec2, real4 vec3) {
	real angle;
	real4 cross1 = computeCross(vec1,vec2);
	real4 cross2 = computeCross(vec2,vec3);
	real dot = cross1.x*cross2.x+cross1.y*cross2.y+cross1.z*cross2.z;
	real4 cross = computeCross(cross1,cross2);

	if (dot != 0.0) dot /= RSQRT(cross1.w*cross2.w);
	if (dot > 1.0) dot = 1.0;
	else if (dot < -1.0) dot = -1.0;

	if (dot > 0.99 || dot < -0.99) {
		angle = ASIN(RSQRT(cross.w/(cross1.w*cross2.w)));
		if (dot < 0.0) angle = M_PI - angle;
	}
	else angle = ACOS(dot);
	return angle;
}

inline __device__ void accessAxisParameter(real4* pos, int4 kparticles, real3* kvecs) {
	real3 vectorX,vectorY,vectorZ,vectemp;
	int axisType,pX,pY,pZ;
	for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom< NUM_ATOMS; atom += gridDim.x*blockDim.x) {
		axisType = kparticles.w;
		pX = kparticles.x;
		pY = kparticles.y;
		pZ = kparticles.z;
		if (axisType == CustomAnisotropicNonbondedForce::NoAxisType) {
			vectorZ = make_real3(0,0,1);
			vectorX = make_real3(1,0,0);
		}
		else {
			vectorZ = pos[pZ] - pos[atom];
			if (axisType == CustomAnisotropicNonbondedForce::ZOnly) {
				if (ABS(vectorZ.x) < 0.866) vectorX = make_real3(1.0,0.0,0.0);
				else vectorX = make_real3(0.0,1.0,0.0);
			}
			else {
				vectorX = pos[pX] - pos[atom];
				if (axisType == CustomAnisotropicNonbondedForce::Bisector) {
					vectemp = cross(vectorZ,vectorX);
					vectorZ += vectorX;
					vector *= 0.5;
					vectorX = vectemp;
				}
				else {
					vectorY = pos[pY] - pos[atom];
					if (axisType == CustomAnisotropicNonbondedForce::ZBisect) {
						vectorX += vectorY;
						vectorX *= 0.5;
					}
					else if (axisType == CustomAnisotropicNonbondedForce::ThreeFold) {
						vectorZ += vectorX + vectorY;
						vectorZ *= 1/3.0;
						if (ABS(vectorZ.x) < 0.866) vectorX = make_real3(1.0,0.0,0.0);
						else vectorX = make_real3(0.0,1.0,0.0);
					}
				}
			}
		}
		vectorY = cross(vectorZ,vectorX);
		kvecs[atom] = vectorX;
		kvecs[atom*3+1] = vectorY;
		kvecs[atom*3+2] = vectorZ;
	}
}


/*
#ifdef USE_CUTOFF
if (!isExcluded && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded) {
#endif

	real tempForce = 0;
	reala switchValue = 1, switchDeriv = 0;
#if USE_SWITCH
	if (r > SWITCH_CUTOFF) {
		real x = r - SWITCH_CUTOFF;
		switchValue = 1+x*x*x*(SWITCH_C3+x*(SWITCH_C4+x*SWITCH_C5));
		switchDeriv = x*x*(3*SWITCH_C3+x*(4*SWITCH_C4+x*5*SWITCH_C5));
	}
#endif
	COMPTE_FORCE
#if USE_SWITCH
	tempForce = tempForce*switchValue - tempEnergy*switchDeriv;
	tempEnergy *= switchValue;
#endif
	dEdR += tempForce*invR;
}
*/
