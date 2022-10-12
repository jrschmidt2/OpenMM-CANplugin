inline __device__ real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = make_real4(vec1.x-vec2.x,vec1.y-vec2.y,vec1.z-vec2.z,0);
    #ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(result);
    #endif
    result.w = result.x*result.x+result.y*result.y+result.z*result.z;
    return result;
}

inline __device__ real computeDot(real4 vec1, real4 vec2) {
    real result = vec1.x*vec2.x+vec1.y*vec2.y+vec1.z*vec2.z;
    return result;
}


inline __device__ real4 computeCross(real4 vec1, real4 vec2) {
    real3 result = make_real3(vec1.y*vec2.z-vec1.z*vec2.y, vec1.z*vec2.x-vec1.x*vec2.z, vec1.x*vec2.y-vec1.y*vec2.x);
    return make_real4(result.x,result.y,result.z,result.x*result.x + result.y*result.y + result.z*result.z);
}


//Angle Calculations
inline __device__ real computeAzim(real4 vec1, real4 vec2) {
    real angle;
    real dot = computeDot(vec1,vec2);
    vec1.w = computeDot(vec1,vec1);
    vec2.w = computeDot(vec2,vec2);
    real cosine = dot/sqrt(vec1.w*vec2.w);
    if (cosine > 0.99 || cosine < -0.99) {
        real4 cp = computeCross(vec1,vec2);
        angle = asin(sqrt(cp.w/(vec1.w*vec2.w)));
        if (cosine < 0.0) angle = M_PI-angle;
    }
    else angle = acos(cosine);
    return angle;
}

inline __device__ real computePolar(real4 vec1, real4 vec2, real4 vec3) {
    real angle;
    real4 cross1 = computeCross(vec1,vec2);
    real4 cross2 = computeCross(vec2,vec3);
    real dot = computeDot(cross1,cross2);
    real4 cross3 = computeCross(cross1,cross2);

    if (dot != 0.0) dot /= sqrt(cross1.w*cross2.w);
    if (dot > 1.0) dot = 1.0;
    else if (dot < -1.0) dot = -1.0;

    if (dot > 0.99 || dot < -0.99) {
        angle = asin(sqrt(cross3.w/(cross1.w*cross2.w)));
        if (dot < 0.0) angle = M_PI - angle;
    }
    else angle = acos(dot);
    return angle;
}

extern "C" __global__ void accessAxisParameter(const real4* __restrict__ posq, const int4* __restrict__ axes, real4* __restrict__ kvecs) {
    real4 vectorX,vectorY,vectorZ,vectemp;
    int axisType,pX,pY,pZ;
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom< NUM_ATOMS; atom += gridDim.x*blockDim.x) {
        __syncthreads();
        axisType = axes[atom].w;
        pX = axes[atom].x;
        pY = axes[atom].y;
        pZ = axes[atom].z;
        if (axisType == 5) {
            vectorZ = make_real4(0.0,0.0,1.0,1.0);
            vectorX = make_real4(1.0,0.0,0.0,1.0);
        }
        else {
            vectorZ = posq[pZ]-posq[atom];
            if (axisType == 4) {
                if (fabs(vectorZ.x) < 0.866f) vectorX = make_real4(1.0,0.0,0.0,1.0);
                else vectorX = make_real4(0.0,1.0,0.0,1.0);
            }
            else {
                vectorY = posq[pY]-posq[atom];
                if (axisType == 1) {
                    vectemp = computeCross(vectorZ,vectorY);
                    vectorZ += vectorY;
                    vectorZ *= 0.5;
                    vectorX = vectemp;
                }
                else {
                    vectorX = posq[pX]-posq[atom];
                    if (axisType == 2) {
                        vectorX += vectorY;
                        vectorX *= 0.5;
                    }
                    else if (axisType == 3) {
                        vectorZ += vectorX+vectorY;
                        vectorZ *= 1/3.0;
                        if (fabs(vectorZ.x) < 0.866f) vectorX = make_real4(1.0,0.0,0.0,1.0);
                        else vectorX = make_real4(0.0,1.0,0.0,1.0);
                    }
                }
            }
        }
        vectorX.w = computeDot(vectorX,vectorX);
        vectorZ.w = computeDot(vectorZ,vectorZ);
        vectorY = computeCross(vectorZ,vectorX);
        kvecs[atom*3] = vectorX;
        kvecs[atom*3+1] = vectorY;
        kvecs[atom*3+2] = vectorZ;
    }
}

extern "C" __global__ void computeCAN(unsigned long long* __restrict__ force, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, 
    const int4* __restrict__ axes, real4* __restrict__ kvecs,
    real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ 
    #ifdef USE_EXCLUSIONS
    , int* __restrict__ exclusions, int* __restrict exclusionStartIdx
    #endif
    PARAMETER_ARGUMENTS) {
    real dEdR=0.0; real dEdTheta1=0.0; real dEdTheta2=0.0; real dEdPhi1=0.0; real dEdPhi2=0.0;
    real4 f1=make_real4(0.0); real4 f2 = make_real4(0.0); real4 f3 = make_real4(0.0); real4 f4 = make_real4(0.0); real4 f0 = make_real4(0.0);
    int4 axs2, axs1;
    real4 kvec1x,kvec1y,kvec1z,kvec2x,kvec2y,kvec2z;
    mixed energy=0.0;


    int ii =  blockIdx.x*blockDim.x+threadIdx.x;
    if (ii < NUM_ATOMS) {
    axs1 = axes[ii]; 
    kvec1x = kvecs[ii*3];
    kvec1y = kvecs[ii*3+1];
    kvec1z = kvecs[ii*3+2];

    for (int jj = ii+1; jj < NUM_ATOMS; jj += 1) {
        __syncthreads();

        bool isExcluded=0;
        #ifdef USE_EXCLUSIONS
        int check = 0;
        int first = exclusionStartIdx[ii];
        int last = exclusionStartIdx[ii+1];
        for (int ex = last-1; ex >= first; ex--) {
            if (exclusions[ex] == jj) check+=1;
            if (check!=0) break;
        }
        if (check!=0) isExcluded = 1;
        #endif
        if (!isExcluded) {
            axs2 = axes[jj];
            kvec2x = kvecs[jj*3];
            kvec2y = kvecs[jj*3+1];
            kvec2z = kvecs[jj*3+2];
            //ixn calc here
            real4 rij = delta(posq[jj],posq[ii],periodicBoxSize,invPeriodicBoxSize,periodicBoxVecX,periodicBoxVecY,periodicBoxVecZ);
            #ifdef USE_CUTOFF
            if (rij.w < CUTOFF_SQUARED) {
            #endif
            real r = sqrt(rij.w);
            COMPUTE_FORCE
            f0 = dEdR*rij/r;
            //write results here
            //distance dependent forces
            atomicAdd(&force[ii],static_cast<unsigned long long>((long long) (f0.x*0x100000000)));
            atomicAdd(&force[ii+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f0.y*0x100000000)));
            atomicAdd(&force[ii+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f0.z*0x100000000)));
            atomicAdd(&force[jj],static_cast<unsigned long long>((long long) (-f0.x*0x100000000)));
            atomicAdd(&force[jj+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f0.y*0x100000000)));
            atomicAdd(&force[jj+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f0.z*0x100000000)));
            //angle expressions here
            //thetas
            if (dEdTheta1*dEdTheta1 > 0.0) {
                real4 thcross = computeCross(-kvec1z,-rij);
                real thcrossL = sqrt(thcross.w);
                if (thcrossL < 1.0e-06) thcrossL = 1.0e-06;
                real termA = dEdTheta1/(kvec1z.w*thcrossL);
                real termC = -dEdTheta1/(rij.w*thcrossL);
                f1 = termA*computeCross(-kvec1z,thcross);
                f3 = termC*computeCross(-rij,thcross);
                f2 = -f1-f3;
                //root atoms
                atomicAdd(&force[ii],static_cast<unsigned long long>((long long) (f2.x*0x100000000)));
                atomicAdd(&force[ii+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f2.y*0x100000000)));
                atomicAdd(&force[ii+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f2.z*0x100000000)));
                atomicAdd(&force[jj],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                atomicAdd(&force[jj+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                atomicAdd(&force[jj+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                __threadfence_block();
                //axis atoms
                if(axs1.w == 4 || axs1.w == 0 || axs1.w == 2) {
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs1.w == 1) {
                    f1 = 0.5*f1;
                    atomicAdd(&force[axs1.y],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs1.w == 3) {
                    f1 = f1/3.0;
                    atomicAdd(&force[axs1.x],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs1.y],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
                }
            }
            if (dEdTheta2*dEdTheta2 > 0.0) {
                real4 thcross = computeCross(-kvec2z,rij);
                real thcrossL = sqrt(thcross.w);
                if (thcrossL < 1.0e-06) thcrossL = 1.0e-06;
                real termA = dEdTheta2/(kvec2z.w*thcrossL);
                real termC = -dEdTheta2/(rij.w*thcrossL);
                f1 = termA*computeCross(-kvec2z,thcross);
                f3 = termC*computeCross(rij,thcross);
                f2 = -f1-f3;
                //root atoms
                atomicAdd(&force[jj],static_cast<unsigned long long>((long long) (f2.x*0x100000000)));
                atomicAdd(&force[jj+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f2.y*0x100000000)));
                atomicAdd(&force[jj+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f2.z*0x100000000)));
                atomicAdd(&force[ii],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                atomicAdd(&force[ii+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                atomicAdd(&force[ii+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                __threadfence_block();
                //axis atoms
                if(axs2.w == 4 || axs2.w == 0 || axs2.w == 2) {
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs2.w == 1) {
                    f1 = 0.5*f1;
                    atomicAdd(&force[axs2.y],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs2.w == 3) {
                    f1 = f1/3.0;
                    atomicAdd(&force[axs2.x],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs2.y],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                __threadfence_block();
                }
            }
            if (dEdPhi1*dEdPhi1 > 0.0) {
                real4 cross2 = computeCross(kvec1z,kvec1z-rij);
                if (cross2.w < 1.0e-06) cross2.w = 1.0e-06;
                real normBC = sqrt(kvec1z.w);
                f1 = dEdPhi1*normBC*kvec1y/kvec1y.w;
                f4 = dEdPhi1*normBC*cross2/cross2.w;
                real dot1 = computeDot(kvec1x,kvec1z)/kvec1z.w;
                real dot2 = computeDot(kvec1z-rij,kvec1z)/kvec1z.w;
                real4 ss = dot1*f1-dot2*f4;
                f2 = f1-ss;
                f3 = f4+ss; 
                atomicAdd(&force[ii],static_cast<unsigned long long>((long long) (-f2.x*0x100000000)));
                atomicAdd(&force[ii+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f2.y*0x100000000)));
                atomicAdd(&force[ii+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f2.z*0x100000000)));
                atomicAdd(&force[jj],static_cast<unsigned long long>((long long) (f4.x*0x100000000)));
                atomicAdd(&force[jj+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f4.y*0x100000000)));
                atomicAdd(&force[jj+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f4.z*0x100000000)));
                __threadfence_block();
                if (axs1.w == 0) {
                    atomicAdd(&force[axs1.x],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs1.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs1.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (-f3.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs1.w == 1) {
                    atomicAdd(&force[axs1.y],static_cast<unsigned long long>((long long) (-0.5*f3.x*0x100000000)));
                    atomicAdd(&force[axs1.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.y*0x100000000)));
                    atomicAdd(&force[axs1.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (-0.5*f3.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.z*0x100000000)));
                    real4 f5 = computeCross(posq[axs1.z]-posq[ii],f1);
                    real4 f6 = computeCross(f1,posq[axs1.y]-posq[ii])+0.5*f1;
                    atomicAdd(&force[axs1.y],static_cast<unsigned long long>((long long) (f5.x*0x100000000)));
                    atomicAdd(&force[axs1.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f5.y*0x100000000)));
                    atomicAdd(&force[axs1.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f5.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (f6.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f6.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f6.z*0x100000000)));
                    atomicAdd(&force[ii],static_cast<unsigned long long>((long long) ((f1.x-f5.x-f6.x)*0x100000000)));
                    atomicAdd(&force[ii+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) ((f1.y-f5.y-f6.y)*0x100000000)));
                    atomicAdd(&force[ii+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) ((f1.z-f5.z-f6.z)*0x100000000)));
                __threadfence_block();
                }
                else if (axs1.w == 2) {
                    atomicAdd(&force[axs1.x],static_cast<unsigned long long>((long long) (0.5*f1.x*0x100000000)));
                    atomicAdd(&force[axs1.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.y*0x100000000)));
                    atomicAdd(&force[axs1.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.z*0x100000000)));
                    atomicAdd(&force[axs1.y],static_cast<unsigned long long>((long long) (0.5*f1.x*0x100000000)));
                    atomicAdd(&force[axs1.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.y*0x100000000)));
                    atomicAdd(&force[axs1.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (-f3.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs1.w == 3) {
                    f3 = -f3/3.0; 
                    atomicAdd(&force[axs1.x],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                    atomicAdd(&force[axs1.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                    atomicAdd(&force[axs1.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                    atomicAdd(&force[axs1.y],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                    atomicAdd(&force[axs1.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                    atomicAdd(&force[axs1.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs1.w == 4) {
                    atomicAdd(&force[axs1.z],static_cast<unsigned long long>((long long) (-f3.x*0x100000000)));
                    atomicAdd(&force[axs1.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.y*0x100000000)));
                    atomicAdd(&force[axs1.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.z*0x100000000)));
                __threadfence_block();
                }
            }
            if (dEdPhi2*dEdPhi2 > 0.0) {
                real4 cross2 = computeCross(kvec2z,kvec2z+rij);
                if (cross2.w < 1.0e-06) cross2.w = 1.0e-06;
                real normBC = sqrt(kvec2z.w);
                f1 = dEdPhi2*normBC*kvec2y/kvec2y.w;
                f4 = dEdPhi2*normBC*cross2/cross2.w;
                real dot1 = computeDot(kvec2x,kvec2z)/kvec2z.w;
                real dot2 = computeDot(kvec2z+rij,kvec2z)/kvec2z.w;
                real4 ss = dot1*f1-dot2*f4;
                f2 = f1-ss;
                f3 = f4+ss;
                atomicAdd(&force[jj],static_cast<unsigned long long>((long long) (-f2.x*0x100000000)));
                atomicAdd(&force[jj+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f2.y*0x100000000)));
                atomicAdd(&force[jj+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f2.z*0x100000000)));
                atomicAdd(&force[ii],static_cast<unsigned long long>((long long) (f4.x*0x100000000)));
                atomicAdd(&force[ii+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f4.y*0x100000000)));
                atomicAdd(&force[ii+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f4.z*0x100000000)));
                __threadfence_block();
                if (axs2.w == 0) {
                    atomicAdd(&force[axs2.x],static_cast<unsigned long long>((long long) (f1.x*0x100000000)));
                    atomicAdd(&force[axs2.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.y*0x100000000)));
                    atomicAdd(&force[axs2.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f1.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (-f3.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs2.w == 1) {
                    atomicAdd(&force[axs2.y],static_cast<unsigned long long>((long long) (-0.5*f3.x*0x100000000)));
                    atomicAdd(&force[axs2.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.y*0x100000000)));
                    atomicAdd(&force[axs2.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (-0.5*f3.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-0.5*f3.z*0x100000000)));
                    real4 f5 = computeCross(posq[axs2.z]-posq[jj],f1);
                    real4 f6 = computeCross(f1,posq[axs2.y]-posq[jj])+0.5*f1;
                    atomicAdd(&force[axs2.y],static_cast<unsigned long long>((long long) (f5.x*0x100000000)));
                    atomicAdd(&force[axs2.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f5.y*0x100000000)));
                    atomicAdd(&force[axs2.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f5.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (f6.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f6.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f6.z*0x100000000)));
                    atomicAdd(&force[jj],static_cast<unsigned long long>((long long) ((f1.x-f5.x-f6.x)*0x100000000)));
                    atomicAdd(&force[jj+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) ((f1.y-f5.y-f6.y)*0x100000000)));
                    atomicAdd(&force[jj+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) ((f1.z-f5.z-f6.z)*0x100000000)));
                __threadfence_block();
                }
                else if (axs2.w == 2) {
                    atomicAdd(&force[axs2.x],static_cast<unsigned long long>((long long) (0.5*f1.x*0x100000000)));
                    atomicAdd(&force[axs2.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.y*0x100000000)));
                    atomicAdd(&force[axs2.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.z*0x100000000)));
                    atomicAdd(&force[axs2.y],static_cast<unsigned long long>((long long) (0.5*f1.x*0x100000000)));
                    atomicAdd(&force[axs2.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.y*0x100000000)));
                    atomicAdd(&force[axs2.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (0.5*f1.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (-f3.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs2.w == 3) {
                    f3 = -f3/3.0;
                    atomicAdd(&force[axs2.x],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                    atomicAdd(&force[axs2.x+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                    atomicAdd(&force[axs2.x+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                    atomicAdd(&force[axs2.y],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                    atomicAdd(&force[axs2.y+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                    atomicAdd(&force[axs2.y+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (f3.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (f3.z*0x100000000)));
                __threadfence_block();
                }
                else if (axs2.w == 4) {
                    atomicAdd(&force[axs2.z],static_cast<unsigned long long>((long long) (-f3.x*0x100000000)));
                    atomicAdd(&force[axs2.z+PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.y*0x100000000)));
                    atomicAdd(&force[axs2.z+2*PADDED_NUM_ATOMS],static_cast<unsigned long long>((long long) (-f3.z*0x100000000)));
                __threadfence_block();
                }
            }
            #ifdef USE_CUTOFF
            }
            #endif
        }
        }
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
