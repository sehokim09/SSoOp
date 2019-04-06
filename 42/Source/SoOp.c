#include "42.h"
#define EXTERN
#include "42GlutGui.h"
#undef EXTERN

double Limit_max(double x, double max)
{
      return(x > max ? max : x);
}
void ECEFToWGS84_LAT(double p[3], double *glat)
{
      double a = 6378137.0;
      double f = 1.0/298.257223563;
      double b = a*(1.0-f);
      double e2 = f*(2.0-f);
      double ep2 = f*(2.0-f)/(1.0-f)/(1.0-f);
      double r,E2,F,G,C,S,P,Q,r0,V,Z0;

      double OneMinusE2,Z1,SpolyG,Qpoly;

      OneMinusE2 = 1.0-e2;

      r = sqrt(p[0]*p[0]+p[1]*p[1]);

      E2 = a*a-b*b;

      Z1 = b*p[2];

      F = 54.0*Z1*Z1;

      G = r*r+OneMinusE2*p[2]*p[2]-e2*E2;

      Z1 = e2*r/G;
      C = Z1*Z1*F/G;

      S = pow(1.0+C+sqrt(C*C+2.0*C),1.0/3.0);

      SpolyG = (S+1.0/S+1.0)*G;
      P = F/(3.0*SpolyG*SpolyG);

      Q = sqrt(1.0+2.0*e2*e2*P);

      Qpoly = 1.0+Q;
      r0 = -P*e2*r/Qpoly+sqrt(0.5*a*a*Qpoly/Q-P*OneMinusE2*p[2]*p[2]/(Q*Qpoly)-0.5*P*r*r);

      Z1 = r-e2*r0;
      Z1 *= Z1;

      V = sqrt(Z1+OneMinusE2*p[2]*p[2]);

      Z1 = b*b/a/V;
      Z0 = Z1*p[2];

      *glat = atan((p[2]+ep2*Z0)/r);
}

long FindSP_MPL(struct SCType *SCT, struct SCType *SCR, struct SpecularType *Sp)
{
	double PosWT[3], PosWR[3], PosWT_UNIT[3], PosWR_UNIT[3];
	double Lat;
	double estCurr[3], corrEst[3], corrEst_UNIT[3], estNext[3], estGradi[3];
	double estCurrNor[3], temp[3], corr[3];
	double ST[3], RT[3];
	double MAGR, MAGT, MAGC;
	double K = 100000.0;
	double sigma = 0.1;
	double M[3][3] = {{0.0}};
	// WGS84
    double a = 6378137.0;
    double f = 1.0/298.257223563;
    double b = a*(1.0-f);
    double r;
    long cnt = 0;
    long i;

    M[0][0] = 2/a/a;
    M[1][1] = 2/a/a;
    M[2][2] = 2/b/b;

	for(i=0;i<3;i++) {
		corr[i] = 100.0;
	}
	/* ECI to ECEF */
	MxV(World[EARTH].CWN,SCT->PosN,PosWT);
	MxV(World[EARTH].CWN,SCR->PosN,PosWR);

	MAGT = CopyUnitV(PosWT,PosWT_UNIT);
	MAGR = CopyUnitV(PosWR,PosWR_UNIT);
    /* Initial Estimate */
	/* J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates \
    to geodetic coordinates," IEEE Transactions on Aerospace and \
    Electronic Systems, vol. 30, pp. 957-961, 1994.*/
	//ECEFToWGS84_LAT(PosWR,&Lat);
    //Lat = asin(PosWR[2]/MAG);
    /* radius of WGS84 ellipsoid at Latitude*/
    //r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
	//SxV(r, PosWR_UNIT, estCurr);

	/* Unconstrained Weighted Initial Estimate */
	for(i=0;i<3;i++){
		estCurr[i] = PosWR_UNIT[i]/MAGR + PosWT_UNIT[i]/MAGT;
	}

	while(cnt < 2000 && MAGV(corr) > sigma) {
	    /* Step 1 */
		for(i=0;i<3;i++) {
			ST[i] = estCurr[i]-PosWT[i];
			RT[i] = estCurr[i]-PosWR[i];
		}
		UNITV(ST);
		UNITV(RT);
		for(i=0;i<3;i++) {
			estGradi[i] = ST[i] + RT[i];
		}

	    /* Step 2 */
		// Tangential correction optimization(long RxNum, long TxNum, struct SpecularType *Sp)
		MxV(M,estCurr,estCurrNor);
		VxV(estCurrNor,estGradi,temp);
		VxV(temp, estCurrNor, estGradi);

		SxV(K, estGradi, estGradi);
		for(i=0;i<3;i++) {
			corrEst[i] = estCurr[i] - estGradi[i];
		}

	    /* Step 3 */
		//ECEFToWGS84_LAT(corrEst,&Lat);
	    MAGC = CopyUnitV(corrEst,corrEst_UNIT);
	    Lat = asin(corrEst[2]/MAGC);
	    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
		SxV(r, corrEst_UNIT, estNext);

	    /* Step 4 */
		for(i=0;i<3;i++) {
			corr[i] = fabs(estNext[i]-estCurr[i]);
			estCurr[i] = estNext[i];
		}
		cnt++;
	}
	/* ECEF to WGS84 */
	ECEFToWGS84(estCurr, &Sp->lat, &Sp->lon, &Sp->alt);

	/* ECEF to ECI */
	MTxV(World[EARTH].CWN,estCurr,Sp->PosN);

	return cnt;
}

void FindSP_UD(double estInit[3], struct SpecularType *Sp)
{
	double Lat;
	double estCurr[3], corrEst[3], corrEst_UNIT[3], estNext[3], estGradi[3];
	double estCurrNor[3], temp[3], corr[3];
	double SurfaceNorVec[3], BistaticVec[3];
	double MAGC;
	double K = 100000.0;
	double sigma = 0.1;
	double M[3][3] = {{0.0}};
	// WGS84
    double a = 6378137.0;
    double f = 1.0/298.257223563;
    double b = a*(1.0-f);
    double r;

    M[0][0] = 2/a/a;
    M[1][1] = 2/a/a;
    M[2][2] = 2/b/b;

	corr[0] = 100.0;
	corr[1] = 100.0;
	corr[2] = 100.0;

	/* Unconstrained Weighted Initial Estimate */
	estCurr[0] = estInit[0];
	estCurr[1] = estInit[1];
	estCurr[2] = estInit[2];

	while(MAGV(corr) > sigma) {
	    /* Step 1: Calculate gradient of path length function */
		MxV(M,estCurr,estCurrNor);
		CopyUnitV(estCurrNor, SurfaceNorVec);
		CopyUnitV(estCurr, BistaticVec);

		estGradi[0] = BistaticVec[0] - SurfaceNorVec[0];
		estGradi[1] = BistaticVec[1] - SurfaceNorVec[1];
		estGradi[2] = BistaticVec[2] - SurfaceNorVec[2];

	    /* Step 2 */
		// Tangential correction optimization
		VxV(estCurrNor,estGradi,temp);
		VxV(temp, estCurrNor, estGradi);

		SxV(K, estGradi, estGradi);

		corrEst[0] = estCurr[0] - estGradi[0];
		corrEst[1] = estCurr[1] - estGradi[1];
		corrEst[2] = estCurr[2] - estGradi[2];

	    /* Step 3 */
	    MAGC = CopyUnitV(corrEst,corrEst_UNIT);
	    Lat = asin(corrEst[2]/MAGC);
	    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
		SxV(r, corrEst_UNIT, estNext);

	    /* Step 4 */
		corr[0] = fabs(estNext[0]-estCurr[0]);
		corr[1] = fabs(estNext[1]-estCurr[1]);
		corr[2] = fabs(estNext[2]-estCurr[2]);
		estCurr[0] = estNext[0];
		estCurr[1] = estNext[1];
		estCurr[2] = estNext[2];
	}

	Sp->PosW[0] = estCurr[0];
	Sp->PosW[1] = estCurr[1];
	Sp->PosW[2] = estCurr[2];

	/* ECEF to WGS84 */
	ECEFToWGS84(estCurr, &Sp->lat, &Sp->lon, &Sp->alt);

	/* ECEF to ECI */
	MTxV(World[EARTH].CWN,estCurr,Sp->PosN);
}

/**********************************************************************/
/* Calculate swath area of FOV. Assume spherical Earth, symmetric FOV */
void FindSwathArea(long RxNum)
{
	double alpha;
	double alpha1;
	double alpha2;
	double CBL[3][3];
	double vl[3] = {0.0, 0.0, 1.0};
	double vb[3];
	double theta, theta1, theta2, thetaMax;
	double halfFOV;
	double ratio;

	ratio = MAGV(SC[RxNum].PosN)/World[EARTH].rad;
	thetaMax = asin(1/ratio);
	halfFOV = FOV[RxNum].Width/2;

	/* DCM from Local to Body */
	MxMT(SC[RxNum].B[0].CN,SC[RxNum].CLN,CBL);
	/* +Z vector of Local frame is expressed in Body frame */
	MxV(CBL,vl,vb);
	/* Angle between +Z vector in Local and Body */
	theta = acos(VoV(vl,vb));

	/* Take into account FOV */
	if(theta > halfFOV){
		theta1 = theta - halfFOV;
		theta2 = Limit_max(theta + halfFOV,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 - alpha1;
	}
	else if(theta < halfFOV){
		theta1 = Limit_max(halfFOV - theta,thetaMax);
		theta2 = Limit_max(halfFOV + theta,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 + alpha1;
	}
	else{
		theta1 = theta;
		theta2 = Limit_max(halfFOV + theta,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 + alpha1;
	}

	Rx[RxNum].SwathWidth = alpha/TwoPi*World[EARTH].rad;
	Rx[RxNum].SwathArea = TwoPi*World[EARTH].rad*World[EARTH].rad*(1-cos(alpha/2));
}

long *FindVisibleTx(long RxNum, long SoOpNum, long ITxStart, long ITxEnd, long *TxNum, long *TxCnt)
{
	static long buffer[25];
	double ToS;
	double Rhat[3], RelPosN[3];
	long ITx, i;
	long cnt = 0;
	double elevation;
	double MinElevation = -16*D2R;


	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(SC[RxNum].PosN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    for(ITx=ITxStart;ITx<ITxEnd+1;ITx++) {
    //   if (SC[ITx].Exists) {
          for(i=0;i<3;i++)
             RelPosN[i] = SC[ITx].PosN[i] - SC[RxNum].PosN[i];
          UNITV(RelPosN);
          ToS = VoV(RelPosN,Rhat);
          elevation = HalfPi - acos(ToS);
          if(elevation > MinElevation){
        	  buffer[cnt++] = ITx;
          }
    //   }
    }
    *TxCnt = cnt;

	TxNum = (long *) calloc(cnt,sizeof(long));

	for(ITx=0; ITx<cnt; ITx++) {
		TxNum[ITx] = buffer[ITx];
	}
	return TxNum;
}

long FindClosestTx(long RxNum, long ITxStart, long ITxEnd)
{
	//double MaxToS = -1.0; /* Bogus */
	double ToS;
	double Rhat[3], RelPosN[3];
	long ITx, i;
	long TxNum = 0;
	double elevation;
	double MinElevation = -20*D2R;
	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(SC[RxNum].PosN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    for(ITx=ITxStart;ITx<ITxEnd+1;ITx++) {
       if (SC[ITx].Exists) {
          for(i=0;i<3;i++)
             RelPosN[i] = SC[ITx].PosN[i] - SC[RxNum].PosN[i];
          UNITV(RelPosN);
          ToS = VoV(RelPosN,Rhat);
          elevation = HalfPi - acos(ToS);
          if(elevation > MinElevation){
        	  MinElevation = elevation;
          //if (ToS > MaxToS) {
           //  MaxToS = ToS;
             TxNum = ITx;
          }
       }
    }
    //printf("ele:%f\n", MinElevation*R2D);
    return TxNum;
}

long *FindBistaticGeo(long RxNum, long SoOpNum, long TxStart, long TxEnd)
{
	long ISp;
	struct SpecularType *Sp;
	long *TxNum = NULL;
	long TxCnt;

	TxNum = FindVisibleTx(RxNum, SoOpNum, TxStart, TxEnd, TxNum, &TxCnt);

	for(ISp=0; ISp<TxCnt; ISp++){
		Sp = &Rx[RxNum].SoOp[SoOpNum].Sp[ISp];
		if(&TxNum+ISp==0) {
			Sp->Exists = 0;
			Sp->PosN[0] = 0.0;
			Sp->PosN[1] = 0.0;
			Sp->PosN[2] = 0.0;
			Sp->lat = 0.0;
			Sp->lon = 0.0;
			Sp->alt = 0.0;
		}
		else {
			Sp->Exists = 1;
			Sp->TxNum = TxNum[ISp];
	        //FindSP_UD(RxNum, TxNum[ISp], Sp);
		}
	}
    return TxNum;
}

void FindPathLength(long RxNum, long SoOpNum, long SpNum, long TxNum)
{
	long i;
	double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3];
	double MagVecSp2Rx, MagVecSp2Tx;

	if(Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].Exists) {
		for(i=0;i<3;i++) {
			VecTx2Rx[i] = SC[TxNum].PosN[i] - SC[RxNum].PosN[i];
			VecSp2Rx[i] = SC[RxNum].PosN[i] - Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PosN[i];
			VecSp2Tx[i] = SC[TxNum].PosN[i] - Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PosN[i];
		}
		MagVecSp2Rx = UNITV(VecSp2Rx);
		MagVecSp2Tx = UNITV(VecSp2Tx);
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenD = MAGV(VecTx2Rx);
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenR = MagVecSp2Rx + MagVecSp2Tx;
	}
	else {
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenD = 0.0;
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenR = 0.0;
	}
}

void CalcRxPwr(long RxNum, long SoOpNum, long SpNum, long TxNum, double Reflectivity)
{
	struct SoOpType *SoOp;
	struct SpecularType *Sp;

	SoOp = &Rx[RxNum].SoOp[SoOpNum];
	Sp = &Rx[RxNum].SoOp[SoOpNum].Sp[SpNum];

	Sp->RxPwrD = SoOp->AntGainS + Tx[TxNum].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
	Sp->RxPwrR = SoOp->AntGainE + Tx[TxNum].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(Reflectivity);
	//Sp->ReflCoef = sqrt(Reflectivity);
	//Sp->EffReflCoef = sqrt(Sp->RxPwrR/Sp->RxPwrD);
}

void CalcSNR(long RxNum, long SoOpNum, long SpNum)
{
	struct SoOpType *SoOp;
	struct SpecularType *Sp;

	SoOp = &Rx[RxNum].SoOp[SoOpNum];
	Sp = &Rx[RxNum].SoOp[SoOpNum].Sp[SpNum];

	SoOp->NoisePwrS = Boltzmann + 10*log10(SoOp->NoiseTemS) + 10*log10(SoOp->Bandwidth);
	SoOp->NoisePwrE = Boltzmann + 10*log10(SoOp->NoiseTemE) + 10*log10(SoOp->Bandwidth);

	Sp->SnrD = SoOp->TotalRxPwrD - SoOp->NoisePwrS;
	Sp->SnrR = SoOp->TotalRxPwrR - SoOp->NoisePwrE;
}

long *CntValidSpecularPt(long RxNum, long SoOpNum, long RoiNum)
{
	long *ValidTxNum;
	long ISp, i;
	long j = 0;
	long cnt = 0;
	double theta;
	double unitR_Roi[3], unitR_Sp[3];

	ValidTxNum = calloc(Rx[RxNum].SoOp[SoOpNum].NSp,sizeof(long));

	/* ECEF to ECI */
	MTxV(World[EARTH].CWN,GroundStation[RoiNum].PosW,unitR_Roi);
	UNITV(unitR_Roi);

	for(ISp=0;ISp<Rx[RxNum].SoOp[SoOpNum].NSp;ISp++){
		if(Rx[RxNum].SoOp[SoOpNum].Sp[ISp].Exists){
			for(i=0;i<3;i++){
				unitR_Sp[i] = Rx[RxNum].SoOp[SoOpNum].Sp[ISp].PosN[i];
			}
			UNITV(unitR_Sp);
			theta = acos(VoV(unitR_Roi,unitR_Sp))*R2D;
			//printf("%f\n", theta);
//			if(theta<30){
//			if(theta<0.00078393){
			if(theta<0.0156){
				ValidTxNum[j] = Rx[RxNum].SoOp[SoOpNum].Sp[ISp].TxNum;
				cnt++;
				j++;
	            if(j == Rx[RxNum].SoOp[SoOpNum].NSp)
	            	break;
			}
		}
	}
	Rx[RxNum].SoOp[SoOpNum].NValidSp = cnt;
	return ValidTxNum;
}

void CalcValidRxPwr(long RoiNum, double Reflectivity)
{
	long IRx, ISoOp, ISp;
	long *ValidTxNum;

	for(IRx=0; IRx<NRx; IRx++){
		for(ISoOp=0; ISoOp<Rx[IRx].NSoOp; ISoOp++){
			ValidTxNum = CntValidSpecularPt(IRx, ISoOp, RoiNum);
			for(ISp=0; ISp<Rx[IRx].SoOp[ISoOp].NSp; ISp++){
				if(ValidTxNum[ISp]==0){
					// do nothing
				}
				else{
					FindPathLength(IRx, ISoOp, ISp, ValidTxNum[ISp]);
					CalcRxPwr(IRx, ISoOp, ISp, ValidTxNum[ISp], Reflectivity);
					Rx[IRx].SoOp[ISoOp].TotalRxPwrD += Rx[IRx].SoOp[ISoOp].Sp[ISp].RxPwrD;
					Rx[IRx].SoOp[ISoOp].TotalRxPwrR += Rx[IRx].SoOp[ISoOp].Sp[ISp].RxPwrR;
				}
			}
		}
	}
}

void ReportAllSNR(double reflectivity)
{
	long IRx, ISoOp, ISp;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;

	for(IRx=0; IRx<NRx; IRx++) {
		for(ISoOp=0; ISoOp<Rx[IRx].NSoOp; ISoOp++) {
			for(ISp=0; ISp<Rx[IRx].SoOp[ISoOp].NSp; ISp++) {
				Sp = &Rx[IRx].SoOp[ISoOp].Sp[ISp];
				SoOp = &Rx[IRx].SoOp[ISoOp];
				if(Sp->Exists) {
					long i;
					double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3];
					double MagVecSp2Rx, MagVecSp2Tx;

					for(i=0;i<3;i++) {
						VecTx2Rx[i] = SC[Sp->TxNum].PosN[i] - SC[IRx].PosN[i];
						VecSp2Rx[i] = SC[IRx].PosN[i] - Sp->PosN[i];
						VecSp2Tx[i] = SC[Sp->TxNum].PosN[i] - Sp->PosN[i];
					}
					MagVecSp2Rx = UNITV(VecSp2Rx);
					MagVecSp2Tx = UNITV(VecSp2Tx);
					Sp->PathLenD = MAGV(VecTx2Rx);
					Sp->PathLenR = MagVecSp2Rx + MagVecSp2Tx;

					Sp->RxPwrD = SoOp->AntGainS + Tx[Sp->TxNum].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
					Sp->RxPwrR = SoOp->AntGainE + Tx[Sp->TxNum].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(reflectivity);

					Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;
					Sp->SnrR = Sp->RxPwrR - SoOp->NoisePwrE;
				}
				else {
					Sp->PathLenD = 0.0;
					Sp->PathLenR = 0.0;
					Sp->RxPwrD = 0.0;
					Sp->RxPwrR = 0.0;
				}
			}
		}
	}
}

void Reflectometry_MUOS(long RxNum)
{
    static char First = TRUE;
    static FILE *Filename, *RxFile;
    struct SCType *RX;
	struct SoOpType *SoOp;
	double PosWR[3], PosWR_UNIT[3];
	double MAGR, lat, lon, alt;
	long ISp;

    if (First) {
        First = FALSE;
      //  Filename = FileOpen(InOutPath,"MUOS.42","w");
      //  RxFile = FileOpen(InOutPath,"Rx.42","w");
    }

	SoOp = &Rx[RxNum].SoOp[0];
	RX = &SC[RxNum];

	/* ECI to ECEF of Receiver's Position */
	MxV(World[EARTH].CWN,RX->PosN,PosWR);
	MAGR = CopyUnitV(PosWR,PosWR_UNIT);
	/* ECEF to WGS84 */
	ECEFToWGS84(PosWR, &lat, &lon, &alt);

	if(AbsTime>689775290 || AbsTime<683726400.4) {
		printf("%lf %lf %lf\n", AbsTime, lat, lon);
	//	fprintf(RxFile, "%lf %lf %lf\n", AbsTime, lat, lon);
	}
/*
	struct SpecularType *Sp;
	long *TxNum = NULL;
	long i;
	double estInit[3], PosWT[3], PosWT_UNIT[3];
	double MAGT;
	struct SCType *TX;

	TxNum = FindVisibleTx(RxNum, 0, 1, 4, TxNum, &SoOp->SpCnt);

	for(ISp=0; ISp<SoOp->SpCnt; ISp++){
		Sp = &SoOp->Sp[ISp];
		Sp->Exists = 1;
		Sp->TxNum = TxNum[ISp];
		TX = &SC[Sp->TxNum];

		// ECI to ECEF of Transmitter's Position
		MxV(World[EARTH].CWN,TX->PosN,PosWT);
		MAGT = CopyUnitV(PosWT,PosWT_UNIT);

		// Unconstrained Weighted Initial Estimate
		for(i=0;i<3;i++){
			estInit[i] = PosWR_UNIT[i]/MAGR + PosWT_UNIT[i]/MAGT;
		}

		FindSP_UD(estInit, Sp);

		double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3], AxisENU[3];
		double MagVecSp2Rx, MagVecSp2Tx, MagRx, MagSp, temp;

		for(i=0;i<3;i++) {
			VecTx2Rx[i] = TX->PosN[i] - RX->PosN[i];
			VecSp2Rx[i] = PosWR[i] - Sp->PosW[i]; // (dx,dy,dz) between Receiver and Specular Point in ECEF frame
			VecSp2Tx[i] = PosWT[i] - Sp->PosW[i]; // (dx,dy,dz) between Transmitter and Specular Point in ECEF frame
		}
		MagVecSp2Rx = UNITV(VecSp2Rx);
		MagVecSp2Tx = UNITV(VecSp2Tx);

		// ECEF(dx,dy,dz From Tx to Sp) to ENU wrt (lat, lon) of Sp
		temp = cos(Sp->lon)*VecSp2Tx[0] + sin(Sp->lon)*VecSp2Tx[1];
		AxisENU[0] = cos(Sp->lon)*VecSp2Tx[1] - sin(Sp->lon)*VecSp2Tx[0];
		AxisENU[1] = cos(Sp->lat)*VecSp2Tx[2] - sin(Sp->lat)*temp;
		AxisENU[2] = cos(Sp->lat)*temp + sin(Sp->lat)*VecSp2Tx[2];

		// ENU to AE
		temp = sqrt(AxisENU[0]*AxisENU[0] + AxisENU[1]*AxisENU[1]);
		Sp->az = atan2(AxisENU[0],AxisENU[1]);
		Sp->elev = atan2(AxisENU[2],temp);

		MagRx = sqrt(RX->PosN[0]*RX->PosN[0] + RX->PosN[1]*RX->PosN[1] + RX->PosN[2]*RX->PosN[2]);
		MagSp = sqrt(Sp->PosN[0]*Sp->PosN[0] + Sp->PosN[1]*Sp->PosN[1] + Sp->PosN[2]*Sp->PosN[2]);

		Sp->lookAng = HalfPi - Sp->elev - acos(VoV(RX->PosN, Sp->PosN)/MagRx/MagSp);

		Sp->PathLenD = MAGV(VecTx2Rx);
		Sp->PathLenR = MagVecSp2Rx + MagVecSp2Tx;

		Sp->b = sqrt(SoOp->Wavelength*alt*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev);
		Sp->a = Sp->b / sin(Sp->elev);
		//temp = sqrt(SoOp->Wavelength*MagVecSp2Rx*MagVecSp2Tx/Sp->PathLenR)/sin(Sp->elev); // Sgphaladi thesis

		Sp->RxPwrD = SoOp->AntGainS + Tx[0].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
		Sp->RxPwrR = SoOp->AntGainE + Tx[0].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(Sm->reflectivity);

		Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;
		Sp->SnrR = Sp->RxPwrR - SoOp->NoisePwrE;

		fprintf(Filename,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
				AbsTime, Sp->lat, Sp->lon, Sp->SnrD, Sp->SnrR, Sp->az, Sp->elev, Sp->a, Sp->b, Sp->lookAng);
	}
	free(TxNum);*/
}

void Reflectometry_Orbcomm(long RxNum)
{
    static char First = TRUE;
    static FILE *Filename;
    struct SCType *RX;
	struct SoOpType *SoOp;
	double PosWR[3], PosWR_UNIT[3];
	double MAGR, lat, lon, alt;
	long ISp;

    if (First) {
        First = FALSE;
        Filename = FileOpen(InOutPath,"Orbcomm.42","w");
    }

	SoOp = &Rx[RxNum].SoOp[1];
	RX = &SC[RxNum];

	/* ECI to ECEF of Receiver's Position */
	MxV(World[EARTH].CWN,RX->PosN,PosWR);
	MAGR = CopyUnitV(PosWR,PosWR_UNIT);
	/* ECEF to WGS84 */
	ECEFToWGS84(PosWR, &lat, &lon, &alt);

	if(lon > -0.523598775598299 && lon < HalfPi) {
		struct SpecularType *Sp;
		long *TxNum = NULL;
		long i;
		double estInit[3], PosWT[3], PosWT_UNIT[3];
		double MAGT;
		struct SCType *TX;

		TxNum = FindVisibleTx(RxNum, 1, 11, 51, TxNum, &SoOp->SpCnt);

		for(ISp=0; ISp<SoOp->SpCnt; ISp++){
			Sp = &SoOp->Sp[ISp];

			Sp->Exists = 1;
			Sp->TxNum = TxNum[ISp];
			TX = &SC[Sp->TxNum];

			/* ECI to ECEF of Transmitter's Position */
			MxV(World[EARTH].CWN,TX->PosN,PosWT);
			MAGT = CopyUnitV(PosWT,PosWT_UNIT);

			/* Unconstrained Weighted Initial Estimate */
			for(i=0;i<3;i++){
				estInit[i] = PosWR_UNIT[i]/MAGR + PosWT_UNIT[i]/MAGT;
			}

			FindSP_UD(estInit, Sp);

			if(Sp->lon > 0.523598775598299 && Sp->lon < 0.524497969228892) {
				double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3], AxisENU[3];
				double MagVecSp2Rx, MagVecSp2Tx, MagRx, MagSp, temp;

				for(i=0;i<3;i++) {
					VecTx2Rx[i] = TX->PosN[i] - RX->PosN[i];
					VecSp2Rx[i] = PosWR[i] - Sp->PosW[i]; /* (dx,dy,dz) between Receiver and Specular Point in ECEF frame */
					VecSp2Tx[i] = PosWT[i] - Sp->PosW[i]; /* (dx,dy,dz) between Transmitter and Specular Point in ECEF frame */
				}
				MagVecSp2Rx = UNITV(VecSp2Rx);
				MagVecSp2Tx = UNITV(VecSp2Tx);

				/*ECEF(dx,dy,dz From Tx to Sp) to ENU wrt (lat, lon) of Sp */
				temp = cos(Sp->lon)*VecSp2Tx[0] + sin(Sp->lon)*VecSp2Tx[1];
				AxisENU[0] = cos(Sp->lon)*VecSp2Tx[1] - sin(Sp->lon)*VecSp2Tx[0];
				AxisENU[1] = cos(Sp->lat)*VecSp2Tx[2] - sin(Sp->lat)*temp;
				AxisENU[2] = cos(Sp->lat)*temp + sin(Sp->lat)*VecSp2Tx[2];

				// ENU to AE
				temp = sqrt(AxisENU[0]*AxisENU[0] + AxisENU[1]*AxisENU[1]);
				Sp->az = atan2(AxisENU[0],AxisENU[1]);
				Sp->elev = atan2(AxisENU[2],temp);

				MagRx = sqrt(RX->PosN[0]*RX->PosN[0] + RX->PosN[1]*RX->PosN[1] + RX->PosN[2]*RX->PosN[2]);
				MagSp = sqrt(Sp->PosN[0]*Sp->PosN[0] + Sp->PosN[1]*Sp->PosN[1] + Sp->PosN[2]*Sp->PosN[2]);

				Sp->lookAng = HalfPi - Sp->elev - acos(VoV(RX->PosN, Sp->PosN)/MagRx/MagSp);

				Sp->PathLenD = MAGV(VecTx2Rx);
				Sp->PathLenR = MagVecSp2Rx + MagVecSp2Tx;

				Sp->b = sqrt(SoOp->Wavelength*alt*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev);
				Sp->a = Sp->b / sin(Sp->elev);
				//temp = sqrt(SoOp->Wavelength*MagVecSp2Rx*MagVecSp2Tx/Sp->PathLenR)/sin(Sp->elev); // Sgphaladi thesis

				Sp->RxPwrD = SoOp->AntGainS + Tx[0].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
				Sp->RxPwrR = SoOp->AntGainE + Tx[0].EIRP + 10*log10(SoOp->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(Sm->reflectivity);

				Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;
				Sp->SnrR = Sp->RxPwrR - SoOp->NoisePwrE;

				fprintf(Filename,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						AbsTime, Sp->lat, Sp->lon, Sp->SnrD, Sp->SnrR, Sp->az, Sp->elev, Sp->a, Sp->b, Sp->lookAng);
			}
			else {
				// do nothing
			}
		}
		free(TxNum);
	}
	else {
		for(ISp=0; ISp<SoOp->NSp; ISp++){
			SoOp->Sp[ISp].Exists = 0;
		}
	}
}

void Reflectometry(long RxNum)
{
    static char First = TRUE;
    static FILE *MuosFile, *OrbcommFile, *GpsFile;
    struct SCType *RX;
	struct SoOpType *SoOp0, *SoOp1, *SoOp2;
	double PosWR[3], PosWR_UNIT[3];
	double MAGR, lat, lon, alt;
	long ISp;

    if (First) {
        First = FALSE;
        MuosFile = FileOpen(InOutPath,"MUOS.42","w");
        GpsFile = FileOpen(InOutPath,"GPS.42","w");
        OrbcommFile = FileOpen(InOutPath,"Orbcomm.42","w");
    }

	SoOp0 = &Rx[RxNum].SoOp[0];
	SoOp1 = &Rx[RxNum].SoOp[1];
	SoOp2 = &Rx[RxNum].SoOp[2];
	RX = &SC[RxNum];

	/* ECI to ECEF of Receiver's Position */
	MxV(World[EARTH].CWN,RX->PosN,PosWR);
	MAGR = CopyUnitV(PosWR,PosWR_UNIT);
	/* ECEF to WGS84 */
	ECEFToWGS84(PosWR, &lat, &lon, &alt);

	if(AbsTime>689775290 || AbsTime<683726400.4) {
		printf("%lf %lf %lf\n", AbsTime, lat, lon);
	//	fprintf(RxFile, "%lf %lf %lf\n", AbsTime, lat, lon);
	}

//	if(lon > Grid->lonRef-1.047197551196598 && lon < Grid->lonRef+1.047197551196598) {
		struct SpecularType *Sp;
		long *TxNum = NULL;
		long i;
		double estInit[3], PosWT[3], PosWT_UNIT[3];
		double MAGT;
		struct SCType *TX;

		/* MUOS */
		TxNum = FindVisibleTx(RxNum, 0, 6, 9, TxNum, &SoOp0->SpCnt);

		for(ISp=0; ISp<SoOp0->SpCnt; ISp++){
			Sp = &SoOp0->Sp[ISp];
			Sp->Exists = 1;
			Sp->TxNum = TxNum[ISp];
			TX = &SC[Sp->TxNum];

			/* ECI to ECEF of Transmitter's Position */
			MxV(World[EARTH].CWN,TX->PosN,PosWT);
			MAGT = CopyUnitV(PosWT,PosWT_UNIT);

			/* Unconstrained Weighted Initial Estimate */
			for(i=0;i<3;i++){
				estInit[i] = PosWR_UNIT[i]/MAGR + PosWT_UNIT[i]/MAGT;
			}

			FindSP_UD(estInit, Sp);

			if((Sp->lon > -3.089232776029963 && Sp->lon < -3.080233447089754)|| /* Lon -177 */
			   (Sp->lon > -1.745329251994330 && Sp->lon < -1.736329923054121)|| /* Lon -100 */
			   (Sp->lon > -0.270526034059121 && Sp->lon < -0.261526705118912)|| /* Lon -15.5 */
			   (Sp->lon > 1.308996938995747 && Sp->lon < 1.317996267935956)|| /* Lon 75 */
			   (Sp->lon > -2.417281014012147 && Sp->lon < -2.408281685071938)|| /* Lon -138.5 */
			   (Sp->lon > -1.007927643026725 && Sp->lon < -0.998928314086516)|| /* Lon -57.75 */
			   (Sp->lon > 0.519235452468313 && Sp->lon < 0.528234781408522)|| /* Lon 29.75 */
			   (Sp->lon > 2.251474735072685 && Sp->lon < 2.260474064012894)|| /* Lon 129 */
			   (Sp->lon > 2.879793265790644 && Sp->lon < 2.888792594730853)) { /* Lon 165 */
				double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3], AxisENU[3];
				double MagVecSp2Rx, MagVecSp2Tx, MagRx, MagSp, temp;

				for(i=0;i<3;i++) {
					VecTx2Rx[i] = TX->PosN[i] - RX->PosN[i];
					VecSp2Rx[i] = PosWR[i] - Sp->PosW[i]; /* (dx,dy,dz) between Receiver and Specular Point in ECEF frame */
					VecSp2Tx[i] = PosWT[i] - Sp->PosW[i]; /* (dx,dy,dz) between Transmitter and Specular Point in ECEF frame */
				}
				MagVecSp2Rx = UNITV(VecSp2Rx);
				MagVecSp2Tx = UNITV(VecSp2Tx);

				/*ECEF(dx,dy,dz From Tx to Sp) to ENU wrt (lat, lon) of Sp */
				temp = cos(Sp->lon)*VecSp2Tx[0] + sin(Sp->lon)*VecSp2Tx[1];
				AxisENU[0] = cos(Sp->lon)*VecSp2Tx[1] - sin(Sp->lon)*VecSp2Tx[0];
				AxisENU[1] = cos(Sp->lat)*VecSp2Tx[2] - sin(Sp->lat)*temp;
				AxisENU[2] = cos(Sp->lat)*temp + sin(Sp->lat)*VecSp2Tx[2];

				// ENU to AE
				temp = sqrt(AxisENU[0]*AxisENU[0] + AxisENU[1]*AxisENU[1]);
				Sp->az = atan2(AxisENU[0],AxisENU[1]);
				Sp->elev = atan2(AxisENU[2],temp);

				MagRx = sqrt(RX->PosN[0]*RX->PosN[0] + RX->PosN[1]*RX->PosN[1] + RX->PosN[2]*RX->PosN[2]);
				MagSp = sqrt(Sp->PosN[0]*Sp->PosN[0] + Sp->PosN[1]*Sp->PosN[1] + Sp->PosN[2]*Sp->PosN[2]);

				Sp->lookAng = HalfPi - Sp->elev - acos(VoV(RX->PosN, Sp->PosN)/MagRx/MagSp);

				Sp->PathLenD = MAGV(VecTx2Rx);
				Sp->PathLenR = MagVecSp2Rx + MagVecSp2Tx;

				Sp->b = sqrt(SoOp0->Wavelength*alt*sin(Sp->elev)+SoOp0->Wavelength*SoOp0->Wavelength) / sin(Sp->elev);
				Sp->a = Sp->b / sin(Sp->elev);
				//temp = sqrt(SoOp->Wavelength*MagVecSp2Rx*MagVecSp2Tx/Sp->PathLenR)/sin(Sp->elev); // Sgphaladi thesis

				Sp->RxPwrD = SoOp0->AntGainS + Tx[0].EIRP + 10*log10(SoOp0->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
				Sp->RxPwrR = SoOp0->AntGainE + Tx[0].EIRP + 10*log10(SoOp0->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(Sm->reflectivity);

				Sp->SnrD = Sp->RxPwrD - SoOp0->NoisePwrS;
				Sp->SnrR = Sp->RxPwrR - SoOp0->NoisePwrE;

				Sp->EffReflCoef = (Sp->RxPwrR-Sp->RxPwrD)*0.5;
				Sp->VarERC = 4*Sp->EffReflCoef - Sp->SnrD - Sp->SnrR;
				Sp->VarSM = 10*log10(Sp->PathLenR*Sp->PathLenR/Sp->PathLenD/Sp->PathLenD) + SoOp0->AntGainS - SoOp0->AntGainE + Sp->VarERC;

				fprintf(MuosFile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						AbsTime, Sp->lat, Sp->lon, Sp->SnrD, Sp->SnrR, Sp->az, Sp->elev, Sp->a, Sp->b, Sp->lookAng, Sp->EffReflCoef, Sp->VarERC, Sp->VarSM);
			}
			else {
				// do nothing
			}
		}
		free(TxNum);

		/* Orbcomm */
		TxNum = FindVisibleTx(RxNum, 1, 11, 51, TxNum, &SoOp1->SpCnt);

		for(ISp=0; ISp<SoOp1->SpCnt; ISp++){
			Sp = &SoOp1->Sp[ISp];
			Sp->Exists = 1;
			Sp->TxNum = TxNum[ISp];
			TX = &SC[Sp->TxNum];

			/* ECI to ECEF of Transmitter's Position */
			MxV(World[EARTH].CWN,TX->PosN,PosWT);
			MAGT = CopyUnitV(PosWT,PosWT_UNIT);

			/* Unconstrained Weighted Initial Estimate */
			for(i=0;i<3;i++){
				estInit[i] = PosWR_UNIT[i]/MAGR + PosWT_UNIT[i]/MAGT;
			}

			FindSP_UD(estInit, Sp);

			if((Sp->lon > -3.089232776029963 && Sp->lon < -3.080233447089754)|| /* Lon -177 */
			   (Sp->lon > -1.745329251994330 && Sp->lon < -1.736329923054121)|| /* Lon -100 */
			   (Sp->lon > -0.270526034059121 && Sp->lon < -0.261526705118912)|| /* Lon -15.5 */
			   (Sp->lon > 1.308996938995747 && Sp->lon < 1.317996267935956)|| /* Lon 75 */
			   (Sp->lon > -2.417281014012147 && Sp->lon < -2.408281685071938)|| /* Lon -138.5 */
			   (Sp->lon > -1.007927643026725 && Sp->lon < -0.998928314086516)|| /* Lon -57.75 */
			   (Sp->lon > 0.519235452468313 && Sp->lon < 0.528234781408522)|| /* Lon 29.75 */
			   (Sp->lon > 2.251474735072685 && Sp->lon < 2.260474064012894)|| /* Lon 129 */
			   (Sp->lon > 2.879793265790644 && Sp->lon < 2.888792594730853)) { /* Lon 165 */
				double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3], AxisENU[3];
				double MagVecSp2Rx, MagVecSp2Tx, MagRx, MagSp, temp;

				for(i=0;i<3;i++) {
					VecTx2Rx[i] = TX->PosN[i] - RX->PosN[i];
					VecSp2Rx[i] = PosWR[i] - Sp->PosW[i]; /* (dx,dy,dz) between Receiver and Specular Point in ECEF frame */
					VecSp2Tx[i] = PosWT[i] - Sp->PosW[i]; /* (dx,dy,dz) between Transmitter and Specular Point in ECEF frame */
				}
				MagVecSp2Rx = UNITV(VecSp2Rx);
				MagVecSp2Tx = UNITV(VecSp2Tx);

				/*ECEF(dx,dy,dz From Tx to Sp) to ENU wrt (lat, lon) of Sp */
				temp = cos(Sp->lon)*VecSp2Tx[0] + sin(Sp->lon)*VecSp2Tx[1];
				AxisENU[0] = cos(Sp->lon)*VecSp2Tx[1] - sin(Sp->lon)*VecSp2Tx[0];
				AxisENU[1] = cos(Sp->lat)*VecSp2Tx[2] - sin(Sp->lat)*temp;
				AxisENU[2] = cos(Sp->lat)*temp + sin(Sp->lat)*VecSp2Tx[2];

				// ENU to AE
				temp = sqrt(AxisENU[0]*AxisENU[0] + AxisENU[1]*AxisENU[1]);
				Sp->az = atan2(AxisENU[0],AxisENU[1]);
				Sp->elev = atan2(AxisENU[2],temp);

				MagRx = sqrt(RX->PosN[0]*RX->PosN[0] + RX->PosN[1]*RX->PosN[1] + RX->PosN[2]*RX->PosN[2]);
				MagSp = sqrt(Sp->PosN[0]*Sp->PosN[0] + Sp->PosN[1]*Sp->PosN[1] + Sp->PosN[2]*Sp->PosN[2]);

				Sp->lookAng = HalfPi - Sp->elev - acos(VoV(RX->PosN, Sp->PosN)/MagRx/MagSp);

				Sp->PathLenD = MAGV(VecTx2Rx);
				Sp->PathLenR = MagVecSp2Rx + MagVecSp2Tx;

				Sp->b = sqrt(SoOp1->Wavelength*alt*sin(Sp->elev)+SoOp1->Wavelength*SoOp1->Wavelength) / sin(Sp->elev);
				Sp->a = Sp->b / sin(Sp->elev);
				//temp = sqrt(SoOp->Wavelength*MagVecSp2Rx*MagVecSp2Tx/Sp->PathLenR)/sin(Sp->elev); // Sgphaladi thesis

				Sp->RxPwrD = SoOp1->AntGainS + Tx[0].EIRP + 10*log10(SoOp1->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
				Sp->RxPwrR = SoOp1->AntGainE + Tx[0].EIRP + 10*log10(SoOp1->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(Sm->reflectivity);

				Sp->SnrD = Sp->RxPwrD - SoOp1->NoisePwrS;
				Sp->SnrR = Sp->RxPwrR - SoOp1->NoisePwrE;

				Sp->EffReflCoef = (Sp->RxPwrR-Sp->RxPwrD)*0.5;
				Sp->VarERC = 4*Sp->EffReflCoef - Sp->SnrD - Sp->SnrR;
				Sp->VarSM = 10*log10(Sp->PathLenR*Sp->PathLenR/Sp->PathLenD/Sp->PathLenD) + SoOp1->AntGainS - SoOp1->AntGainE + Sp->VarERC;

				fprintf(OrbcommFile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						AbsTime, Sp->lat, Sp->lon, Sp->SnrD, Sp->SnrR, Sp->az, Sp->elev, Sp->a, Sp->b, Sp->lookAng, Sp->EffReflCoef, Sp->VarERC, Sp->VarSM);
			}
			else {
				// do nothing
			}
		}
		free(TxNum);

		/* GPS */
		TxNum = FindVisibleTx(RxNum, 2, 52, 82, TxNum, &SoOp2->SpCnt);

		for(ISp=0; ISp<SoOp2->SpCnt; ISp++){
			Sp = &SoOp2->Sp[ISp];
			Sp->Exists = 1;
			Sp->TxNum = TxNum[ISp];
			TX = &SC[Sp->TxNum];

			/* ECI to ECEF of Transmitter's Position */
			MxV(World[EARTH].CWN,TX->PosN,PosWT);
			MAGT = CopyUnitV(PosWT,PosWT_UNIT);

			/* Unconstrained Weighted Initial Estimate */
			for(i=0;i<3;i++){
				estInit[i] = PosWR_UNIT[i]/MAGR + PosWT_UNIT[i]/MAGT;
			}

			FindSP_UD(estInit, Sp);

			if((Sp->lon > -3.089232776029963 && Sp->lon < -3.080233447089754)|| /* Lon -177 */
			   (Sp->lon > -1.745329251994330 && Sp->lon < -1.736329923054121)|| /* Lon -100 */
			   (Sp->lon > -0.270526034059121 && Sp->lon < -0.261526705118912)|| /* Lon -15.5 */
			   (Sp->lon > 1.308996938995747 && Sp->lon < 1.317996267935956)|| /* Lon 75 */
			   (Sp->lon > -2.417281014012147 && Sp->lon < -2.408281685071938)|| /* Lon -138.5 */
			   (Sp->lon > -1.007927643026725 && Sp->lon < -0.998928314086516)|| /* Lon -57.75 */
			   (Sp->lon > 0.519235452468313 && Sp->lon < 0.528234781408522)|| /* Lon 29.75 */
			   (Sp->lon > 2.251474735072685 && Sp->lon < 2.260474064012894)|| /* Lon 129 */
			   (Sp->lon > 2.879793265790644 && Sp->lon < 2.888792594730853)) { /* Lon 165 */
				double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3], AxisENU[3];
				double MagVecSp2Rx, MagVecSp2Tx, MagRx, MagSp, temp;

				for(i=0;i<3;i++) {
					VecTx2Rx[i] = TX->PosN[i] - RX->PosN[i];
					VecSp2Rx[i] = PosWR[i] - Sp->PosW[i]; /* (dx,dy,dz) between Receiver and Specular Point in ECEF frame */
					VecSp2Tx[i] = PosWT[i] - Sp->PosW[i]; /* (dx,dy,dz) between Transmitter and Specular Point in ECEF frame */
				}
				MagVecSp2Rx = UNITV(VecSp2Rx);
				MagVecSp2Tx = UNITV(VecSp2Tx);

				/*ECEF(dx,dy,dz From Tx to Sp) to ENU wrt (lat, lon) of Sp */
				temp = cos(Sp->lon)*VecSp2Tx[0] + sin(Sp->lon)*VecSp2Tx[1];
				AxisENU[0] = cos(Sp->lon)*VecSp2Tx[1] - sin(Sp->lon)*VecSp2Tx[0];
				AxisENU[1] = cos(Sp->lat)*VecSp2Tx[2] - sin(Sp->lat)*temp;
				AxisENU[2] = cos(Sp->lat)*temp + sin(Sp->lat)*VecSp2Tx[2];

				// ENU to AE
				temp = sqrt(AxisENU[0]*AxisENU[0] + AxisENU[1]*AxisENU[1]);
				Sp->az = atan2(AxisENU[0],AxisENU[1]);
				Sp->elev = atan2(AxisENU[2],temp);

				MagRx = sqrt(RX->PosN[0]*RX->PosN[0] + RX->PosN[1]*RX->PosN[1] + RX->PosN[2]*RX->PosN[2]);
				MagSp = sqrt(Sp->PosN[0]*Sp->PosN[0] + Sp->PosN[1]*Sp->PosN[1] + Sp->PosN[2]*Sp->PosN[2]);

				Sp->lookAng = HalfPi - Sp->elev - acos(VoV(RX->PosN, Sp->PosN)/MagRx/MagSp);

				Sp->PathLenD = MAGV(VecTx2Rx);
				Sp->PathLenR = MagVecSp2Rx + MagVecSp2Tx;

				Sp->b = sqrt(SoOp2->Wavelength*alt*sin(Sp->elev)+SoOp2->Wavelength*SoOp2->Wavelength) / sin(Sp->elev);
				Sp->a = Sp->b / sin(Sp->elev);
				//temp = sqrt(SoOp->Wavelength*MagVecSp2Rx*MagVecSp2Tx/Sp->PathLenR)/sin(Sp->elev); // Sgphaladi thesis

				Sp->RxPwrD = SoOp2->AntGainS + Tx[0].EIRP + 10*log10(SoOp2->Wavelength) - 20*log10(Sp->PathLenD) - 21.984197280441926;
				Sp->RxPwrR = SoOp2->AntGainE + Tx[0].EIRP + 10*log10(SoOp2->Wavelength) - 20*log10(Sp->PathLenR) - 21.984197280441926 + 10*log10(Sm->reflectivity);

				Sp->SnrD = Sp->RxPwrD - SoOp2->NoisePwrS;
				Sp->SnrR = Sp->RxPwrR - SoOp2->NoisePwrE;

				Sp->EffReflCoef = (Sp->RxPwrR-Sp->RxPwrD)*0.5;
				Sp->VarERC = 4*Sp->EffReflCoef - Sp->SnrD - Sp->SnrR;
				Sp->VarSM = 10*log10(Sp->PathLenR*Sp->PathLenR/Sp->PathLenD/Sp->PathLenD) + SoOp2->AntGainS - SoOp2->AntGainE + Sp->VarERC;

				fprintf(GpsFile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						AbsTime, Sp->lat, Sp->lon, Sp->SnrD, Sp->SnrR, Sp->az, Sp->elev, Sp->a, Sp->b, Sp->lookAng, Sp->EffReflCoef, Sp->VarERC, Sp->VarSM);
			}
			else {
				// do nothing
			}
		}
		free(TxNum);
//	}
//	else {
//		for(ISp=0; ISp<SoOp0->NSp; ISp++)
//			SoOp0->Sp[ISp].Exists = 0;
//		for(ISp=0; ISp<SoOp1->NSp; ISp++)
//			SoOp1->Sp[ISp].Exists = 0;
//		for(ISp=0; ISp<SoOp2->NSp; ISp++)
//			SoOp2->Sp[ISp].Exists = 0;
//	}
}

void ReportSoOp(void)
{
    if (OutFlag) {
    	Reflectometry_MUOS(0);
    	/*Reflectometry_MUOS(1, 0.1);
    	Reflectometry_MUOS(2, 0.1);
    	Reflectometry_MUOS(3, 0.1);
    	Reflectometry_MUOS(4, 0.1);
    	Reflectometry_MUOS(5, 0.1);
*/
    	Reflectometry_Orbcomm(0);
  /*  	Reflectometry_Orbcomm(1, 0.1);
    	Reflectometry_Orbcomm(2, 0.1);
    	Reflectometry_Orbcomm(3, 0.1);
    	Reflectometry_Orbcomm(4, 0.1);
    	Reflectometry_Orbcomm(5, 0.1);*/
    }
}

void ReportAttitude(void)
{
	  static FILE *timefile,*AbsTimeFile;
	  static FILE *PosNfile,*VelNfile;
	  static FILE *qbnfile,*wbnfile, *qbnfile2,*wbnfile2, *qbnfile3,*wbnfile3;
	  static FILE *PosWfile,*VelWfile;
	  static FILE *PosRfile,*VelRfile;
	  static FILE *Hvnfile,*KEfile;
	  static FILE *RPYfile;
	  static FILE *Hwhlfile;
	  static FILE *MTBfile;
	  static FILE *ProjAreaFile;
	  static FILE *EclipseFile;
	  static FILE *altFile;
	  static long outFlagAttitude = 0;
	  static long outFlagCnt = 0;

	  static char First = TRUE;
	  double CBL[3][3],Roll,Pitch,Yaw;
	  double PosW[3],VelW[3],PosR[3],VelR[3];

	  if (First) {
		 First = FALSE;
		 timefile = FileOpen(InOutPath,"time.42","w");
		 AbsTimeFile = FileOpen(InOutPath,"AbsTime.42","w");
		 PosNfile = FileOpen(InOutPath,"PosN.42","w");
		 VelNfile = FileOpen(InOutPath,"VelN.42","w");
		 PosWfile = FileOpen(InOutPath,"PosW.42","w");
		 VelWfile = FileOpen(InOutPath,"VelW.42","w");
		 PosRfile = FileOpen(InOutPath,"PosR.42","w");
		 VelRfile = FileOpen(InOutPath,"VelR.42","w");
		 qbnfile = FileOpen(InOutPath,"qbn.42","w");
		 wbnfile = FileOpen(InOutPath,"wbn.42","w");
		 qbnfile2 = FileOpen(InOutPath,"qbn2.42","w");
		 wbnfile2 = FileOpen(InOutPath,"wbn2.42","w");
		 qbnfile3 = FileOpen(InOutPath,"qbn3.42","w");
		 wbnfile3 = FileOpen(InOutPath,"wbn3.42","w");
		 Hvnfile = FileOpen(InOutPath,"Hvn.42","w");
		 KEfile = FileOpen(InOutPath,"KE.42","w");
		 ProjAreaFile = FileOpen(InOutPath,"ProjArea.42","w");
		 RPYfile = FileOpen(InOutPath,"RPY.42","w");
		 Hwhlfile = FileOpen(InOutPath,"Hwhl.42","w");
		 MTBfile = FileOpen(InOutPath,"MTB.42","w");
		 EclipseFile = FileOpen(InOutPath, "EclipseFlag.42","w");
		 altFile = FileOpen(InOutPath, "alt.42","w");
	  }

	  if(OutFlag) {
		  outFlagCnt++;
		  if(outFlagCnt == 10)
		  {
			  outFlagAttitude = 1;
			  outFlagCnt = 0;
		  }
	  }

	  if (outFlagAttitude) {
		  outFlagAttitude = 0;
		 fprintf(timefile,"%lf\n",SimTime);
		 fprintf(AbsTimeFile,"%lf\n",AbsTime);
			fprintf(PosNfile,"%le %le %le\n",
			   SC[0].PosN[0],SC[0].PosN[1],SC[0].PosN[2]);
			fprintf(VelNfile,"%le %le %le\n",
			   SC[0].VelN[0],SC[0].VelN[1],SC[0].VelN[2]);
			MxV(World[EARTH].CWN,SC[0].PosN,PosW);
			MxV(World[EARTH].CWN,SC[0].VelN,VelW);
			fprintf(PosWfile,"%18.12le %18.12le %18.12le\n",
			   PosW[0],PosW[1],PosW[2]);
			fprintf(VelWfile,"%18.12le %18.12le %18.12le\n",
			   VelW[0],VelW[1],VelW[2]);
			MxV(Rgn[Orb[SC[0].RefOrb].Region].CN,SC[0].PosR,PosR);
			MxV(Rgn[Orb[SC[0].RefOrb].Region].CN,SC[0].VelR,VelR);
			fprintf(PosRfile,"%le %le %le\n",
			   PosR[0],PosR[1],PosR[2]);
			fprintf(VelRfile,"%le %le %le\n",
			   VelR[0],VelR[1],VelR[2]);
			fprintf(qbnfile,"%le %le %le %le\n",
			   SC[0].B[0].qn[0],SC[0].B[0].qn[1],SC[0].B[0].qn[2],SC[0].B[0].qn[3]);
			fprintf(wbnfile,"%le %le %le\n",
			   SC[0].B[0].wn[0],SC[0].B[0].wn[1],SC[0].B[0].wn[2]);
			fprintf(Hvnfile,"%18.12le %18.12le %18.12le\n",
			   SC[0].Hvn[0],SC[0].Hvn[1],SC[0].Hvn[2]);
			fprintf(KEfile,"%18.12le\n",FindTotalKineticEnergy(&SC[0]));
	   //     fprintf(ProjAreaFile,"%18.12le %18.12le\n",
	   //        FindTotalProjectedArea(&SC[0],ZAxis),
	   //        FindTotalUnshadedProjectedArea(&SC[0],ZAxis));
			MxMT(SC[0].B[0].CN,SC[0].CLN,CBL);
			C2A(123,CBL,&Roll,&Pitch,&Yaw);
			fprintf(RPYfile,"%lf %lf %lf\n",Roll*R2D,Pitch*R2D,Yaw*R2D);
			fprintf(Hwhlfile,"%lf %lf %lf\n",SC[0].Whl[0].H,SC[0].Whl[1].H,SC[0].Whl[2].H);
			fprintf(MTBfile,"%lf %lf %lf\n",SC[0].MTB[0].M,SC[0].MTB[1].M,SC[0].MTB[2].M);

			fprintf(EclipseFile,"%ld\n", SC[0].Eclipse);
			if(SC[0].Nb>1)
			{
				fprintf(qbnfile2,"%le %le %le %le\n",
				   SC[0].B[1].qn[0],SC[0].B[1].qn[1],SC[0].B[1].qn[2],SC[0].B[1].qn[3]);
				fprintf(wbnfile2,"%le %le %le\n",
				   SC[0].B[1].wn[0],SC[0].B[1].wn[1],SC[0].B[1].wn[2]);
				fprintf(qbnfile3,"%le %le %le %le\n",
				   SC[0].B[2].qn[0],SC[0].B[2].qn[1],SC[0].B[2].qn[2],SC[0].B[2].qn[3]);
				fprintf(wbnfile3,"%le %le %le\n",
				   SC[0].B[2].wn[0],SC[0].B[2].wn[1],SC[0].B[2].wn[2]);
			}
			fprintf(altFile, "%le\n",
			   MAGV(SC[0].PosN)-World[Orb[SC[0].RefOrb].World].rad);

			CmdReport();
	  }
	  if (CleanUpFlag) {
		 fclose(timefile);
	  }

}
void RepeatGroundTrack(struct OrbitType *O)
{
	double nPeriod, nDays, nOrbits;
	double tolerance, error1, error2;
	double lon, dlon;
	double wE = 7.2921151467E-5;

	tolerance = 0.1*D2R;
	nOrbits = 0;
	nDays = 0;
	error1 = 1*D2R;
	error2 = 1*D2R;
	/* ECI to ECEF of Receiver's Position */
//	MxV(World[EARTH].CWN,S->PosN,PosWR);
	/* ECEF to WGS84 */
//	ECEFToWGS84(PosWR, &lat, &lon, &alt);

	lon = 0;

	nPeriod = TwoPi / (O->MeanMotion + O->ArgPdot);
	dlon = nPeriod * (wE - O->RAANdot);

	printf("dlon: %lf rad, %lf deg\n", dlon, dlon*R2D);

	while(1)
	{
		lon = lon + dlon;
		if(lon >= TwoPi)
			lon = lon - TwoPi;
		if(lon < 0)
			lon = -lon;
		nOrbits = nOrbits + 1;
		if(lon <= tolerance)
			break;
	}
	nDays = nOrbits*nPeriod;
	printf("nDays: %lf days = %lf sec \n", nDays/86400, nDays);
	printf("nOrbits: %lf\n", nOrbits);
	printf("nPeriod: %lf min\n", nPeriod/60);
	printf("tol: %lf\n", tolerance);
	printf("lon: %lf\n", lon);

//		AdvanceTime();
//		EnckeRK4(S);

//		O->RAAN = O->RAAN0 + O->RAANdot*(SimTime - 0.5/O->MeanMotion*sin(2.0*O->ArgP+2.0*O->anom));
//		O->ArgP = O->ArgP0 + O->ArgPdot*SimTime;

//		Eph2RV(O->MuPlusJ2,O->SLR,O->ecc,
//				O->inc,O->RAAN,O->ArgP,
//				AbsTime+DTSIM-O->tp,
//				O->PosN,O->VelN,&O->anom);

//		Ephemerides_Earth();
//		GravPertForce(S);

		/* ECI to ECEF of Receiver's Position */
//		MxV(World[EARTH].CWN,S->PosN,PosWR);
		/* ECEF to WGS84 */
//		ECEFToWGS84(PosWR, &lat, &lon, &alt);
}

void InitFOVs(void)
{
      FILE *infile;
      char junk[120],newline;
      char response[120],response1[120],response2[120];
      double Ang1,Ang2,Ang3;
      long Seq;
      long i;

      infile = FileOpen(InOutPath,"Inp_FOV.txt","r");
      fscanf(infile,"%[^\n] %[\n]",junk,&newline);
      fscanf(infile,"%[^\n] %[\n]",junk,&newline);
      fscanf(infile,"%ld %[^\n] %[\n]",&Nfov,junk,&newline);
      FOV = (struct FovType *) calloc(Nfov,sizeof(struct FovType));
      for(i=0;i<Nfov;i++) {
         fscanf(infile,"%[^\n] %[\n]",junk,&newline);
         fscanf(infile,"\"%[^\"]\" %[^\n] %[\n]",
            FOV[i].Label,junk,&newline);
         fscanf(infile,"%ld %lf %[^\n] %[\n]",
            &FOV[i].Nv,&FOV[i].Length,junk,&newline);
         fscanf(infile,"%lf %lf %[^\n] %[\n]",
            &FOV[i].Width,&FOV[i].Height,junk,&newline);
         if (FOV[i].Width >= 180.0) {
            printf("FOV[%ld] Width >= 180 deg.  This is not allowed.  Bailing out.\n",i);
            exit(1);
         }
         if (FOV[i].Height >= 180.0) {
            printf("FOV[%ld] Width >= 180 deg.  This is not allowed.  Bailing out.\n",i);
            exit(1);
         }
         FOV[i].Width *= D2R;
         FOV[i].Height *= D2R;
         fscanf(infile,"%f %f %f %f %[^\n] %[\n]",
            &FOV[i].Color[0],&FOV[i].Color[1],&FOV[i].Color[2],&FOV[i].Color[3],
            junk,&newline);
         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
         FOV[i].Type = DecodeString(response);
         fscanf(infile,"%s %s %[^\n] %[\n]",
            response1,response2,junk,&newline);
         FOV[i].NearExists = DecodeString(response1);
         FOV[i].FarExists = DecodeString(response2);
         fscanf(infile,"%ld %ld %[^\n] %[\n]",
            &FOV[i].SC,&FOV[i].Body,junk,&newline);
         if (FOV[i].SC >= Nsc) {
            printf("FOV[%ld].SC is out of range.\n",i);
            exit(1);
         }
         if (SC[FOV[i].SC].Exists && FOV[i].Body >= SC[FOV[i].SC].Nb) {
            printf("FOV[%ld].Body is out of range.\n",i);
            exit(1);
         }
         FOV[i].RefOrb = SC[FOV[i].SC].RefOrb;
         if (!SC[FOV[i].SC].Exists) {
            FOV[i].NearExists = FALSE;
            FOV[i].FarExists = FALSE;
         }

         fscanf(infile,"%lf %lf %lf %[^\n] %[\n]",
            &FOV[i].pb[0],&FOV[i].pb[1],&FOV[i].pb[2],junk,&newline);
         fscanf(infile,"%lf %lf %lf %ld %[^\n] %[\n]",
            &Ang1,&Ang2,&Ang3,&Seq,junk,&newline);
            A2C(Seq,Ang1*D2R,Ang2*D2R,Ang3*D2R,FOV[i].CB);
      }
      fclose(infile);
}

void InitSoOp(void)
{
    FILE *infile;
    char junk[120],newline;
	long IRx, ITx, ISoOp, i;

/* .. Read from file Inp_SoOp.txt */
    infile=FileOpen(InOutPath,"Inp_SoOp.txt","r");

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

/* .. 1km Grid .. */
    Grid = (struct GridType *) calloc(1,sizeof(struct GridType));
	fscanf(infile,"%lf %[^\n] %[\n]",&Grid->lonRef,junk,&newline);
	fscanf(infile,"%lf %[^\n] %[\n]",&Grid->latMax,junk,&newline);
	fscanf(infile,"%lf %[^\n] %[\n]",&Grid->latMin,junk,&newline);
	Grid->lonRef *= D2R;
	Grid->latMax *= D2R;
	Grid->latMin *= D2R;
	fscanf(infile,"%[^\n] %[\n]",junk,&newline);

/* .. Reflectivity of Soil Moisture */
    Sm = (struct SoilType *) calloc(1,sizeof(struct SoilType));
	fscanf(infile,"%lf %[^\n] %[\n]",&Sm->reflectivity,junk,&newline);

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
/* .. Number of Receivers */
    fscanf(infile,"%ld %[^\n] %[\n]",&NRx,junk,&newline);
    Rx = (struct RxType *) calloc(NRx,sizeof(struct RxType));
	if (Rx == NULL) {
	 printf("SoOp Rx calloc returned null pointer.  Bailing out!\n");
	 exit(1);
	}
/* .. Initialization for Receivers */
    for (IRx=0;IRx<NRx;IRx++){
        fscanf(infile,"%[^\n] %[\n]",junk,&newline);
        /* .. Assign Spacecraft .. */
        fscanf(infile,"%ld %[^\n] %[\n]",&Rx[IRx].SC,junk,&newline);
        /* .. Number of SoOp Configurations .. */
   		fscanf(infile,"%ld %[^\n] %[\n]",&Rx[IRx].NSoOp,junk,&newline);
   	    Rx[IRx].SoOp = (struct SoOpType *) calloc(Rx[IRx].NSoOp,sizeof(struct SoOpType));
        if (Rx[IRx].SoOp == NULL) {
           printf("Rx[IRx].SoOp calloc returned null pointer.  Bailing out!\n");
           exit(1);
        }
        for (ISoOp=0;ISoOp<Rx[IRx].NSoOp;ISoOp++){
            fscanf(infile,"%[^\n] %[\n]",junk,&newline);
        	/* .. Number of Specular Points */
       		fscanf(infile,"%ld %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].NSp,junk,&newline);
       	    Rx[IRx].SoOp[ISoOp].Sp = (struct SpecularType *) calloc(Rx[IRx].SoOp[ISoOp].NSp,sizeof(struct SpecularType));
            if (Rx[IRx].SoOp[ISoOp].Sp == NULL) {
               printf("Rx[IRx].SoOp[ISoOp].Sp calloc returned null pointer.  Bailing out!\n");
               exit(1);
            }
       		/* .. Sprite File Name for Specular Point .. */
			fscanf(infile,"%s %[^\n] %[\n]",Rx[IRx].SoOp[ISoOp].SpriteFileName,junk,&newline);
       		/* .. Rx Frequency .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].Freq,junk,&newline);
       		/* .. Bandwidth of Channel .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].Bandwidth,junk,&newline);
			/* .. Sky-view Antenna Gain for Direct Signal .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].AntGainS,junk,&newline);
			/* .. Earth-view Antenna Gain for Reflected Signal .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].AntGainE,junk,&newline);
       		/* .. Noise Temperature of Sky-view Channel .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].NoiseTemS,junk,&newline);
       		/* .. Noise Temperature of Earth-view Channel .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].NoiseTemE,junk,&newline);

			Rx[IRx].SoOp[ISoOp].Bandwidth *= 1000; // [kHz] to [Hz]
			Rx[IRx].SoOp[ISoOp].Wavelength = Clight / Rx[IRx].SoOp[ISoOp].Freq;
			Rx[IRx].SoOp[ISoOp].NoisePwrS = 10*log10(Rx[IRx].SoOp[ISoOp].NoiseTemS) + 10*log10(Rx[IRx].SoOp[ISoOp].Bandwidth) - 228.6;
			Rx[IRx].SoOp[ISoOp].NoisePwrE = 10*log10(Rx[IRx].SoOp[ISoOp].NoiseTemE) + 10*log10(Rx[IRx].SoOp[ISoOp].Bandwidth) - 228.6;
        }
    }

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
/* .. Number of Transmitters */
    fscanf(infile,"%ld %[^\n] %[\n]",&NTx,junk,&newline);
    Tx = (struct TxType *) calloc(NTx,sizeof(struct TxType));
	if (Tx == NULL) {
	 printf("SoOp Tx calloc returned null pointer.  Bailing out!\n");
	 exit(1);
	}
/* .. Initialization for Transmitters */
	for (ITx=0;ITx<NTx;ITx++){
		/* .. Tx Frequency .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&Tx[ITx].Freq,junk,&newline);
		/* .. EIRP .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&Tx[ITx].EIRP,junk,&newline);
	}
	fclose(infile);
}

void RxFsw(struct SCType *S)
{
	struct AcType *AC;
	struct AcCfsCtrlType *C;
	struct AcJointType *G;
	double Hb[3],HxB[3];
	long i;

	AC = &S->AC;
	C = &AC->CfsCtrl;
	G = &AC->G[0];

	if (C->Init) {
		C->Init = 0;
		for(i=0;i<3;i++) FindPDGains(AC->MOI[i][i],0.1,0.7,&C->Kr[i],&C->Kp[i]);
		C->Kunl = 1.0E6;
		FindPDGains(100.0,0.2,1.0,&G->AngRateGain[0],&G->AngGain[0]);
		G->MaxAngRate[0] = 1.0*D2R;
		G->MaxTrq[0] = 10.0;
	}

	/* .. Sensor Processing */
	GyroProcessing(AC);
	MagnetometerProcessing(AC);
	CssProcessing(AC);
	FssProcessing(AC);
	//StarTrackerProcessing(AC);
	GpsProcessing(AC);

/* .. Commanded Attitude */
	/* Find qrn, wrn and joint angle commands */
    ThreeAxisAttitudeCommand(S);

    /* .. Attitude Control */
    QxQT(AC->qbn,AC->Cmd.qrn,AC->qbr);
    RECTIFYQ(AC->qbr);
	for(i=0;i<3;i++) {
		C->therr[i] = Limit(2.0*AC->qbr[i],-0.05,0.05);
		C->werr[i] = AC->wbn[i] - AC->Cmd.wrn[i];
		AC->Tcmd[i] = Limit(-C->Kr[i]*C->werr[i] - C->Kp[i]*C->therr[i],-0.1,0.1);
	}

	/* .. Momentum Management */
	for(i=0;i<3;i++) Hb[i] = AC->MOI[i][i]*AC->wbn[i] + AC->Whl[i].H;
	VxV(Hb,AC->bvb,HxB);
	for(i=0;i<3;i++) AC->Mcmd[i] = C->Kunl*HxB[i];

	/* .. Solar Array Steering */
//	G->Cmd.Ang[0] = atan2(AC->svb[0],AC->svb[2]);

	/* .. Actuator Processing */
	WheelProcessing(AC);
	MtbProcessing(AC);
}
