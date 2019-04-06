
#ifndef __SoOp_H
#define __SoOp_H__
/*
** #ifdef __cplusplus
** namespace _42 {
** using namespace Kit;
** #endif
*/

#define Clight 299.792458 // speed of light *10^6 [m/s]
#define Boltzmann -228.6 // Boltzmann's constant [dBW/K/Hz]

#define SC_RX1 0
#define SC_RX2 1
#define SC_RX3 2
#define SC_RX4 3
#define SC_RX5 4
#define SC_RX6 5

#define SC_MUOS1 1
#define SC_MUOS2 2
#define SC_MUOS3 3
#define SC_MUOS4 4
#define SC_MUOS5 10

#define SC_FM4 11
#define SC_FM5 12
#define SC_FM6 13
#define SC_FM7 14
#define SC_FM8 15
#define SC_FM9 16
#define SC_FM10 17
#define SC_FM11 18
#define SC_FM12 19
#define SC_FM13 20
#define SC_FM14 21
#define SC_FM15 22
#define SC_FM16 23
#define SC_FM18 24
#define SC_FM19 25
#define SC_FM20 26
#define SC_FM21 27
#define SC_FM23 28
#define SC_FM27 29
#define SC_FM30 30
#define SC_FM31 31
#define SC_FM32 32
#define SC_FM34 33
#define SC_FM35 34
#define SC_FM36 35
#define SC_FM103 36
#define SC_FM104 37
#define SC_FM105 38
#define SC_FM106 39
#define SC_FM107 40
#define SC_FM108 41
#define SC_FM109 42
#define SC_FM110 43
#define SC_FM112 44
#define SC_FM113 45
#define SC_FM114 46
#define SC_FM115 47
#define SC_FM116 48
#define SC_FM117 49
#define SC_FM118 50
#define SC_FM119 51

#define SC_PRN1 52
#define SC_PRN2 53
#define SC_PRN3 54
#define SC_PRN5 55
#define SC_PRN6 56
#define SC_PRN7 57
#define SC_PRN8 58
#define SC_PRN9 59
#define SC_PRN10 60
#define SC_PRN11 61
#define SC_PRN12 62
#define SC_PRN13 63
#define SC_PRN14 64
#define SC_PRN15 65
#define SC_PRN16 66
#define SC_PRN17 67
#define SC_PRN18 68
#define SC_PRN19 69
#define SC_PRN20 70
#define SC_PRN21 71
#define SC_PRN22 72
#define SC_PRN23 73
#define SC_PRN24 74
#define SC_PRN25 75
#define SC_PRN26 76
#define SC_PRN27 77
#define SC_PRN28 78
#define SC_PRN29 79
#define SC_PRN30 80
#define SC_PRN31 81
#define SC_PRN32 82

struct SpecularType{
	long Exists;
	long TxNum; // Transmitter Number
	double PosN[3]; // Position of Specular Point expressed in ECI frame [m]
	double PosW[3]; // Position of Specular Point expressed in ECEF frame [m]
	double lat; // latitude of specular point expressed in WGS84 [rad]
	double lon; // longitude of specular point expressed in WGS84 [rad]
	double alt; // altitude of specular point expressed in WGS84 [m]
	double az; // Azimuth angle [rad]
	double elev; // Elevation angle [rad]
//	double slantRng; // Slant Range [m]
	double lookAng; // Antenna Look Angle [rad]
	double a; // 1st Fresnel Zone Semi-major Axis [m]
	double b; // 1st Fresnel Zone Semi-minor Axis [m]
	long cnt; // the number of iteration to find specular point
	double PathLenD; // Direct Path Length [m]
	double PathLenR; // Specular Reflection Path Length [m]
	double RxPwrD; // Received Signal Power from Direct Signal [dBW]
	double RxPwrR; // Received Signal Power from Reflected Signal [dBW]
	double SnrD; // SNR of Direct Signal [dBW]
	double SnrR; // SNR of Reflected Signal [dBW]
	double EffReflCoef; // Effective Reflection Coefficient [dB]
	double VarSM; // Variance of Soil Moisture Estimate [dB]
	double VarERC; // Variance of Effective Reflection Coefficient [dB]
};

struct SoOpType{
	long NSp; // The Maximum Number of Specular Points
	long SpCnt; // The Number of Specular points
	long NValidSp; // Number of Valid Specular Points
	struct SpecularType *Sp;
	char SpriteFileName[40];
    unsigned int SpriteTexTag;
	double Freq; // [MHz]
	double Wavelength; // [m]
	double Bandwidth; // [Hz]
	double AntGainS; // Antenna Gain of Sky-view Antenna [dB]
	double AntGainE; // Antenna Gain of Earth-view Antenna [dB]
	double NoiseTemS; // Noise Temperature of Sky-view Antenna [K]
	double NoiseTemE; // Noise Temperature of Earth-view Antenna [K]
	double NoisePwrS; // Noise Power arriving at Sky-view Antenna [dBW]
	double NoisePwrE; // Noise Power arriving at Earth-view Antenna [dBW]
	double TotalRxPwrD; // Total Power of Received Direct Signals [dB]
	double TotalRxPwrR; // Total Power of Received Reflected Signals [dB]
};

struct RxType{
	long SC; // Assigned Spacecraft
    long NSoOp;   /* Number of SoOp1 Configurations */
	struct SoOpType *SoOp;
	double SwathWidth; // [m]
	double SwathArea; // [m^2]
	double corIntTime; // Coherent Integration Time [sec]
};

struct TxType{
	double Freq; // [MHz]
	double Wavelength; // [m]
	double EIRP; // [EIRP]
	double Bandwidth; // [kHz]
	double NChannel; // the number of channels
};

struct SoilType{
	double reflectivity; // Gamma
};

struct GridType{
	double lonRef;
	double latMax;
	double latMin;
};
/* Functions for Signals of Opportunity */
long FindSP_MPL(struct SCType *SCT, struct SCType *SCR, struct SpecularType *Sp);
void FindSP_UD(double estInit[3], struct SpecularType *Sp);
void FindSwathArea(long RxNum);
long *FindVisibleTx(long RxNum, long SoOpNum, long ITxStart, long ITxEnd, long *TxNum, long *TxCnt);
long FindClosestTx(long RxNum, long ITxStart, long ITxEnd);
long *FindBistaticGeo(long RxNum, long SoOpNum, long TxStart, long TxEnd);
void FindPathLength(long RxNum, long SoOpNum, long SpNum, long TxNum);
void CalcRxPwr(long RxNum, long SoOpNum, long SpNum, long TxNum, double Reflectivity);
long *CntValidSpecularPt(long RxNum, long SoOpNum, long RoiNum);
void CalcValidRxPwr(long RoiNum, double Reflectivity);

void ReportSoOp(void);
void ReportAttitude(void);

void Reflectometry_MUOS(long RxNum);
void Reflectometry_Orbcomm(long RxNum);
void Reflectometry(long RxNum);

void RepeatGroundTrack(struct OrbitType *O);

void InitFOVs(void);
void InitSoOp(void);

void RxFsw(struct SCType *S);

void GyroProcessing(struct AcType *AC);
void MagnetometerProcessing(struct AcType *AC);
void CssProcessing(struct AcType *AC);
void FssProcessing(struct AcType *AC);
void StarTrackerProcessing(struct AcType *AC);
void GpsProcessing(struct AcType *AC);
void AccelProcessing(struct AcType *AC);
void WheelProcessing(struct AcType *AC);
void MtbProcessing(struct AcType *AC);
void ThreeAxisAttitudeCommand(struct SCType *S);
/////////////////////////////////////////////
/*
** #ifdef __cplusplus
** }
** #endif
*/
#endif /* __SoOp_H__ */
