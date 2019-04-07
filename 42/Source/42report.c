/*    This file is distributed with 42,                               */
/*    the (mostly harmless) spacecraft dynamics simulation            */
/*    created by Eric Stoneking of NASA Goddard Space Flight Center   */

/*    Copyright 2010 United States Government                         */
/*    as represented by the Administrator                             */
/*    of the National Aeronautics and Space Administration.           */

/*    No copyright is claimed in the United States                    */
/*    under Title 17, U.S. Code.                                      */

/*    All Other Rights Reserved.                                      */


#include "42.h"

/* #ifdef __cplusplus
** namespace _42 {
** using namespace Kit;
** #endif
*/

/*********************************************************************/
double FindTotalProjectedArea(struct SCType *S,double VecN[3])
{
      struct BodyType *B;
      struct GeomType *G;
      struct PolyType *P;
      double ProjArea = 0.0;
      double VecB[3],VoN;
      long Ib,Ipoly;

      for(Ib=0;Ib<S->Nb;Ib++) {
         B = &S->B[Ib];

         /* Transform Direction Vector from N to B */
         MxV(B->CN,VecN,VecB);

         G = &Geom[B->GeomTag];
         for(Ipoly=0;Ipoly<G->Npoly;Ipoly++) {
            P = &G->Poly[Ipoly];
            VoN = VoV(VecB,P->Norm);
            if (VoN > 0.0) {
               ProjArea += VoN*P->Area;
            }
         }
      }
      return(ProjArea);
}
/*********************************************************************/
double FindTotalUnshadedProjectedArea(struct SCType *S,double VecN[3])
{
      struct BodyType *B;
      struct GeomType *G;
      struct PolyType *P;
      double ProjArea = 0.0;
      double VecB[3],VoN;
      long Ib,Ipoly;

      FindUnshadedAreas(S,VecN);

      for(Ib=0;Ib<S->Nb;Ib++) {
         B = &S->B[Ib];

         /* Transform Direction Vector from N to B */
         MxV(B->CN,VecN,VecB);

         G = &Geom[B->GeomTag];
         for(Ipoly=0;Ipoly<G->Npoly;Ipoly++) {
            P = &G->Poly[Ipoly];
            VoN = VoV(VecB,P->Norm);
            if (VoN > 0.0) {
               ProjArea += VoN*P->UnshadedArea;
            }
         }
      }
      return(ProjArea);
}
/*********************************************************************/
void PotatoReport(void)
{
      static FILE *Bodywnfile,*Bodyqnfile;
      static FILE *Bodyvnfile,*Bodypnfile;
      struct SCType *S;
      struct BodyType *B;
      long Ib;
      static long First = 1;

      if (First) {
         First = 0;
         Bodywnfile = FileOpen(InOutPath,"Bodywn.42","w");
         Bodyqnfile = FileOpen(InOutPath,"Bodyqn.42","w");
         Bodyvnfile = FileOpen(InOutPath,"Bodyvn.42","w");
         Bodypnfile = FileOpen(InOutPath,"Bodypn.42","w");
      }

      if (OutFlag) {
         S = &SC[0];
         for(Ib=0;Ib<S->Nb;Ib++) {
            B = &S->B[Ib];
            fprintf(Bodywnfile,"%lf %lf %lf ",
               B->wn[0],B->wn[1],B->wn[2]);
            fprintf(Bodyqnfile,"%lf %lf %lf %lf ",
               B->qn[0],B->qn[1],B->qn[2],B->qn[3]);
            fprintf(Bodyvnfile,"%lf %lf %lf ",
               B->vn[0],B->vn[1],B->vn[2]);
            fprintf(Bodypnfile,"%lf %lf %lf ",
               B->pn[0],B->pn[1],B->pn[2]);
         }
         fprintf(Bodywnfile,"\n");
         fprintf(Bodyqnfile,"\n");
         fprintf(Bodyvnfile,"\n");
         fprintf(Bodypnfile,"\n");
      }
}
/*********************************************************************/
void MagReport(void)
{
      static FILE *magfile;
      static long First = 1;
      
      if (First) {
         First = 0;
         magfile = FileOpen(InOutPath,"MagBVB.42","wt");
      }
      
      fprintf(magfile,"%le %le %le %le %le %le %le %le %le \n",
         SC[0].bvb[0],SC[0].bvb[1],SC[0].bvb[2],
         SC[0].MAG[0].Field,SC[0].MAG[1].Field,SC[0].MAG[2].Field,
         SC[0].AC.bvb[0],SC[0].AC.bvb[1],SC[0].AC.bvb[2]);
      
}
/*********************************************************************/
void GyroReport(void)
{
      static FILE *gyrofile;
      static long First = 1;
      
      if (First) {
         First = 0;
         gyrofile = FileOpen(InOutPath,"Gyro.42","wt");
      }
      
      fprintf(gyrofile,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le \n",
         SC[0].B[0].wn[0],SC[0].B[0].wn[1],SC[0].B[0].wn[2],
         SC[0].Gyro[0].TrueRate,SC[0].Gyro[1].TrueRate,SC[0].Gyro[2].TrueRate,
         SC[0].Gyro[0].Bias,SC[0].Gyro[1].Bias,SC[0].Gyro[2].Bias,
         SC[0].Gyro[0].Angle,SC[0].Gyro[1].Angle,SC[0].Gyro[2].Angle,
         SC[0].Gyro[0].MeasRate,SC[0].Gyro[1].MeasRate,SC[0].Gyro[2].MeasRate,
         SC[0].AC.wbn[0],SC[0].AC.wbn[1],SC[0].AC.wbn[2]);
      
}
/*********************************************************************/
// add new function ////
void CmdReport(void)
{
    static FILE *cmdfile, *cmd2file, *WhlCmdFile, *MtqCmdFile, *acqbnFile;
    static long First = 1;

    if (First) {
       First = 0;
       cmdfile = FileOpen(InOutPath,"Cmd.42","wt");
       cmd2file = FileOpen(InOutPath,"Cmd2.42","wt");
       acqbnFile = FileOpen(InOutPath,"acqbn.42","wt");
       WhlCmdFile = FileOpen(InOutPath, "WhlCmd.42", "wt");
       MtqCmdFile = FileOpen(InOutPath, "MtqCmd.42", "wt");
    }
    fprintf(acqbnFile,"%le %le %le %le\n",
       SC[0].AC.qbn[0], SC[0].AC.qbn[1], SC[0].AC.qbn[2], SC[0].AC.qbn[3]);

    fprintf(cmdfile,"%le %le %le %le\n",
       SC[0].AC.Cmd.qrn[0], SC[0].AC.Cmd.qrn[1], SC[0].AC.Cmd.qrn[2], SC[0].AC.Cmd.qrn[3]);
    fprintf(cmd2file,"%le %le %le \n",
       SC[0].AC.Cmd.wrn[0], SC[0].AC.Cmd.wrn[1], SC[0].AC.Cmd.wrn[2]);
    fprintf(WhlCmdFile,"%le %le %le \n",
       SC[0].AC.Tcmd[0], SC[0].AC.Tcmd[1], SC[0].AC.Tcmd[2]);
    fprintf(MtqCmdFile,"%le %le %le \n",
       SC[0].AC.Mcmd[0], SC[0].AC.Mcmd[1], SC[0].AC.Mcmd[2]);
}

void WGS84Report(double PosWT1[3], double PosWT2[3], double PosWR[3])
{
	double LatT1, LngT1, AltT1;
	double LatT2, LngT2, AltT2;
	double LatR, LngR, AltR;

    static FILE *WGS84T1file, *WGS84T2file;
    static FILE *WGS84Rfile;
    static long First = 1;

    /* Find Lng, Lat, Alt from PosW */
    ECEFToWGS84(PosWT1,&LatT1,&LngT1,&AltT1);
    ECEFToWGS84(PosWT2,&LatT2,&LngT2,&AltT2);
    ECEFToWGS84(PosWR,&LatR,&LngR,&AltR);

    if (First) {
       First = 0;
       WGS84T1file = FileOpen(InOutPath,"WGS84T1.42","wt");
       WGS84T2file = FileOpen(InOutPath,"WGS84T2.42","wt");
       WGS84Rfile = FileOpen(InOutPath,"WGS84R.42","wt");
    }

    fprintf(WGS84T1file,"%le %le %le \n",
    		 LatT1, LngT1, AltT1);
    fprintf(WGS84T2file,"%le %le %le \n",
    		 LatT2, LngT2, AltT2);
    fprintf(WGS84Rfile,"%le %le %le \n",
    		 LatR, LngR, AltR);
}
////////////////////////


void Report(void)
{
      static FILE *timefile,*AbsTimeFile;
      static FILE **xfile, **ufile, **xffile, **uffile;
      static FILE **ConstraintFile;
      static FILE *PosNfile,*VelNfile,*qbnfile,*wbnfile /*added*/, *qbnfile2, *wbnfile2, *qbnfile3, *wbnfile3/**/;
      static FILE *PosWfile,*VelWfile, /*added*/ *PosNT1file, *PosNT2file, *PosWT1file, *PosWT2file; ///
      static FILE *PosRfile,*VelRfile;
      static FILE *Hvnfile,*KEfile;
      static FILE *RPYfile;
      static FILE *Hwhlfile;
      static FILE *MTBfile;
      static FILE *ProjAreaFile;
////////////////////////////////////////////////// added ///////////////////////////////////////////////
      static FILE *EclipseFile;
      static FILE *altFile;
      static FILE *pathLenFile1;
      static FILE *pathLenFile2;
      struct SCType *SCT1, *SCT2;
      struct SpecularType *Sp1, *Sp2;
      static FILE *SP1File, *SP2File;
      double PosWT1[3], PosWT2[3];
      static FILE *SwathFile, *ClosestTxFile;
      long *TxNum1, *TxNum2;
      long *TxNum;
      long ISp;
      long *ValidSpNum;
      double zerovec[3] = {0.0, 0.0, 0.0};
//////////////////////////////////////////////////////////////////////////////////////////////////////
      static char First = TRUE;
      long Isc,i;
      struct DynType *D;
      double CBL[3][3],Roll,Pitch,Yaw;
      double PosW[3],VelW[3],PosR[3],VelR[3];
      char s[40];
   //   double ZAxis[3] = {0.0,0.0,1.0};

      if (First) {
         First = FALSE;
         timefile = FileOpen(InOutPath,"time.42","w");
         AbsTimeFile = FileOpen(InOutPath,"AbsTime.42","w");

         ufile = (FILE **) calloc(Nsc,sizeof(FILE *));
         xfile = (FILE **) calloc(Nsc,sizeof(FILE *));
         uffile = (FILE **) calloc(Nsc,sizeof(FILE *));
         xffile = (FILE **) calloc(Nsc,sizeof(FILE *));
         ConstraintFile = (FILE **) calloc(Nsc,sizeof(FILE *));
         for(Isc=0;Isc<Nsc;Isc++) {
            if (SC[Isc].Exists) {
               sprintf(s,"u%02li.42",Isc);
               ufile[Isc] = FileOpen(InOutPath,s,"w");
               sprintf(s,"x%02li.42",Isc);
               xfile[Isc] = FileOpen(InOutPath,s,"w");
               if (SC[Isc].FlexActive) {
                  sprintf(s,"uf%02li.42",Isc);
                  uffile[Isc] = FileOpen(InOutPath,s,"w");
                  sprintf(s,"xf%02li.42",Isc);
                  xffile[Isc] = FileOpen(InOutPath,s,"w");
               }
               if (SC[Isc].ConstraintsRequested) {
                  sprintf(s,"Constraint%02li.42",Isc);
                  ConstraintFile[Isc] = FileOpen(InOutPath,s,"w");
               }
            }
         }
         PosNfile = FileOpen(InOutPath,"PosN.42","w");
         PosNT1file = FileOpen(InOutPath,"PosNT1.42","w");
         PosNT2file = FileOpen(InOutPath,"PosNT2.42","w");
         VelNfile = FileOpen(InOutPath,"VelN.42","w");
         PosWfile = FileOpen(InOutPath,"PosW.42","w");
         VelWfile = FileOpen(InOutPath,"VelW.42","w");
         PosRfile = FileOpen(InOutPath,"PosR.42","w");
         VelRfile = FileOpen(InOutPath,"VelR.42","w");
         qbnfile = FileOpen(InOutPath,"qbn.42","w");
         wbnfile = FileOpen(InOutPath,"wbn.42","w");
         Hvnfile = FileOpen(InOutPath,"Hvn.42","w");
         KEfile = FileOpen(InOutPath,"KE.42","w");
         ProjAreaFile = FileOpen(InOutPath,"ProjArea.42","w");
         RPYfile = FileOpen(InOutPath,"RPY.42","w");
         Hwhlfile = FileOpen(InOutPath,"Hwhl.42","w");
         MTBfile = FileOpen(InOutPath,"MTB.42","w");
//////////////////////////////////////////////////////////////// Added//////////////////////////////////////////////////////
         qbnfile2 = FileOpen(InOutPath,"qbn2.42","w");
         wbnfile2 = FileOpen(InOutPath,"wbn2.42","w");
         qbnfile3 = FileOpen(InOutPath,"qbn3.42","w");
         wbnfile3 = FileOpen(InOutPath,"wbn3.42","w");
         EclipseFile = FileOpen(InOutPath, "EclipseFlag.42","w");
         altFile = FileOpen(InOutPath, "alt.42","w");
         SP1File = FileOpen(InOutPath, "SP1.42", "w");
         SP2File = FileOpen(InOutPath, "SP2.42", "w");
         PosWT1file = FileOpen(InOutPath,"PosWT1.42","w");
         PosWT2file = FileOpen(InOutPath,"PosWT2.42","w");
         SwathFile = FileOpen(InOutPath,"Swath.42","w");
  //       ClosestTxFile = FileOpen(InOutPath,"ClosestTx.42","w");
         pathLenFile1 = FileOpen(InOutPath,"PathLength1.42","w");
         pathLenFile2 = FileOpen(InOutPath,"PathLength2.42","w");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
      if (OutFlag) {
         fprintf(timefile,"%lf\n",SimTime);
         fprintf(AbsTimeFile,"%lf\n",AbsTime);
         for(Isc=0;Isc<Nsc;Isc++) {
            if (SC[Isc].Exists) {
               D = &SC[Isc].Dyn;
               for(i=0;i<D->Nu;i++) fprintf(ufile[Isc],"% le ",D->u[i]);
               fprintf(ufile[Isc],"\n");
               for(i=0;i<D->Nx;i++) fprintf(xfile[Isc],"% le ",D->x[i]);
               fprintf(xfile[Isc],"\n");
               if (SC[Isc].FlexActive) {
                  for(i=0;i<D->Nf;i++) fprintf(uffile[Isc],"% le ",D->uf[i]);
                  fprintf(uffile[Isc],"\n");
                  for(i=0;i<D->Nf;i++) fprintf(xffile[Isc],"% le ",D->xf[i]);
                  fprintf(xffile[Isc],"\n");
               }
               if (SC[Isc].ConstraintsRequested) {
                  for(i=0;i<D->Nc;i++)
                     fprintf(ConstraintFile[Isc],"% le ",
                             D->GenConstraintFrc[i]);
                  fprintf(ConstraintFile[Isc],"\n");
               }
            }
         }
         if (SC[0].Exists) {
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
//////////////////////////////// added /////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
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
            
/////////////////////////////////////// added /////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //MagReport();
            //GyroReport();

////////////////////////////////////////////////////////// added ///////////////////////////////////////////////////////
            CmdReport();

            TxNum1 = FindBistaticGeo(0, 0, 6, 9);
            TxNum2 = FindBistaticGeo(0, 1, 11, 51);

  //          fprintf(ClosestTxFile, "%ld %ld\n",
  //          		TxNum1, TxNum2);
            /*
            for(ISp=0; ISp<Rx[0].SoOp[0].NSp;ISp){
                if(TxNum1[ISp] == 0) {
                    fprintf(PosNT1file,"%le %le %le\n",
                    		zerovec[0], zerovec[1], zerovec[2]);
                    MxV(World[EARTH].CWN,zerovec,PosWT1);
                }
                else {
                	SCT1 = &SC[TxNum1[ISp]];
                    fprintf(PosNT1file,"%le %le %le\n",
                    	SCT1->PosN[0],SCT1->PosN[1],SCT1->PosN[2]);
                    MxV(World[EARTH].CWN,SCT1->PosN,PosWT1);
                }
            }

            for(ISp=0; ISp<Rx[0].SoOp[1].NSp;ISp){
                if(TxNum2[ISp] == 0) {
                    fprintf(PosNT2file,"%le %le %le\n",
                    		zerovec[0], zerovec[1], zerovec[2]);
                    MxV(World[EARTH].CWN,zerovec,PosWT2);
                }
                else {
                    SCT2 = &SC[TxNum2[ISp]];
                    fprintf(PosNT2file,"%le %le %le\n",
                    	SCT2->PosN[0],SCT2->PosN[1],SCT2->PosN[2]);
                    MxV(World[EARTH].CWN,SCT2->PosN,PosWT2);
                }
            }*/


            if(TxNum1[0] == 0) {
                fprintf(PosNT1file,"%le %le %le\n",
                		zerovec[0], zerovec[1], zerovec[2]);
                MxV(World[EARTH].CWN,zerovec,PosWT1);
            }
            else {
                SCT1 = &SC[TxNum1[0]];
                fprintf(PosNT1file,"%le %le %le\n",
                	SCT1->PosN[0],SCT1->PosN[1],SCT1->PosN[2]);
                MxV(World[EARTH].CWN,SCT1->PosN,PosWT1);
            }

            if(TxNum2[0] == 0) {
                fprintf(PosNT2file,"%le %le %le\n",
                		zerovec[0], zerovec[1], zerovec[2]);
                MxV(World[EARTH].CWN,zerovec,PosWT2);
            }
            else {
                SCT2 = &SC[TxNum2[0]];
                fprintf(PosNT2file,"%le %le %le\n",
                	SCT2->PosN[0],SCT2->PosN[1],SCT2->PosN[2]);
                MxV(World[EARTH].CWN,SCT2->PosN,PosWT2);
            }

            fprintf(PosWT1file,"%18.12le %18.12le %18.12le\n",
               PosWT1[0],PosWT1[1],PosWT1[2]);
            fprintf(PosWT2file,"%18.12le %18.12le %18.12le\n",
               PosWT2[0],PosWT2[1],PosWT2[2]);

            WGS84Report(PosWT1, PosWT2, PosW);

            Sp1 = &Rx[0].SoOp[0].Sp[0];
            Sp2 = &Rx[0].SoOp[1].Sp[0];

            fprintf(SP1File, "%f %f %f %f %f %f %lu\n",
            		Sp1->PosN[0], Sp1->PosN[1], Sp1->PosN[2], Sp1->lat, Sp1->lon, Sp1->alt, Sp1->cnt);
            fprintf(SP2File, "%f %f %f %f %f %f %lu\n",
            		Sp2->PosN[0], Sp2->PosN[1], Sp2->PosN[2], Sp2->lat, Sp2->lon, Sp2->alt, Sp2->cnt);

            FindSwathArea(0);

            fprintf(SwathFile, "%f %f\n",
            		Rx[0].SwathArea, Rx[0].SwathWidth);
            //FindPathLength(0,0,0,TxNum1[0]);
            fprintf(pathLenFile1, "%f %f\n", Sp1->PathLenR);
            //FindPathLength(1,0,0,TxNum2[0]);
            fprintf(pathLenFile2, "%f\n", Sp2->PathLenR);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

         }

      }

      if (!strcmp(InOutPath,"./Potato/")) PotatoReport();
      

      if (CleanUpFlag) {
         fclose(timefile);
      }

}

/* #ifdef __cplusplus
** }
** #endif
*/
