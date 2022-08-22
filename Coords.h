#ifndef COORDS_H
#define COORDS_H

#include <math.h>
//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>      /* Needed only for _O_RDWR definition */
//#include <io.h>
//#include <windows.h>
#include <time.h>

#define PI  (3.141592653589793238462643)
#define r1 (PI/180)
#define mVC (149.896229) // [м/мкс]
#define Ff  (1.0/298.257) //Сжатие Земли
#define Ee  (2.0*Ff-Ff*Ff) //Эксцентриситет Земли в квадрате
#define Ua (7.0) //Угол между меридианом и большой осью антенны
#define Ub (10.0) //Угол между нормалью раскрыва рупора и зенитом
#define ca (cos(Ua*r1))
#define sa (sin(Ua*r1))
#define cb (cos(Ub*r1))
#define sb (sin(Ub*r1))
#define HR (246.0) //Длинна рупора в метрах
#define R0  (6378140.0) // Большая полуось Земли в метрах


/*Расчет артангенса угла*/
double ugol(double, double);
void L360(double *A);
double MJD(int y, int mes, int d);
void S0TIM(double MJDD, double *STG0, double *DTETA, double *RKPRIM);
//Расчет нутации
void NUT(double MJS, double *DKSI, double *DEPS, double *EPS);
void TimeTo(double Dd);
void FoundLineCtruct(double d,double ep1,double ep2,double fq1,double fq2,double *z1,double *z2);
void CentrDNR(double freq,double *ep0);
void CentrDNRPoly(double freq,double *ep0);
void EkvatorCoordFreqRadara(int year,int mes,int day,double t,double freq,double RA[],double DEC[]);

#endif
