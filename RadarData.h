#ifndef RADARDATA_H
#define RADARDATA_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <complex>

#include "RadarTime.h"
#include "Coords.h"

#define OFFSET_CONST 1010

#define PI (3.141592653589793238462643)
#define Rad (PI/180)
#define R0 (6378140.0) // Большая полуось Земли в метрах
#define Ff (1.0/298.257) //Сжатие Земли
#define Ee (2.0*Ff-Ff*Ff) //Эксцентриситет Земли в квадрате
#define Ua (7.0) //Угол между меридианом и большой осью антенны
#define Ub (10.0) //Угол между нормалью раскрыва рупора и зенитом
#define nr_lg 103.257
#define nr_lat 52.8814     //радар
#define ir_lg 104.27
#define ir_lat 52.25      //Иркутск
#define tr_lg 103.07611
#define tr_lat 51.8092  //Торы
#define graph_koef 75
#define lat_off 49
#define lng_off 97

using namespace std;

#pragma pack(push,1)
union RadarRecordHeader{ //#pragma pack(push/pop)  for structure alignment
	char Total[9];
	struct{
		char Name[4];
		char Type;
		char Length[4];
	}Char;
	struct{
		unsigned Name:32;
		unsigned Type:8;
		unsigned Length:32;
	}Int;
};
union RadarVariable{
	char Total[4];
	struct{
		short Name;
		short Value;
	}Short;
};
union RadarFreq{
	int Total;
	struct{
		short Lo;
		short Hi;
	}Short;
};
#pragma pack(pop)

struct RadarChannelInfo{
	unsigned Channel;
	RadarTime Time;
	unsigned PulseLength;
	unsigned PulseFreq;
	unsigned AddressOfData;
	unsigned LengthOfData;
	unsigned DecimationFreq;
	unsigned TotalDelay;
};

bool operator< ( RadarChannelInfo a, RadarChannelInfo b );

class RadarData{

	private:
		string		DataDirection;
		ifstream	inp;
		ofstream	out;

	public:

		bool	ValidFile;
		short	FirstDelay;
		short	NumberAll;
		short	SampleFreq;
		short	Offset;
		unsigned Sample_number;
		unsigned Decimation;

		map < double, map < int, RadarChannelInfo > > Realizations;

		RadarData(string, string);
		~RadarData();

		void ReadQuadratures(map < double, map < int, RadarChannelInfo > >::iterator, map< int, vector < complex < double > > >&, vector< double >&);
};

class Containers{
	private:
		map<RadarChannelInfo, vector<vector<double> > > Interpol;
	public:
		void getC(vector<double>&, RadarChannelInfo, size_t);
		void setC(vector<double>, RadarChannelInfo, size_t);
		void checkC(RadarChannelInfo, size_t);

		Containers();
		Containers(vector<RadarChannelInfo>,  size_t);
		Containers(const Containers &other);
		~Containers();
};

/* -----------------------------------------Range to height---------------------------------------*/

double AzimMaxDNR(double freqkHz);
void Found_H_Fi(double z, double R, double &H, double &Fi);
void Found_FI_LG_H__for_Ep_Gam_R(double ep, double gam, double R, double FI, double LG, double HN,\
	 double &fi, double &lg, double &H);
void Topoc_to_Ant(double el, double az, double *gam, double *ep);
void Ant_to_Topoc(double ep, double gam, double *el, double *az);

#endif
