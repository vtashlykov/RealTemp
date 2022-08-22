#include "RadarData.h"

bool operator< ( RadarChannelInfo a, RadarChannelInfo b )
{
	std::pair <unsigned, unsigned> aFreqAndLength, bFreqAndLength;
	aFreqAndLength=std::make_pair(a.PulseFreq, a.PulseLength);
	bFreqAndLength=std::make_pair(b.PulseFreq, b.PulseLength);
	return std::make_pair(a.Channel, aFreqAndLength) < std::make_pair(b.Channel, bFreqAndLength) ;
}

RadarData::RadarData(string FileName, string Direction){
	DataDirection=Direction;
	if(DataDirection=="in"){
		inp.open(FileName.c_str(), ios::in|ios::binary);
		ValidFile=inp.is_open();
		unsigned int i=0, CurrentAddress;
		RadarRecordHeader Head;
		//cout<<"Head Size : "<<sizeof(Head)<<" byte."<<endl;
	    while(Head.Char.Name[0]!='O'&&Head.Char.Name[1]!='R'&&Head.Char.Name[2]!='D'&&Head.Char.Name[3]!='A'&&Head.Int.Type!=3){
			inp.seekg(i);
			inp.read(Head.Total,9);
			if(Head.Char.Name[0]!='O'&&Head.Char.Name[1]!='R'&&Head.Char.Name[2]!='D'&&Head.Char.Name[3]!='A'&&Head.Int.Type!=3){
				cout<<i<<" "<<inp.tellg()<<" Name:";for(int ii=0;ii<4;ii++)cout<<Head.Char.Name[ii];cout<<". Type:"<<Head.Int.Type<<". Length:"<<Head.Int.Length<<"."<<endl;
			}
			i++;
	    }
		RadarVariable Var;
		//cout<<"Var Size : "<<sizeof(Var)<<" byte."<<endl;
		i=0;
		while(i<Head.Int.Length){
			inp.read(Var.Total,4);
			switch(Var.Short.Name){
				case 3:
					NumberAll=Var.Short.Value;
				break;
				case 5:
					FirstDelay=Var.Short.Value;
				break;
				case 30:
					SampleFreq=Var.Short.Value;
				break;
				case 33:
					Offset=Var.Short.Value;
				break;
			}
			i+=4;
		}
		short Channel, LengthSt1Long, LengthSt1Short, LengthSt2Long, LengthSt2Short;
		double Time, NewTime=-1;
		RadarFreq FreqSt1Long, FreqSt1Short, FreqSt2Long, FreqSt2Short;
		RadarChannelInfo ChInfo;
		map < int, RadarChannelInfo > RealizationInfo;
		int j=0;
		while(!inp.eof()){
		    j++;
			CurrentAddress=inp.tellg();
			inp.read(Head.Total,9);
			if(!inp.eof()){
				switch(Head.Int.Type){
					case 1:
						i=0;
						while(i<Head.Int.Length){
							inp.read(Var.Total,4);
							switch(Var.Short.Name){
								case 7:
									Channel=Var.Short.Value;
								break;
								case 9:
									ChInfo.Time.Year=Var.Short.Value;
								break;
								case 10:
									ChInfo.Time.Month=Var.Short.Value>>8;
									ChInfo.Time.Day=Var.Short.Value&0x00FF;
								break;
								case 11:
									ChInfo.Time.Min=Var.Short.Value>>8;
									ChInfo.Time.Hour=Var.Short.Value&0x00FF;
								break;
								case 12:
									ChInfo.Time.Sec=Var.Short.Value;
								break;
								case 13:
									ChInfo.Time.mSec=Var.Short.Value;
								break;
								case 14:
									FreqSt1Long.Short.Lo=Var.Short.Value;
								break;
								case 15:
									FreqSt1Long.Short.Hi=Var.Short.Value;
								break;
								case 16:
									FreqSt1Short.Short.Lo=Var.Short.Value;
								break;
								case 17:
									FreqSt1Short.Short.Hi=Var.Short.Value;
								break;
								case 18:
									FreqSt2Long.Short.Lo=Var.Short.Value;
								break;
								case 19:
									FreqSt2Long.Short.Hi=Var.Short.Value;
								break;
								case 20:
									FreqSt2Short.Short.Lo=Var.Short.Value;
								break;
								case 21:
									FreqSt2Short.Short.Hi=Var.Short.Value;
								break;
								case 22:
									LengthSt1Long=Var.Short.Value;
								break;
								case 23:
									LengthSt1Short=Var.Short.Value;
								break;
								case 24:
									LengthSt2Long=Var.Short.Value;
								break;
								case 25:
									LengthSt2Short=Var.Short.Value;
								break;
							}
							i+=4;
						}
						switch(Channel){
							case 0:
								ChInfo.PulseLength=LengthSt1Long;
								ChInfo.PulseFreq=FreqSt1Long.Total;
								ChInfo.Channel=0;
							break;
							case 1:
								ChInfo.PulseLength=LengthSt1Short;
								ChInfo.PulseFreq=FreqSt1Short.Total;
								ChInfo.Channel=1;
							break;
							case 2:
								ChInfo.PulseLength=LengthSt1Long;
								ChInfo.PulseFreq=FreqSt1Long.Total;
								ChInfo.Channel=2;
							break;
							case 3:
								ChInfo.PulseLength=LengthSt1Short;
								ChInfo.PulseFreq=FreqSt1Short.Total;
								ChInfo.Channel=3;
							break;
						}
						inp.read(Head.Total,9);
						//cout<<Head.Int.Type<<" "<<Head.Int.Length<<endl;
						ChInfo.LengthOfData=Head.Int.Length;
						//cout<<inp.tellg()<<" "<<j<<endl;
						ChInfo.AddressOfData=inp.tellg();
						Time=ChInfo.Time.GetYear();
						Sample_number=ChInfo.LengthOfData/4;
						Decimation=NumberAll/Sample_number;
						ChInfo.DecimationFreq=SampleFreq/Decimation;
						ChInfo.TotalDelay=FirstDelay-Offset-OFFSET_CONST;
						if(NewTime==-1) NewTime=Time;
						if(Time==NewTime){
							RealizationInfo[Channel]=ChInfo;
							/*cout<<" "<<Channel<<" ";
							cout<<RealizationInfo[Channel].Time.Year<<" ";
							cout<<RealizationInfo[Channel].Time.Month<<" ";
							cout<<RealizationInfo[Channel].Time.Day<<" ";
							cout<<RealizationInfo[Channel].Time.Hour<<" ";
							cout<<RealizationInfo[Channel].Time.Min<<" ";
							cout<<RealizationInfo[Channel].Time.Sec<<" ";
							cout<<RealizationInfo[Channel].Time.mSec<<endl;*/
						}
						else{
							Realizations[NewTime]=RealizationInfo;
							NewTime=Time;
							RealizationInfo.clear();
							RealizationInfo[Channel]=ChInfo;
							/*cout<<"  "<<Channel<<" ";
							cout<<RealizationInfo[Channel].Time.Year<<" ";
							cout<<RealizationInfo[Channel].Time.Month<<" ";
							cout<<RealizationInfo[Channel].Time.Day<<" ";
							cout<<RealizationInfo[Channel].Time.Hour<<" ";
							cout<<RealizationInfo[Channel].Time.Min<<" ";
							cout<<RealizationInfo[Channel].Time.Sec<<" ";
							cout<<RealizationInfo[Channel].Time.mSec<<endl;*/
						}
						//cout<<setprecision(12)<<Time<<" "<<NewTime<<setprecision(3)<<" "<<Realizations.size()<<endl;
						inp.seekg(ChInfo.AddressOfData+Head.Int.Length);
					break;
					case 2:
						inp.seekg(CurrentAddress+9+Head.Int.Length);
					break;
					case 3:
						i=0;
						while(i<Head.Int.Length){
							inp.read(Var.Total,4);
							switch(Var.Short.Name){
								case 3:
									NumberAll=Var.Short.Value;
								break;
								case 5:
									FirstDelay=Var.Short.Value;
								break;
								case 30:
									SampleFreq=Var.Short.Value;
								break;
							}
							i+=4;
						}
					break;
				}
			}
		}
		Realizations[NewTime]=RealizationInfo;
		//cout<<setprecision(12)<<Time<<" "<<NewTime<<setprecision(3)<<" "<<Realizations.size()<<endl;
		/*for(auto it=Realizations.begin(); it!=Realizations.end();it++){
			cout<<setprecision(12)<<it->first<<setprecision(3)<<" "<<it->second[3].Time.Day<<endl;
		}*/
		inp.clear();
		inp.seekg(0);
	}
	if(DataDirection=="out"){
		out.open(FileName.c_str(),ios::out|ios::binary);
	}
}

RadarData::~RadarData(){
	if(DataDirection=="in"){inp.close();}
	if(DataDirection=="out"){out.close();}
	this->Realizations.clear();
}

Containers::Containers()
{

}

Containers::Containers(vector<RadarChannelInfo> Info, size_t nheight)
{
	for(unsigned i=0; i<Info.size(); i++)
	{
		this->Interpol[Info[i]].resize(nheight);
		cout<<"________Check container "<<this<<" : "<<Info[i].PulseFreq<<"\t"<<Info[i].PulseLength<<"\t"\
			<<this->Interpol[Info[i]].size()<<"\n";
	}
}

Containers::Containers(const Containers &other)
{
	map<RadarChannelInfo, vector<vector<double> > > Copy;
	Copy=other.Interpol;
	this->Interpol=Copy;
}

Containers::~Containers()
{
	this->Interpol.clear();
}

void Containers::getC(vector<double> &P, RadarChannelInfo info, size_t height)
{
	P=this->Interpol[info][height];
	// cout<<"________Check container "<<this<<" : "<<this->Interpol[info].size()<<"\t"\
	// 	<<this->Interpol[info][height].size()<<"\t"<<P.size()<<"\n";
}

void Containers::setC(vector<double> P, RadarChannelInfo info, size_t height)
{
	// this->Interpol[info][height]=P;
	// cout<<"________Check container "<<this<<" : "<<this->Interpol[info].size()<<"\n";//\
		<<this->Interpol[info][height].size()<<"\t"<<P.size()<<"\n";

	// for(int i=0; i<this->Interpol[info][height].size(); i++)
	// 	cout<<Interpol[info][height][i]<<"\t";
	// cout<<endl;
}

void Containers::checkC(RadarChannelInfo info, size_t height)
{
	cout<<"________Check container "<<this<<" : "<<this->Interpol[info].size()<<"\n";//\
		<<this->Interpol[info][height].size()<<"\n";
}


void RadarData::ReadQuadratures(map < double, map < int, RadarChannelInfo > >::iterator Global, map< int, \
	vector < complex < double > > > &Quadratures, vector< double > &Powers)
{
	Quadratures.clear();
	complex < double > ComplexPoint;
	for(auto Local=Global->second.begin(); Local != Global->second.end(); Local++){
		inp.seekg(Local->second.AddressOfData);
		char Data[Local->second.LengthOfData];
		inp.read(Data,Local->second.LengthOfData);
		double Pwr=0;
		unsigned int j=0;
		for(unsigned int i=0; i<Local->second.LengthOfData/2; i+=2){
			// ComplexPoint.real(-1.0*double(*(short*)(Data+i)));//clarify -1.0
			ComplexPoint.real(double(*(short*)(Data+i)));
			ComplexPoint.imag(double(*(short*)(Data+i+Local->second.LengthOfData/2)));
			Quadratures[Local->first].push_back(ComplexPoint);
			Pwr+=norm(ComplexPoint);
			j++;
		}
		Powers.push_back(Pwr);
	}
}

/* -----------------------------------------Range to height---------------------------------------*/

double AzimMaxDNR(double freqkHz)
{//Возвращаем угол в радианах
	double ep0;
	ep0=185634.9983-4.877867876*freqkHz+4.803667396E-005*freqkHz*freqkHz-2.102600271E-010*freqkHz*freqkHz*freqkHz+3.453540782E-016*freqkHz*freqkHz*freqkHz*freqkHz;
	ep0+=0.2;
	ep0*=Rad;
	return ep0;
}

void Found_H_Fi(double z, double R, double &H, double &Fi)
{
    int j;
    double C, fi, x;
    fi=atan(z/R);
    for(j=1; j<=15; j++)
	{
        C=1.0/sqrt(1.0-Ee*sin(fi)*sin(fi));
        x=z+R0*C*Ee*sin(fi);
        fi=atan(x/R);
    }
    C=1.0/sqrt(1.0-Ee*sin(fi)*sin(fi));
    x=R/cos(fi)-R0*C;

    Fi=fi;
    H=x;
}

void Found_FI_LG_H__for_Ep_Gam_R(double ep, double gam, double R, double FI, double LG, double HN, double &fi, double &lg, double &H)
{
	double RX, RY, RZ, ZN, RR, delt, th;
	double el, az;

	ZN=R0/sqrt(1.0-Ee*sin(FI)*sin(FI));
	RX=(ZN+HN)*cos(FI)*cos(LG);
	RY=(ZN+HN)*cos(FI)*sin(LG);
	RZ=(ZN*(1.0-Ee)+HN)*sin(FI);

	el=asin(cos(Ub*Rad-gam)*cos(ep));
	az=1.5*PI-Ua*Rad-atan2(sin(ep),cos(ep)*sin(Ub*Rad-gam));
	delt=asin(sin(el)*sin(FI)+cos(az)*cos(el)*cos(FI));
	th=-atan2(-sin(az)*cos(el),sin(el)*cos(FI)-cos(az)*cos(el)*sin(FI));

	RX+=R*cos(delt)*cos(LG+th);
	RY+=R*cos(delt)*sin(LG+th);
	RZ+=R*sin(delt);
	RR=sqrt(RX*RX+RY*RY);
	lg=atan2(RY, RX);
	Found_H_Fi(RZ, RR, H, fi);
}

void Topoc_to_Ant(double el, double az, double *gam, double *ep)
{
	double Az, El, gam0;

	Az=az+Ua*Rad;
	El=el;
	*ep=asin(-cos(Az)*cos(El));
	gam0=Ub*Rad-atan2(-sin(Az)*cos(El),sin(El));
//	*gam=asin(sin(gam0)*cos(*ep));
}

void Ant_to_Topoc(double ep, double gam, double *el, double *az)
{
	double gam0;

	gam0=asin(sin(gam)/cos(ep));
	*el=asin(cos(Ub*Rad-gam0)*cos(ep));
	*az=1.5*PI-Ua*Rad-atan2(sin(ep), cos(ep)*sin(Ub*Rad-gam0));
}


/*void RadarData::ChangeFreq(void){
	map < double, map < int, RadarChannelInfo > >::iterator Global=Realizations.begin();
	int previous=Global->second[0].PulseFreq;
	while(Global != Realizations.end()){
		if(Global->second[0].PulseFreq-previous==300){
			Global->second[1].PulseFreq+=1;
			Global->second[3].PulseFreq+=1;
		}
		previous=Global->second[0].PulseFreq;
		Global++;
	}
}

void RadarData::ConvoluteChannel(map < double, map < int, RadarChannelInfo > >::iterator Global, string &Sequence, int Channel){
	double PSKBinDuration=double(Global->second[Channel].PulseLength)/double(Sequence.size());
	double DataTimeStep=NumberAll/(Global->second[Channel].LengthOfData/4)*1000/SampleFreq;
	int LocalTime=0;
	int PulsePhaseLength=int(Global->second[Channel].PulseLength/DataTimeStep);
	double PulsePhase[PulsePhaseLength];
	string::iterator SequenseIt=Sequence.begin();
	int Sign;
	if (*SequenseIt == '0') Sign=1;
	else Sign=-1;
	SequenseIt++;
	for(int i=0; i < PulsePhaseLength; i++){
		if(LocalTime>=PSKBinDuration){
			if(*SequenseIt == '0') Sign=1;
			else Sign=-1;
			SequenseIt++;
			LocalTime=LocalTime-PSKBinDuration;
		}
		PulsePhase[i]=Sign;
		LocalTime+=4;
	}
	vector < complex < double > >::iterator ScanIt=QData[Channel].begin();
	while(ScanIt!=QData[Channel].end()){
		complex < double > Product;
		vector <  complex < double > >::iterator ConvIt=ScanIt;
		for(int i=0; i < PulsePhaseLength && ConvIt!=QData[Channel].end(); i++){
			Product+=*ConvIt*PulsePhase[i];
			ConvIt++;
		}
		*ScanItProduct;
		ScanIt++;
	}
}

double RadarData::MultiplayChannels(vector < complex < double > > &MultiplicationResult, int FirstChannel, int SecondChannel){
	double CorrelationCoefficient=0;
	double FirstPwr=0, SecondPwr=0;
	MultiplicationResult.clear();
	vector < complex < double > >::iterator FirstIt=QData[FirstChannel].begin();
	vector < complex < double > >::iterator SecondIt=QData[SecondChannel].begin();
	while(FirstIt!=QData[FirstChannel].end()||SecondIt!=QData[SecondChannel].end()){
		FirstPwr+=norm(*FirstIt);
		SecondPwr+=norm(*SecondIt);
		complex < double > Mult=*FirstIt*conj(*SecondIt);
		MultiplicationResult.push_back(Mult);
		CorrelationCoefficient+=norm(Mult);
		FirstIt++;
		SecondIt++;
	}
	return CorrelationCoefficient/sqrt(FirstPwr*SecondPwr);
}
*/
