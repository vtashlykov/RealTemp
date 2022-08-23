#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <stdlib.h>
#include <complex>
#include <vector>
#include <string>
#include <stdio.h>
#include <dirent.h>
#include <iomanip>
#include <map>
#include <omp.h>
#include <utility>
#include <thread>
#include <chrono>
#include "RadarData.h"
#include "InverseProblemsolver.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>

#define TLENGTH 1024
#define LIGHT_VELOCITY 0.3

using namespace std;

int sel_ISE(const dirent *d){
	string name=d->d_name;
	string subname;
	if ( name == "." || name == "..") return 0;
	size_t pos=name.find(".");
	if ( pos == string::npos ) return 0;
	subname=name.substr(pos);
	if(subname==".ISE") return 1;
	else return 0;
}

int sel_dat(const dirent *d){
	string name=d->d_name;
	string subname;
	if ( name == "." || name == "..") return 0;
	size_t pos=name.find(".");
	if ( pos == string::npos ) return 0;
	subname=name.substr(pos);
	if(subname==".dat") return 1;
	else return 0;
}

double LT_transform(double H){
	H+=7.0;
	if(H>=24.0) H-=24.0;
	return double(H);
}

void UT_transform(int &DoY, double &H){
	H-=7.0;
	if(H<0) 
	{
		H+=24.0;
		DoY-=1;
	}
}

template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element)
{
    std::pair<bool, int > result;
    // Find given element in vector
    auto it=std::find(vecOfElements.begin(), vecOfElements.end(), element);
    if (it != vecOfElements.end())
    {
        result.second=distance(vecOfElements.begin(), it);
        result.first=true;
    }
    else
    {
        result.first=false;
        result.second=-1;
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma warning(disable:4756)
int main(int argc, char* argv[])
{
	gsl_set_error_handler_off();
	double init_time=omp_get_wtime();
	
	if(argc!=2)
	{
		cout<<"usage: $PATH\\RealTemp ParFile"<<endl;
		exit(1);
	}

	int day_0, day_n;
	int Y, M, D, DoY;
	string First_Day, Last_Day, YM;
	vector<string> Date;
	string cmd, Path, Dir;

	ifstream inp;
	string line;
	size_t found=0;
	unsigned L_700=0, L_900=0;
	RadarChannelInfo info_0_700;
	
	inp.open(argv[1]);
	getline(inp, line);
	getline(inp, line);
	found=line.find_first_of("=");
	Y=atoi(line.substr(found+2, 4).c_str());
	printf("Y %d\n", Y);
	
	getline(inp, line);
	found=line.find_first_of("=");
	Dir=line.substr(found+2);
	// Dir.erase(Dir.end()-1);
	printf("Dir %s\n", Dir.c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	M=atoi(line.substr(found+2, 2).c_str());
	D=atoi(line.substr(found+4, 2).c_str());
	day_0=getDoy(Y, M, D);
	printf("First_Day %d\n", day_0);
	
	getline(inp, line);
	found=line.find_first_of("=");
	M=atoi(line.substr(found+2, 2).c_str());
	D=atoi(line.substr(found+4, 2).c_str());
	day_n=getDoy(Y, M, D);
	printf("Last_Day %d\n", day_n);

	for(int d=day_0; d<=day_n; d++)
	{
		string s;
		getDate(Y, d, M, D);
		s=to_string(Y);
		if(M<10)
			s=s+"0"+to_string(M);
		else
			s=s+to_string(M);
		if(D<10)
			s=s+"0"+to_string(D);
		else
			s=s+to_string(D);
		Date.push_back(s);
	}
	
	getline(inp, line);
	found=line.find_first_of("=");
	const unsigned Integration_Time=atoi(line.substr(found+2, 3).c_str());
	printf("Integration_Time %s\n", line.substr(found+2, 3).c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	const unsigned Number_of_samples=atoi(line.substr(found+2).c_str());
	printf("Number_of_samples %s\n", line.substr(found+2).c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	const unsigned Ground_Clutter_End=atoi(line.substr(found+2, 3).c_str());
	printf("Ground_Clutter_End %s\n", line.substr(found+2, 3).c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	const unsigned Ground_Clutter_Start=atoi(line.substr(found+2, 2).c_str());
	printf("Ground_Clutter_Start %s\n", line.substr(found+2, 2).c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	const unsigned Lag_Truncate=atoi(line.substr(found+2, 2).c_str());
	printf("Lag_Truncate %s\n", line.substr(found+2, 2).c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	unsigned Number_of_threads=atoi(line.substr(found+2).c_str());
	printf("Number_of_threads %s\n", line.substr(found+2).c_str());

	getline(inp, line);
	found=line.find_first_of("=");
	string Path_to_result_files=line.substr(found+2).c_str();
	printf("Path_to_result_files %s\n", line.substr(found+2).c_str());
	
	getline(inp, line);
	found=line.find_first_of("=");
	unsigned UT0=atoi(line.substr(found+2).c_str());
	printf("UT0 %d\n", UT0);
	
	getline(inp, line);
	found=line.find_first_of("=");
	unsigned NUMBER_OF_INTERVALS=atoi(line.substr(found+2).c_str());
	printf("NUMBER_OF_INTERVALS %d\n", NUMBER_OF_INTERVALS);

	double Time_Step=6.0;

	inp.close();

	printf("%f!\n", Time_Step);

	for(size_t d=0; d<Date.size(); d++)
	{
		// ostringstream t;
		// t<<(d/10)-(d/100)*10;
		// t<<(d%10);
		// Date[d]=Y+M+string(t.str());

		//string Path=  "/mnt/data/orda/experiments/filtered/201301/20130115/";//52
		// string Path="/mnt/data/orda/experiments/"+Y+"/"+YM+"/"+Date+"/";//48
		// string Path="/mnt/data2/orda/experiments/4archive/2019/test/"+Date+"/";//67
		Path=Dir+"/"+Date[d]+"/";
		// Path=Dir+"/"+Y+"/"+YM+"/"+Date[d]+"/";
		// Path=Dir+"/"+Y+"/"+YM+"/filtered/"+Date[d]+"/";
		// Path=Dir+"/"+Y+"/filter2/"+Date+"/";
		// cout<<Date<<"\t"<<t.str()<<"\n";
		cout<<Path<<"\n";
		string out_Path=Path_to_result_files+Date[d].substr(0, 6)+"/";
		string LastFileName=Date[d]+"_0000_000_0000_002_000.ISE";
		string FinalFileName=Date[d]+"_2359_000_0000_002_000.ISE";
		dirent **entry;
		unsigned NumberOfFiles=0;
		NumberOfFiles=scandir(Path.c_str(), &entry, sel_ISE, alphasort);

		string dir;
		cmd="mkdir "+out_Path;
		system(cmd.c_str());

		dir=out_Path+"CorMatrices/"+Date[d].substr(0, 6)+"/";
		// int NumberOfDates=scandir(out_Path.c_str(), &outry, NULL, alphasort);
		// cmd="rmdir --ignore-fail-on-non-empty "+dir;
		cmd="mkdir "+dir;
		system(cmd.c_str());
		dir=out_Path+"CorMatrices/"+Date[d]+"/";
		// cmd="rmdir --ignore-fail-on-non-empty "+dir;
		// system(cmd.c_str());
		cmd="mkdir "+dir;
		system(cmd.c_str());

		// cmd="rm /home/tashlykov/test/ConfigFile_*";
		// system(cmd.c_str());
		// cmd="rm /home/tashlykov/test/List_*";
		// system(cmd.c_str());

		cout<<"Total number of files is "<<NumberOfFiles<<"."<<endl;
		unsigned start=0;
		cout<<"File # "<<start<<" : "<<entry[start]->d_name<<endl;
		unsigned stop=0;
		while(entry[stop]->d_name!=FinalFileName && stop<NumberOfFiles-1) stop++;
		cout<<"File # "<<stop<<" "<<entry[stop]->d_name<<endl;
		RadarTime Time;

		vector<map<RadarChannelInfo, unsigned> > Nsample, Nsample2;
		vector<map<RadarChannelInfo, double> > Te_mean, Tr_mean, Norm, Aconst, Bconst;
		vector<map<RadarChannelInfo, vec> > Hist, Chi, Chi1, UT_hist, RawPower, Noise, mHist, dHist;
		vector<map<RadarChannelInfo, vector<vec> > > Err, Err1, InterPol, InterPol1;
		vector<map<RadarChannelInfo, vec> > MHist, DHist, height_s, P_model, Ne, P_exp, P_conv;
		vector<map<RadarChannelInfo, vecc> > Kr;
		vector<vector<RadarChannelInfo> > Info;
		vector<map<RadarChannelInfo, bool> > valid;
		vector<map<RadarChannelInfo, vector<vecc> > > Matrix, R1;
		valid.resize(NUMBER_OF_INTERVALS);
		Nsample.resize(NUMBER_OF_INTERVALS);
		Te_mean.resize(NUMBER_OF_INTERVALS);
		Tr_mean.resize(NUMBER_OF_INTERVALS);
		Norm.resize(NUMBER_OF_INTERVALS);
		Nsample2.resize(NUMBER_OF_INTERVALS);
		Hist.resize(NUMBER_OF_INTERVALS);
		Err.resize(NUMBER_OF_INTERVALS);
		Err1.resize(NUMBER_OF_INTERVALS);
		InterPol.resize(NUMBER_OF_INTERVALS);
		InterPol1.resize(NUMBER_OF_INTERVALS);
		Chi.resize(NUMBER_OF_INTERVALS);
		Chi1.resize(NUMBER_OF_INTERVALS);
		UT_hist.resize(NUMBER_OF_INTERVALS);
		RawPower.resize(NUMBER_OF_INTERVALS);
		Noise.resize(NUMBER_OF_INTERVALS);
		Info.resize(NUMBER_OF_INTERVALS);
		MHist.resize(NUMBER_OF_INTERVALS);
		DHist.resize(NUMBER_OF_INTERVALS);
		mHist.resize(NUMBER_OF_INTERVALS);
		dHist.resize(NUMBER_OF_INTERVALS);
		Kr.resize(NUMBER_OF_INTERVALS);
		height_s.resize(NUMBER_OF_INTERVALS);
		P_model.resize(NUMBER_OF_INTERVALS);
		Ne.resize(NUMBER_OF_INTERVALS);
		P_exp.resize(NUMBER_OF_INTERVALS);
		P_conv.resize(NUMBER_OF_INTERVALS);
		Aconst.resize(NUMBER_OF_INTERVALS);
		Bconst.resize(NUMBER_OF_INTERVALS);
		Matrix.resize(NUMBER_OF_INTERVALS);
		R1.resize(NUMBER_OF_INTERVALS);

		unsigned CurrentNumberOfFile;
		vector<int> ThreadIsBuisy;
		map<int, vector<string> > ListMap;
		ThreadIsBuisy.assign(Number_of_threads, 0);

		omp_set_num_threads(Number_of_threads);
		printf("Parallel calculations start.\n");
		#pragma omp parallel
		{
			CurrentNumberOfFile=0;
			ofstream out, out_info;
			double iter_start=omp_get_wtime(), iter_end, phi;
			// map<RadarChannelInfo, vecc > Phi, Amp;
			// map<RadarChannelInfo, vec> Chisq, Chisqa;
			string FileName, List;

			vecc Q, Q0, Qu;
			// map<RadarChannelInfo, vecc> Qprev;
			vector<double> ChannelsPowers;
			map<int, vecc> TimeQuadratures;
			double UT;
			RadarChannelInfo info;
			vector<unsigned> count;
			
			// schedule(static, NUMBER_OF_INTERVALS/Number_of_threads)
			
			#pragma omp for schedule(dynamic)
			for(size_t it=0; it<NUMBER_OF_INTERVALS; it++)
			{
				int Time_it=(it+UT0)*10;
				unsigned Pid=omp_get_thread_num();
				string list;
				count.resize(2);
				for(size_t k=start; k<=stop; k++)
				{
					FileName=Path+entry[k]->d_name;
					size_t found=FileName.find_first_of("_");
					int minute=atoi(FileName.substr(found+3, 2).c_str());
					int hour=atoi(FileName.substr(found+1, 2).c_str());
					// printf("Quadratures analysing for %s\t%d:%d by thread #%d\n", \
						// FileName.c_str(), hour, minute, Pid);
					if(count[0]==Number_of_samples)
						count[0]=0;
					if(count[1]==Number_of_samples)
						count[1]=0;
					if(count[0]==2000)
						count[0]=0;
					if(count[1]==2000)
						count[1]=0;
					if((hour*60+minute>=Time_it-int(Integration_Time) and \
						hour*60+minute<Time_it+int(Integration_Time)))
					{
						printf("Quadratures analysing for %s\t%d:%d - %d:%d by thread #%d\n", \
							FileName.c_str(), hour, minute, Time_it/60, Time_it%60, Pid);
						RadarData *MyData=new RadarData(FileName.c_str(), "in");
						if(!MyData->ValidFile)
							printf("File not found.\n");
						for(auto N=MyData->Realizations.begin(); N!=MyData->Realizations.end(); N++)
						{
							for(int ch=0; ch<2; ch++)
							{
								info=N->second[2*ch];
								MyData->ReadQuadratures(N, TimeQuadratures, ChannelsPowers);
								Time.PutYear(N->first);
								UT=Time.GetHour();
								
								if(!isnormal(info.PulseLength) or !isnormal(info.PulseFreq) or !isnormal(info.LengthOfData))
									fprintf(stderr, "Wrong signal info!	PulseLength=%d, PulseFreq=%d,\
										DataLength=%d, info.DecimationFreq=%d\n", info.PulseLength, \
										info.PulseFreq, info.LengthOfData, info.DecimationFreq);

								double deltat=1000.0/double(info.DecimationFreq);
								unsigned L=floor(double(info.PulseLength)/(1000.0/double(info.DecimationFreq)));
								unsigned Tlength=TimeQuadratures[2*ch].size();

								if(Nsample[it].count(info)==0)
								{
									Nsample[it][info]=0;
									Matrix[it][info].resize(Tlength-L, vecc (L));
									R1[it][info].resize(Tlength-L, vecc (L));
									Noise[it][info].resize(L);
									Te_mean[it][info]=0.0;
									Tr_mean[it][info]=0.0;
									Norm[it][info]=0.0;
									vec Zero;
									MHist[it][info].resize(Number_of_samples);
									DHist[it][info].resize(Number_of_samples);
									UT_hist[it][info].resize(Number_of_samples);
									Kr[it][info].resize(Number_of_samples);
									Zero.resize(2000);
									mHist[it][info]=Zero;
									dHist[it][info]=Zero;
								}
									
								Q.resize(Tlength);
								Qu.resize(Tlength);
								if(Nsample[it][info]<Number_of_samples)
								{
									if(ch==0)
									{
										if(List.empty())
											List=FileName;
										else if(list!=FileName)
											List+="\n"+FileName;
										list=FileName;
										if(Nsample[it][info]==0)
										{
											ListMap[Pid].push_back(FileName);
										}
									}
									size_t j=0;
									for(auto t=TimeQuadratures[2*ch].begin(); \
										t!=TimeQuadratures[2*ch].end(); t++)
									{
										Q[j]=complex<double> (real(*t), imag(*t));
										j++;
									}

									if(ch==0)
										Qu=Q;
									double Qm=0.0, Qd=0.0;
									complex<double> num, den1, den2;
									num=complex<double> (0.0, 0.0);
									den1=complex<double> (0.0, 0.0);
									den2=complex<double> (0.0, 0.0);
									for(size_t t=0; t<Tlength; t++)
									{
										double r=double(t)*deltat*LIGHT_VELOCITY/2.0+\
											double(info.TotalDelay)*0.15;
										if(r>250.0 and r<550.0)
										{
											Qm+=real(Q[t]);
											Qd+=real(Q[t])*real(Q[t]);
											num+=Q[t]*conj(Qu[t]);
											den1+=norm(Q[t]);
											den2+=norm(Qu[t]);
										}
									}
									Qm/=double(Tlength-Ground_Clutter_End);
									Qd/=double(Tlength-Ground_Clutter_End);
									MHist[it][info][Nsample[it][info]]=Qm;
									DHist[it][info][Nsample[it][info]]=sqrt(Qd);
									Kr[it][info][Nsample[it][info]]=num/sqrt(den1*den2);
									if(abs(Qm)<20.0)
									{
										j=round(abs(Qm)*100.0);
										mHist[it][info][j]++;
									}
									if(sqrt(Qd)<20.0)
									{
										j=round(sqrt(Qd)*100.0);
										dHist[it][info][j]++;
									}

									Nsample[it][info]++;
									// Qt[info].push_back(Q);
									UT_hist[it][info][Nsample[it][info]]=UT;

									for(size_t t=0; t<Tlength-L; t++)
									{
										for(size_t tau=0; tau<L; tau++)
											Matrix[it][info][t][tau]+=Q[t]*conj(Q[t+tau]);
									}
								}
							}
						}
						delete MyData;
					}
				}
#pragma omp critical
				{	
					printf("\nThread %d List: ", Pid);
					for(unsigned i=0; i<Number_of_threads; i++)
					{
						for(size_t j=0; j<ListMap[Pid].size(); j++)
						{	
							if(i!=Pid and ListMap[Pid][j]==ListMap[i][j])
							{
								ThreadIsBuisy[i]++;
								printf("%s\t", ListMap[i][j].c_str());
							}
						}
					}
					printf("\nThreadIsBuisy %d!\n", ThreadIsBuisy[Pid]);
				}

				iter_end=omp_get_wtime();
				// #pragma omp flush(CurrentNumberOfFile)
				printf("Processing for %d:%d UT by thread #%d.\nProcessing time is %g seconds.\n%f %% completed.\n\n",
					Time_it/60, Time_it%60, Pid, iter_end-iter_start, 100.0*double(it)/double(NUMBER_OF_INTERVALS));

				string current_time;
				ostringstream tt;
				double hour=floor(Time_it/60), min=Time_it%60;
				tt<<floor(hour/10);
				tt<<hour-floor(hour/10)*10;
				tt<<floor(min/10);
				tt<<min-floor(min/10)*10;
				current_time=tt.str();

				for(auto& t: Matrix[it])
				{
					RadarChannelInfo info=t.first;
					if(isnormal(info.PulseLength) and isnormal(info.PulseFreq) \
						and isnormal(info.LengthOfData) and Nsample[it][info]>0)
							Info[it].push_back(info);
				}

				// for(size_t m=0; m<Info[it].size(); m++)
				// {
				// 	RadarChannelInfo info=Info[it][m];
				// 	double deltat=1000.0/double(info.DecimationFreq);
				// 	unsigned L=floor(double(info.PulseLength)/deltat);
				// 	unsigned Tlength=info.LengthOfData/4;
				// 	printf("______________Time %d; Channel %d; Frequency %d; Pulse length %d; DataLength %d; DeltaT %d, N %d for thread No %d\n",\
				// 		it, info.Channel, info.PulseFreq, info.PulseLength, info.LengthOfData, int(1000.0/double(info.DecimationFreq)), Nsample[it][info], omp_get_thread_num());
				// 	printf("Matrix size [%d, %d, %d] for thread No %d\n", \
				// 		Matrix.size(), Matrix[it][info].size(), Matrix[it][info][0].size(), omp_get_thread_num());
				// 	printf("%d_%d\n", L, Tlength);	
				// }

				for(size_t m=0; m<Info[it].size(); m++)
				{
					RadarChannelInfo info=Info[it][m];
					double deltat=1000.0/double(info.DecimationFreq);
					unsigned L=floor(double(info.PulseLength)/deltat);
					unsigned Tlength=info.LengthOfData/4;

					for(size_t t=0; t<Tlength-L; t++)
					{
						for(size_t tau=0; tau<L; tau++)
							Matrix[it][info][t][tau]/=double(Nsample[it][info]);
					}
					for(size_t tau=0; tau<L; tau++)
						Noise[it][info][tau]=real(Matrix[it][info][Tlength/2][tau]);

					for(size_t t=0; t<Tlength-L; t++)
					{
						for(size_t tau=0; tau<L; tau++)
						{
							for(size_t r=0; r<L-tau; r++)
							{
								if(t+r<Tlength-L)
									R1[it][info][t][tau]+=Matrix[it][info][t+r][tau];
							}
							R1[it][info][t][tau]/=double(L-tau);
						}
					}

					InterPol[it][info].resize(Tlength-L);
					InterPol1[it][info].resize(Tlength-L);
					RawPower[it][info].resize(Tlength-L);
					Err[it][info].resize(Tlength-L);
					Err1[it][info].resize(Tlength-L);
					Chi[it][info].resize(Tlength-L);
					Chi1[it][info].resize(Tlength-L);
					for(size_t t=0; t<Tlength-L; t++)
					{
						InterPol[it][info][t].resize(POLY_DEGREE);
						Err[it][info][t].resize(POLY_DEGREE);
						InterPol1[it][info][t].resize(POLY_DEGREE);
						Err1[it][info][t].resize(POLY_DEGREE);
					}
					// if(omp_get_thread_num()==0)
					// 	printf("L=%d\nMatrix size is [%d:%d]\n", L, CM[info].size(), CM[info][0].size());
					
					// #pragma omp_critical
					{
						printf("Interpolaton... by thread %d ", omp_get_thread_num());
						if(L>0)
						{
							const size_t n=L;
							// double y[n];
							vec y(n);
							y.assign(n, 0.0);
							if(n<POLY_DEGREE)
								printf("n<POLY_DEGREE; Time %d; Channel %d; Frequency %d; Pulse length %d; DataLength %d\n",\
									it, info.Channel, info.PulseFreq, info.PulseLength, info.LengthOfData);

							// #pragma omp_set_barrier							
							for(size_t t=0; t<Tlength-L; t++)
							{
								for(size_t tau=0; tau<n-2*Lag_Truncate; tau++)
									y[tau]=real(Matrix[it][info][t][tau+Lag_Truncate]);
								for(size_t tau=n-2*Lag_Truncate; tau<n; tau++)
									y[tau]=0.0;
								ReACF_interpolator(y, Lag_Truncate, InterPol[it][info][t], \
									Err[it][info][t], Chi[it][info][t]);
								RawPower[it][info][t]=real(Matrix[it][info][t][0]);

								// CM_print=true;

								for(size_t tau=0; tau<n-2*Lag_Truncate; tau++)
									y[tau]=real(R1[it][info][t][tau+Lag_Truncate]);
								for(size_t tau=n-2*Lag_Truncate; tau<n; tau++)
									y[tau]=0.0;
								ReACF_interpolator(y, Lag_Truncate, InterPol1[it][info][t], \
									Err1[it][info][t], Chi1[it][info][t]);
							}
						}
						else
						{
							printf("Wrong signal!\n\n");
						}	
						printf("done %d!\n", Pid);
					}

//==========================================================================================================\
============================================================================================================\
============================================================================================================
					bool CM_print=false;
					if(CM_print)
					{
						dir=out_Path+"CorMatrices/"+Date[d]+"/";
						string CorMatrix_FileName=dir+Date[d]+"_CorMatrix_"+current_time+"_"+to_string(info.Channel)+\
						+"_"+to_string(info.PulseFreq)+"_"+to_string(info.PulseLength)+".dat";
						out.open(CorMatrix_FileName, std::ofstream::trunc);
						for(size_t t=0; t<Tlength-L; t++)
						{
							for(size_t tau=0; tau<L; tau++)
							{
								double A[POLY_DEGREE+1];
								vec p=InterPol[it][info][t];
								std::copy(p.begin(), p.end(), A);
								out<<info.Channel<<"\t"<<info.PulseFreq<<"\t"<<info.PulseLength<<"\t";
								out<<double(t)*(1000.0/double(info.DecimationFreq))*LIGHT_VELOCITY/2.0\
									+double(info.TotalDelay)*0.15<<"\t"<<tau*(1000.0/double(info.DecimationFreq))<<"\t";
								out<<real(Matrix[it][info][t][tau])<<"\t";
								out<<polyval(double(tau), p)<<"\t";
								out<<real(R1[it][info][t][tau])<<"\t";
								p=InterPol1[it][info][t];
								out<<polyval(double(tau), p)<<"\n";
							}
						}
						out.close();
					}
				}

				string cmd;
				string dir1="/home/tashlykov/test/", dir2="current/", dir3="archive/";
				if(Nsample[it][info]>1000.0)
				{
					FileName="List_"+to_string(omp_get_thread_num())+".txt";
					out.open(FileName);
					out<<List<<"\n";
					out.close();
					
					inp.open("ConfigFile.cfg");
					string line, str;
					line.clear();
					list.clear();
					size_t found=0, found1=0;
					while(!inp.eof())
					{
						getline(inp, line);
						found=line.find("PathArchive");
						if(found!=-1)
						{
							found1=line.find_first_of("/");
							dir1=line.substr(found1, line.size()-found1);
						}
						found=line.find("PathCurrent");
						if(found!=-1)
						{
							found1=line.find_first_of("/");
							dir2=line.substr(found1, line.size()-found1);
						}
						list+=line+"\n";
					}
					inp.close();
					printf("PathArchive %s\nPathCurrent %s\n", dir1.c_str(), dir2.c_str());
					
					found=list.find("EndFile");
					line=list.substr(found);
					// printf("Pid %d: %s\n\n", omp_get_thread_num(), line.c_str());
					found=line.find("1");
					str=to_string(omp_get_thread_num());
					line.replace(found, line.size()-found, str);
					// printf("\t\t\t\tPid %d: %s\n", omp_get_thread_num(), line.c_str());
					found=list.find("EndFile");
					list.replace(found, list.size()-found, line);
					FileName="ConfigFile_"+to_string(omp_get_thread_num())+".cfg";
					out.open(FileName);
					out<<list<<"\n";
					out.close();
					inp.open(FileName);
					if(!inp.is_open())
						printf("Filing config failed!\n");
					inp.close();
					
					if(it==UT0)
					{
						cmd="rm "+dir1+"current_*";
						system(cmd.c_str());
						
						cmd="rm "+dir2+"*.dat";
						system(cmd.c_str());
					}

					cmd="./ise -s List_"+to_string(omp_get_thread_num())\
						+".txt ConfigFile_"+to_string(omp_get_thread_num())+".cfg";
					printf("\n\t%s. Nsample %d.\n\n", cmd.c_str(), Nsample[it][info]);

					double t0=omp_get_wtime();
					int status;
					if(ThreadIsBuisy[Pid]>0)
					{
						for(int i=1; i<ThreadIsBuisy[Pid]; i++)
							this_thread::sleep_for(std::chrono::seconds(20));
					}
					status=system(cmd.c_str());
					{
						if(status<0)
							std::cerr<<"Error: "<<strerror(errno) << '\n';
						else
						{
							if(WIFEXITED(status))
								std::cout<<"Program returned normally, exit code "<<WEXITSTATUS(status) << '\n';
							else
								std::cout<<"Program exited abnormaly\n";
						}
						printf("Pid %d : Processing time %f!\n", \
							omp_get_thread_num(), omp_get_wtime()-t0);		
					}
				}

				for(size_t m=0; m<Info[it].size(); m++)
				{
					RadarChannelInfo info=Info[it][m];
					if(info.PulseLength<900)
					{
						FileName=dir1+"current_"+to_string(info.PulseFreq+300)+"_"\
							+to_string(info.Channel+1)+"_"+to_string(omp_get_thread_num())+"_dat.dat";
						inp.open(FileName);
						valid[it][info]=inp.is_open();
						if(!valid[it][info])
							fprintf(stderr, "File %s does not exist!\n", FileName.c_str());
						else
						{
							vector<string> buf;
							vec H, N, Pm, Pe;
							getline(inp, line);
							while(!inp.eof())
							{
								getline(inp, line);
								// printf("line       %s\n", line.c_str());
								buf=LineSplit(line, " ");
								if(buf.size()>1)
								{		
									H.push_back(double(round(atof(buf[0].c_str()))*10.0)/10.0);
									Pm.push_back(atof(buf[2].c_str()));
									N.push_back(atof(buf[3].c_str()));
									Pe.push_back(atof(buf[4].c_str()));
								}
							}
							inp.close();
							printf("FileName %s; H size %d; Freq %d, Ch %d, ut %d, pid # %d\n",\
								FileName.c_str(), H.size(), info.PulseFreq, info.Channel, it, omp_get_thread_num());
							if(H.size()==0)
								printf("\n\t============Header not found in %s============\n\n", FileName.c_str());
							
							height_s[it][info]=H;
							P_model[it][info]=Pm;
							Ne[it][info]=N;
							P_exp[it][info]=Pe;

							// else if(it>1)
							// {
							// 	height_s[it][info]=height_s[it-1][info];
							// 	P_model[it][info].resize(height_s[it][info].size());
							// 	Ne[it][info].resize(height_s[it][info].size());
							// 	P_exp[it][info].resize(height_s[it][info].size());
							// }
							// else
							// {
							// 	height_s[it][info].resize(2000);
							// 	P_model[it][info].resize(2000);
							// 	Ne[it][info].resize(2000);
							// 	P_exp[it][info].resize(2000);
							// }
						}
					}
				}
			}
			TimeQuadratures.clear();
			vecc().swap(Q);
			vecc().swap(Q0);
			vecc().swap(Qu);
			#pragma omp barrier
		}

		for(size_t u=0; u<NUMBER_OF_INTERVALS; u++)
		{
			Matrix[u].clear();
			R1[u].clear();
			UT_hist[u].clear();
			Hist[u].clear();
		}

		vector<map<RadarChannelInfo, vector<vecc> > >().swap(Matrix);
		vector<map<RadarChannelInfo, vector<vecc> > >().swap(R1);
		vector<map<RadarChannelInfo, vec> >().swap(UT_hist);
		vector<map<RadarChannelInfo, vec> >().swap(Hist);

		// printf("Matrix size is %d, capacity is %d, address %p\n", \
			Matrix.size(), Matrix.capacity(), &Matrix);

		printf("CorMatrice are done for Pid %d!\n", omp_get_thread_num());

		map<RadarChannelInfo, bool> Xinfo;
		for(size_t u=0; u<NUMBER_OF_INTERVALS; u++)
		{
			for(size_t i=0; i<Info[u].size(); i++)
			{
				if(Info[u][i].PulseLength!=0 and Info[u][i].PulseFreq!=0 and Info[u][i].LengthOfData!=0)
					Xinfo[Info[u][i]]=true;
			}
		}

		RadarChannelInfo info_0_155500_700, info_2_155500_700, info_0_155500_900, info_2_155500_900;
		RadarChannelInfo info_0_159500_700, info_2_159500_700, info_0_159500_900, info_2_159500_900;

		printf("Arranging for info file!\n");
		dir=out_Path;
		// printf("____________dir %s\n", dir.c_str());
		string FileName=dir+Date[d]+"_info.dat";
		ofstream out;
		out.open(FileName, std::ofstream::trunc);
		out<<"Channel_#\t"<<"Frequency\t"<<"Pulse_length\t"<<"Length_of_data\t";
		out<<"Decimation_frequency\t"<<"Total_delay\t"<<"UT\t"<<"Number_of_samples\n";
		for(size_t m=0; m<Xinfo.size(); m++)
		{
			auto j=Xinfo.begin();
			std::advance(j, m);
			RadarChannelInfo info=j->first;
			for(size_t it=0; it<NUMBER_OF_INTERVALS; it++)
			{
				if(Nsample[it].count(info)!=0)
				{
					out<<info.Channel<<"\t"<<info.PulseFreq<<"\t"<<info.PulseLength<<"\t";
					out<<info.LengthOfData<<"\t"<<info.DecimationFreq<<"\t"<<info.TotalDelay<<"\t";
					out<<double(it+UT0)/double(Time_Step)<<"\t";
					out<<Nsample[it][info]<<"\n";
					// out<<real(Matrix[it][info][200][10])<<"\t";
					// out<<InterPol[it][info][200][0]<<"\n";
				}
			}

			if(info.PulseLength==900)
				L_900=900/unsigned(1000.0/double(info.DecimationFreq));
			if(info.PulseLength==700)
				L_700=700/unsigned(1000.0/double(info.DecimationFreq));
		}
		out.close();
		if(L_900==0)
			L_900=L_700;

		printf("Number of signals %d!\n", Xinfo.size());

		for(size_t i=0; i<Xinfo.size(); i++)
		{
			auto j=Xinfo.begin();
			std::advance(j, i);
			RadarChannelInfo info=j->first;
			FileName=dir+Date[d]+"_hist_"+to_string(info.Channel)+"_"+\
				to_string(info.PulseFreq)+"_"+to_string(info.PulseLength)+".dat";
			out.open(FileName, std::ofstream::trunc);
			out<<"UT\t"<<"HistMean\t"<<"HistDistr\t"<<"abs(Kr)\t"<<"arg(Kr)\n";
			printf("%s\n", FileName.c_str());
			// printf("%d %d %d!!!\n", UT_hist[0].size(), MHist[0].size(), DHist[0].size());

			for(size_t it=0; it<NUMBER_OF_INTERVALS; it++)
			{
				// printf("%d %d %d!!!\n", UT_hist[it][info].size(), MHist[it][info].size(), DHist[it][info].size());
				if(Nsample[it].count(info)!=0)
				{
					for(size_t j=0; j<Nsample[it][info]; j++)
					{
						// out<<double(it+UT0)/double(Time_Step)<<"\t";
						// out<<UT_hist[it][info][j]<<"\t";
						out<<j<<"\t";
						out<<MHist[it][info][j]<<"\t"<<DHist[it][info][j]<<"\t";
						out<<abs(Kr[it][info][j])<<"\t"<<arg(Kr[it][info][j])<<"\n";
					}
				}
			}
			out.close();
		}

		// printf("info_0_900 %d %d %d by #%d\n", info_0_900.Channel, info_0_900.PulseFreq, info_0_900.PulseLength, omp_get_thread_num());
		// printf("info_2_900 %d %d %d by #%d\n", info_2_900.Channel, info_2_900.PulseFreq, info_2_900.PulseLength, omp_get_thread_num());
		
		// omp_set_num_threads(Xinfo.size());
		omp_set_num_threads(1);
		#pragma omp parallel
		{
			ofstream out, sec;
			#pragma omp for
			for(size_t m=0; m<Xinfo.size(); m++)
			{
				auto j=Xinfo.begin();
				std::advance(j, m);
				RadarChannelInfo info=j->first;
				unsigned Tlength=info.LengthOfData/4;
				double deltat=1000.0/double(info.DecimationFreq);
				unsigned L=floor(double(info.PulseLength)/deltat);
				printf("Tlength = %d, deltat = %f\n", Tlength, deltat);
				vec height(Tlength-L_900), height1(Tlength-L_900);
				for(size_t r=0; r<Tlength-L_900; r++)
				{
					double range=double(r)*(1000.0/double(info.DecimationFreq))*LIGHT_VELOCITY/2.0\
						+double(info.TotalDelay)*0.15;
					double Eps=AzimMaxDNR(double(info.PulseFreq)), fi, lg;
					Found_FI_LG_H__for_Ep_Gam_R(Eps, 0.0, range-double(info.PulseLength)*0.15, \
						nr_lat*Rad, nr_lg*Rad, 0.0, fi, lg, height[r]);
					Found_FI_LG_H__for_Ep_Gam_R(Eps, 0.0, range-double(info.PulseLength)*0.15/2.0, \
						nr_lat*Rad, nr_lg*Rad, 0.0, fi, lg, height1[r]);
				}

				string FileName=dir+Date[d]+"_temps_"+to_string(info.Channel)+"_"+\
					to_string(info.PulseFreq)+"_"+to_string(info.PulseLength)+".dat";
				out.open(FileName, std::ofstream::trunc);
				FileName=dir+Date[d]+"_tempsec_"+to_string(info.Channel)+"_"+\
					to_string(info.PulseFreq)+"_"+to_string(info.PulseLength)+".dat";
				sec.open(FileName, std::ofstream::trunc);
				printf("Arranging for temp file! %d %d %d\n", info.Channel, info.PulseFreq, info.PulseLength);
				// printf("InterPol size [%d, %d, %d]\n", InterPol1.size(), InterPol1[it][info].size(), InterPol1[it][info][0].size());
				// printf("Err size [%d, %d, %d]\n", Err1.size(), Err1[it][info].size(), Err1[it][info][0].size());
				// printf("Chi size [%d, %d]\n", Chi1.size(), Chi1[it][info].size());
				out<<"Date\t"<<"H:M\t";
				out<<"Height\t"<<"S/N\t"<<"Te\t"<<"Tr\t"<<"t_0\t"<<"t_min\t"<<"a_min\t"<<"Noise\t"<<"Misfit\t";
				out<<"Height_acf\t"<<"S/N_acf\t"<<"Te_acf\t"<<"Tr_acf\t"<<"Te_fit\t"<<"Tr_fit\t";
				out<<"t_0_acf\t"<<"t_min_acf\t"<<"a_min_acf\t"<<"Noise_acf\t"<<"Misfit_acf\t";
				out<<"Height_EST\t"<<"S/N_EST\t"<<"Te_EST\t"<<"Tr_EST\t"<<"t_0_EST\t"<<"t_min_EST\t"<<"a_min_EST\t";
				out<<"Ne\t"<<"P_exp\t"<<"P_model\t"<<"P_conv\t"<<"\n";
				for(size_t u=0; u<NUMBER_OF_INTERVALS; u++)
				{
					double ut=double(u+UT0)/double(Time_Step);
					if(Nsample[u][info]>100)
					{
						Norm[u][info]=0.0;
						vec Te, Tr, W, Te_fit, Tr_fit, Te_err, Tr_err;
						for(size_t r=0; r<Tlength-L_900; r++)
						{
							double t_0, t_min, a_min, te, tr;
							double SNR=1.0;
							SNR=(RawPower[u][info][r]/InterPol1[u][info][r][0]-1.0);
							ACF_pardet(t_0, t_min, a_min, L, Lag_Truncate, InterPol1[u][info][r]);
							Get_Temperatures(te, tr, t_0, t_min, a_min);
							Te.push_back(te);
							Tr.push_back(tr);
							if(height[r]<170.0 or height[r]>400.0 or te<400.0 or te>3500.0 or tr<1.0 or tr>3.5) 
								W.push_back(0.0);
							else
								W.push_back(SNR);
							if(SNR>0.02 and height1[r]>200.0 and height1[r]<400.0)
							{
								Te_mean[u][info]+=te;
								Tr_mean[u][info]+=tr;
								Norm[u][info]+=1.0;
							}
						}
						B_spline_fit(Te_fit, Te_err, height, Te, W);
						B_spline_fit(Tr_fit, Tr_err, height, Tr, W);
						
						sec<<Date[d]<<"\t"<<floor(ut)<<":"<<round(60.0*(ut-floor(ut)))<<"\t";
						sec<<Te_mean[u][info]/Norm[u][info]<<"\t"<<Tr_mean[u][info]/Norm[u][info]<<"\n";

						printf("Heights scaling!\n");
						if(info.PulseLength>700 or height_s[u][info].size()==0)
							valid[u][info]=false;
						if(valid[u][info])
						{
							vector<bool> isempty;
							isempty.resize(height.size());
							double r=height[0];
							printf("Ne size %d, H size %d, F %d, Ch %d, ut %d, pid # %d!\n",\
								 Ne[u][info].size(), height_s[u][info].size(), info.PulseFreq, info.Channel, u, omp_get_thread_num());
							while(r<height_s[u][info][0])
							{
								// printf("Ne size %d, r = %f!\n", Ne[u][info].size(), r);
								Ne[u][info].insert(Ne[u][info].begin(), 0.0/0.0);
								P_model[u][info].insert(P_model[u][info].begin(), 0.0/0.0);
								P_exp[u][info].insert(P_exp[u][info].begin(), 0.0/0.0);
								r+=0.6;
							}

							r=height_s[u][info].back();
							cout<<height_s[u][info].back()<<"\t"<<height.back()<<endl;
							while(Ne[u][info].size()<height.size())
							{
								// printf("Ne size %d, r = %f!\n", Ne[u][info].size(), r);
								Ne[u][info].push_back(0.0/0.0);
								P_model[u][info].push_back(0.0/0.0);
								P_exp[u][info].push_back(0.0/0.0);	
								r+=0.6;
							}
							printf("Ne size %d, r = %f!\n", Ne[u][info].size(), r);
							for(size_t t=0; t<P_model[u][info].size(); t++)
							{
								if(isnormal(P_model[u][info][t]))
									isempty[t]=false;
								else
									isempty[t]=true;
								// isempty[t]=false;
							}
							P_conv[u][info]=Radeq(P_model[u][info], height, isempty, double(L));

							int n=0;
							double sumP, sumF, sumPF, sumF2;
							sumP=sumF=sumPF=sumF2=0.0;
							// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							size_t TEND=736;
							for(size_t r=0; r<TEND; r++)
							{
								// double N=RawPower[u][info][r]-real(InterPol[u][info][r][0]);
								// if(r>Ground_Clutter_End+L/2 and height[r]<height_s[u][info].back()-2*L)
								if(!isempty[r])
								{
									double SNR=1.0/(RawPower[u][info][r]/InterPol[u][info][r][0]-1.0);
									sumF+=P_conv[u][info][r];
									sumP+=SNR;
									sumPF+=P_conv[u][info][r]*SNR;
									sumF2+=P_conv[u][info][r]*P_conv[u][info][r];
									n++;
								}
							}
							// fprintf(stderr, "SUM: %f, %f, %f, %f\n", sumF, sumP, sumPF, sumF2);
							Aconst[u][info]=(sumPF-sumP*sumF/double(n))/\
								(sumF2-sumF*sumF/double(n));
							Bconst[u][info]=(sumP-Aconst[u][info]*sumF)/double(n);
							vector<bool>().swap(isempty);
						}
						printf("Done!\n");

						for(size_t r=0; r<Tlength-L_900; r++)
						{
							double t_0=0.0, t_min=0.0, a_min=0.0, te=0.0, tr=0.0;
							double SNR=0.0, N=0.0;
							// ACF_pardet(t_0, t_min, a_min, L, InterPol[u][info_2_900][r]);
							// Get_Temperatures(Te2, Tr2, t_0, t_min, a_min);

							out<<Date[d]<<"\t"<<floor(ut)<<":"<<round(60.0*(ut-floor(ut)))<<"\t";
							SNR=1.0/(RawPower[u][info][r]/InterPol[u][info][r][0]-1.0);
							N=RawPower[u][info][r]-InterPol[u][info][r][0];
							ACF_pardet(t_0, t_min, a_min, L, Lag_Truncate, InterPol[u][info][r]);
							t_0*=deltat;
							t_min*=deltat;
							Get_Temperatures(te, tr, t_0, t_min, a_min);
							if(N<1.2*Noise[u][info][0] and N>0.8*Noise[u][info][0])
								out<<height[r]<<"\t"<<SNR<<"\t"<<te<<"\t"<<tr<<"\t"<<t_0<<"\t"<<t_min<<"\t"<<a_min<<"\t";
							else
								out<<height[r]<<"\t"<<SNR<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t";
							out<<N<<"\t"<<Chi[u][info][r]<<"\t";

							SNR=1.0/(RawPower[u][info][r]/InterPol1[u][info][r][0]-1.0);
							ACF_pardet(t_0, t_min, a_min, L, Lag_Truncate, InterPol1[u][info][r]);
							t_0*=deltat;
							t_min*=deltat;
							Get_Temperatures(te, tr, t_0, t_min, a_min);
							if(N<1.2*Noise[u][info][0] and N>0.8*Noise[u][info][0])
							{
								out<<height1[r]<<"\t"<<SNR<<"\t"<<te<<"\t"<<tr<<"\t";
								out<<Te_fit[r]<<"\t"<<Tr_fit[r]<<"\t"<<t_0<<"\t"<<t_min<<"\t"<<a_min<<"\t";
							}
							else
								out<<height1[r]<<"\t"<<W[r]<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t";
							out<<N<<"\t"<<Chi1[u][info][r]<<"\t";

							RadarChannelInfo info1=info;
							RadarChannelInfo info2=info;
							info1.PulseLength=900;
							info1.PulseLength=700;
							vec p1=InterPol[u][info1][r];
							vec p2=InterPol[u][info2][r];
							vec p(InterPol[u][info1][r].size());
							for(size_t i=0; i<p.size(); i++)
								p[i]=p1[i]-p2[i];
							if(InterPol[u].count(info1)>0 and InterPol[u].count(info2)>0)
							{
								SNR=1.0/((RawPower[u][info1][r]-RawPower[u][info2][r])/(InterPol[u][info1][r][0]-InterPol[u][info2][r][0])-1.0);
								ACF_pardet(t_0, t_min, a_min, 0, Lag_Truncate, p);
								t_0*=deltat;
								t_min*=deltat;
								Get_Temperatures(te, tr, t_0, t_min, a_min);
								if(SNR>0.01)
								{
									out<<height[r]<<"\t"<<SNR<<"\t"<<te<<"\t"<<tr<<"\t";
									out<<t_0<<"\t"<<t_min<<"\t"<<a_min<<"\n";
								}
								else
									out<<height[r]<<"\t"<<SNR<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t";
							}
							else
							{
								out<<height[r]<<"\t"<<SNR<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t";
								printf("Wrong pulse length!\n\n");
							}
							if(valid[u][info])
							{
								out<<Ne[u][info][r]<<"\t"<<P_exp[u][info][r]<<"\t"<<P_model[u][info][r]<<"\t";
								// if(r>Ground_Clutter_End+L/2 and height[r]<height_s[u][info].back()-2*L)
									// out<<P_conv[u][info][r]<<"\t"<<Aconst[u][info]<<"\t"<<Bconst[u][info]<<"\n";
									out<<Aconst[u][info]*P_conv[u][info][r]+Bconst[u][info]<<"\n";
								// else out<<"nan\n";
							}
							else
								out<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\n";
						}
					}
					else
					{
						for(size_t r=0; r<Tlength-L_900; r++)
						{
							out<<Date[d]<<"\t"<<floor(ut)<<":"<<round(60.0*(ut-floor(ut)))<<"\t";
							out<<height[r]<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t";
							out<<height1[r]<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t";
							out<<"nan\t"<<"nan\t"<<"nan\t";
							out<<height[r]<<"\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\t"<<"nan\n";
						}
					}					
				}
				out.close();
				sec.close();
				printf("Done! %d %d %d\n\n", info.Channel, info.PulseFreq, info.PulseLength);

				#pragma omp_barrier
				#pragma omp_critical
				{
					printf("Arranging for polynom file! %d %d %d\n", info.Channel, info.PulseFreq, info.PulseLength);
					FileName=dir+Date[d]+"_polynom_"+to_string(info.Channel)+"_"+\
					to_string(info.PulseFreq)+"_"+to_string(info.PulseLength)+".dat";;
					out.open(FileName, std::ofstream::trunc);
					out<<"Date\t"<<"UT\t"<<"Range\t";
					for(size_t p=0; p<POLY_DEGREE; p++)
						out<<"Polynom_coef_value\t"<<"Polynom_coef_err\t";
					out<<"Polynom_misfit_value\n";
					double deltat=1000.0/double(info.DecimationFreq);
					unsigned L=floor(double(info.PulseLength)/deltat);	
					for(size_t u=0; u<NUMBER_OF_INTERVALS; u++)
					{
						// printf("%d\n", u);
						if(Nsample[u][info]>100)
						{	
							for(size_t t=0; t<Tlength-L; t++)
							{
								out<<Date[d]<<"\t"<<double(u+UT0)/double(Time_Step)<<"\t";
								out<<double(t)*deltat*LIGHT_VELOCITY/2.0+double(info.TotalDelay)*0.15<<"\t";
								for(size_t p=0; p<POLY_DEGREE; p++)
									out<<InterPol1[u][info][t][p]<<"\t"<<Err1[u][info][t][p]<<"\t";
								out<<Chi1[u][info][t]<<"\n";
							}
						}
						else
						{
							for(size_t t=0; t<Tlength-L; t++)
							{
								out<<Date[d]<<"\t"<<double(u+UT0)/double(Time_Step)<<"\t";
								out<<double(t)*deltat*LIGHT_VELOCITY/2.0+double(info.TotalDelay)*0.15<<"\t";
								for(size_t p=0; p<POLY_DEGREE; p++)
									out<<"nan\t"<<"nan\t";
								out<<"nan\n";
							}
						}
					}
					out.close();
					printf("Done! %d %d %d\n", info.Channel, info.PulseFreq, info.PulseLength);
				}
			}
			#pragma omp barrier
		}

		for(size_t u=0; u<NUMBER_OF_INTERVALS; u++)
		{
			Chi[u].clear();
			Chi1[u].clear();
			RawPower[u].clear();
			Nsample[u].clear();
			Nsample2[u].clear();
			Err[u].clear();
			Err1[u].clear();
			InterPol[u].clear();
			InterPol1[u].clear();
			Info[u].clear();
		}

		ListMap.clear();
		vector<map<RadarChannelInfo, vec> >().swap(Chi);
		vector<map<RadarChannelInfo, vec> >().swap(Chi1);
		vector<map<RadarChannelInfo, vec> >().swap(RawPower);
		vector<map<RadarChannelInfo, vec> >().swap(Noise);
		vector<map<RadarChannelInfo, vec> >().swap(UT_hist);
		vector<map<RadarChannelInfo, vec> >().swap(Hist);
		vector<map<RadarChannelInfo, vec> >().swap(mHist);
		vector<map<RadarChannelInfo, vec> >().swap(dHist);
		vector<map<RadarChannelInfo, vector<vec> > >().swap(Err);
		vector<map<RadarChannelInfo, vector<vec> > >().swap(Err1);
		vector<map<RadarChannelInfo, vector<vec> > >().swap(InterPol);
		vector<map<RadarChannelInfo, vector<vec> > >().swap(InterPol1);
		vector<map<RadarChannelInfo, vector<vecc> > >().swap(Matrix);
		vector<map<RadarChannelInfo, vector<vecc> > >().swap(R1);
		vector<map<RadarChannelInfo, unsigned> >().swap(Nsample);
		vector<map<RadarChannelInfo, unsigned> >().swap(Nsample2);
		vector<vector<RadarChannelInfo> >().swap(Info);
		vector<map<RadarChannelInfo, double> >().swap(Te_mean);
		vector<map<RadarChannelInfo, double> >().swap(Tr_mean);
		vector<map<RadarChannelInfo, double> >().swap(Norm);
		vector<map<RadarChannelInfo, double> >().swap(Aconst);
		vector<map<RadarChannelInfo, double> >().swap(Bconst);
		vector<map<RadarChannelInfo, vecc> >().swap(Kr);
		vector<map<RadarChannelInfo, vec> >().swap(MHist);
		vector<map<RadarChannelInfo, vec> >().swap(DHist);
		vector<map<RadarChannelInfo, vec> >().swap(height_s);
		vector<map<RadarChannelInfo, vec> >().swap(P_model);
		vector<map<RadarChannelInfo, vec> >().swap(P_exp);
		vector<map<RadarChannelInfo, vec> >().swap(P_conv);
		vector<map<RadarChannelInfo, vec> >().swap(Ne);
		vector<map<RadarChannelInfo, bool> >().swap(valid);
	}

	printf("Processing time %f!\n", omp_get_wtime()-init_time);
}
