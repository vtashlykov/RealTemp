#ifndef RADARTIME_H
#define RADARTIME_H

#include<iostream>
using namespace std;

class RadarTime{
	private:
		int monthes[12];
		void InitMonthes();

	public:
		RadarTime();
		bool Empty;
		int Year;
		int Month;
		int Day;
		int Hour;
		int Min;
		int Sec;
		int mSec;
		int GetDaysInYear();
		double GetYear();
		void PutYear(double);
		double GetDayOfYear();
		void PutDayOfYear(double);
		double GetHour();
		void PutHour(double);
};

int getNumberOfDays(int month, int year);

int getDoy(int year, int month, int day);

void getDate(int Year, int doY, int &month, int &day);

double Ap_to_Kp(double Ap);

#endif
