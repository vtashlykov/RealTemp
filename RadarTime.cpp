#include"RadarTime.h"

RadarTime::RadarTime(){
	Empty=true;
}

void RadarTime::InitMonthes(){
	monthes[0]=31;
	if(Year%4==0) monthes[1]=29; else monthes[1]=28; 
	monthes[2]=31;
	monthes[3]=30;
	monthes[4]=31;
	monthes[5]=30;
	monthes[6]=31;
	monthes[7]=31;
	monthes[8]=30;
	monthes[9]=31;
	monthes[10]=30;
	monthes[11]=31;
}

int RadarTime::GetDaysInYear(){
	InitMonthes();
	int i=0, days=0;
	for(i=0;i<12;i++)days+=monthes[i];
	return days;
}

double RadarTime::GetHour(){
	return double(Hour)+double(Min)/60.0+double(Sec)/(60.0*60.0)+double(mSec)/(60.0*60.0*1000);
}

void RadarTime::PutHour(double H){

}

double RadarTime::GetDayOfYear(){
	InitMonthes();
	int i=0, doy=0;
	for(i=0;i<Month-1;i++){
		doy=doy+monthes[i];
	}
	doy=doy+Day-1;
	return double(doy)+double(Hour)/24.0+double(Min)/(24.0*60.0)+double(Sec)/(24.0*60.0*60.0)+double(mSec)/(24.0*60.0*60.0*1000);
}

void RadarTime::PutDayOfYear(double D){

}

double RadarTime::GetYear(){
	int days=GetDaysInYear();
	double doy=GetDayOfYear();
	return double(Year)+doy/double(days);
//	return double(Year)+doy/double(days)+double(Hour)/(double(days)*24.0)+double(Min)/(double(days)*24.0*60.0)+double(Sec)/(double(days)*24.0*3600.0)+double(mSec)/(double(days)*24.0*3600.0*1000.0);
}

void RadarTime::PutYear(double Y){
	Year=int(Y);
	int days=GetDaysInYear();
	double doy=((Y-int(Y))*double(days));
	int m=0,d=0;
	while(m<12&&d<(int(doy)+1)){
		d+=monthes[m];
		m++;
	}
	Month=m;
	if(m>0)d=d-monthes[m-1];
	else d=0;
	Day=int(doy)-d+1;
	Hour=int((doy-int(doy))*24.0);
	Min=int((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0);
	Sec=int((((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0)-int((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0))*60.0);
	mSec=int(((((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0)-int((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0))*60.0
					-int((((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0)-int((((doy-int(doy))*24.0)-int((doy-int(doy))*24.0))*60.0))*60.0)
					)*1000.0);
	Empty=false;
}


int getNumberOfDays(int month, int year)
{
	if( month == 2)
	{
		if((year%400==0) or (year%4==0 && year%100!=0))
			return 29;
		else
			return 28;
	}
	else if(month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month==12)
		return 31;
	else
		return 30;
}

int getDoy(int year, int month, int day)
{
    int d=0;
    for(int i=1; i<month; i++)
        d+=getNumberOfDays(i, year);
    return d+day;
}

void getDate(int Year, int doY, int &month, int &day)
{
	day=doY;
    int d=31, m=1;
    while(d<doY)
    {
        d+=getNumberOfDays(m, Year);
		day-=getNumberOfDays(m, Year);
		m++;
    }
	month=m;
}

double Ap_to_Kp(double Ap)
{
    double Kp;
    double ap_table[]={0, 2, 3, 4, 5, 6, 7, 9, 12,\
        15, 18, 22, 27, 32, 39, 48, 56, 67,\
		80, 94, 111, 132, 154, 179, 207, 236, 300};
    for(int i=0; i<26; i++)
    {
        if(Ap>=ap_table[i] and Ap<ap_table[i+1])
			Kp=double(i)/3.0;
		else Kp=0.0/0.0;
    }
    if(Ap>=ap_table[26]) Kp=9.0;
    return Kp;
}