#include "InverseProblemsolver.h"

bool operator< ( Plasma_pars a, Plasma_pars b )
{
	return std::make_pair(a.Te, a.Tr) < std::make_pair(b.Te, b.Tr) ;
}

double Pulse(double t, double T)
{
	if(t>=0.0 and t-T<0.0) return 1.0;
	else return 0.0;
}

double Ne_profile(double x, void * params)
{
	struct Ne_pars Ch=*(struct Ne_pars *) params;
	double z;
	if(x<=Ch.Hmf2)
		z=(x-Ch.Hmf2)/Ch.H_b;
	else
		z=(x-Ch.Hmf2)/Ch.H_t;
	return Ch.Nmf2*exp(1.0-z-exp(-z));
}

vec Faraday(vec height, std::vector<bool> isempty, void * params)
{
	struct Ne_pars Ch=*(struct Ne_pars *) params;
	vec Far;
	double Omega=0.0, error, P_max=0.0;
	size_t neval, n;
	n=height.size();
	Far.resize(n);
	gsl_function F;
	F.function=&Ne_profile;
	F.params=&Ch;
	for(size_t r=0; r<n; r++)
	{
		gsl_integration_qng(&F, 0.0, height[r], 0.0, 1.0e-4, &Omega, &error, &neval);
		Far[r]=Ne_profile(height[r], &Ch)*pow(cos(1.0e-3*Omega+double(Ch.Phi)*M_PI/180.0)/height[r], 2.0);
        if(Far[r]>P_max and !isempty[r])
            P_max=Far[r];
    }
    for(size_t r=0; r<n; r++)
        Far[r]/=P_max;
	return Far;
}

vec Power_v(vec Ne)
{
	vec P;
	double P_max=0.0, Omega=0.0;
	P.resize(Ne.size());
	for(size_t r=0; r<Ne.size(); r++)
	{
		for(size_t z=0; z<r; z++)
			Omega+=Ne[z];
		double h=100.0+double(r)*10.0;
		P[r]=Ne[r]*pow(cos(1.0e-3*Omega)/h, 2.0);
	}
	return P;
}

vec Radeq(vec Far, vec height, std::vector<bool> isempty, double Length)
{
	vec P;
    double Pmax=0.0;
	size_t n=height.size();
	P.resize(n);
	for(size_t r=0; r<n; r++)
	{
		for(size_t z=0; z<n; z++)
		{
			if(!isempty[z] and isnormal(Far[z]))
				P[r]+=Far[z]*Pulse(height[z]-height[r], Length);
		}
        if(P[r]>Pmax and !isempty[r])
            Pmax=P[r];
	}
    for(size_t r=0; r<n; r++)
        P[r]/=Pmax;
	// fprintf(stderr, "Convolution Length %f\n", Length);
	// for(size_t r=0; r<n; r++)
	// 	cerr<<r<<"\t"<<height[r]<<"\t"<<Far[r]<<"\t"<<P[r]<<"\t"<<isempty[r]<<"\n";
	return P;
}

double Misfit(vec f, vec Y)
{
	double Chisq=0.0;
	for(size_t i=0; i<f.size(); i++)
		Chisq+=pow(f[i]-Y[i], 2.0);
	return Chisq/double(f.size());
}

/* ------------------Nonlinear least squares fitting------------------*/

double ObjFunc(const vec &x, vec &grad, void *f_data)
{
    vec data, height;
	double pulselength;
	std::vector<bool> isempty;
	OF_data* D;
	D=(OF_data*)f_data;
	data=D->data;
	height=D->height;
	isempty=D->isempty;
	pulselength=D->pulselength;

    double chisq=0.0;
    struct Ne_pars Ch={x[0], x[1], x[2], x[3], x[4]};
	vec Far=Faraday(height, isempty, &Ch);
	vec P_model=Radeq(Far, height, isempty, pulselength);
	// vec P_model=Pol;

	double sumP, sumF, sumPF, sumF2, A, B;
	sumP=sumF=sumPF=sumF2=0.0;
	int n=0;
	// printf("Data size is %d:%d:%d for process #%d!\n", D->height.size(), D->isempty.size(), D->data.size(), omp_get_thread_num());
	for(size_t r=0; r<data.size(); r++)
	{
		if(!isempty[r] and isnormal(data[r]))
		{
			sumF+=P_model[r];
			sumP+=data[r];
			sumPF+=P_model[r]*data[r];
			sumF2+=P_model[r]*P_model[r];
			n++;
		}
	}
	A=(sumPF-sumP*sumF/double(n))/(sumF2-sumF*sumF/double(n));
	B=(sumP-A*sumF)/double(n);

    for(size_t r=0; r<data.size(); r++)
	{
		if(!isempty[r] and isnormal(data[r]))
			chisq+=pow((data[r]-B)/A-P_model[r], 2.0);
	}

	// printf("Const: %f\t%f\t%f\t%f\t%f\t%f\t%f\n", sumF, sumP, sumPF, sumF2, A, B, sqrt(chisq)/n);

	// if(D->Global_solver_on)
	// {
	// 	fprintf(stderr, "%d\t0\t%f\t%f\t%f\t%f\t%f\t%f\n",\
	// 		omp_get_thread_num(), x[0], x[1], x[2], x[3], x[4], sqrt(chisq)/n);
	// }

	*(D->out)<<x[0]<<"\t"<<x[1]<<"\t"<<x[2]<<"\t"<<x[3]<<"\t"<<x[4]<<"\t"<<sqrt(chisq)/n<<"\n";

	/* std::string e="Infinite misfit!\n";
	if(!std::isnormal(chisq))
		throw e; */
    return sqrt(chisq)/n;
}

double GenObjFunc(const vec &x, vec &grad, void *f_data)
{
	std::vector<OF_data>* D;
	D=(std::vector<OF_data>*)f_data;
	OF_data data=D->at(0);
	double MisFit=0.0;
	for(size_t i=0; i<D->size(); i++)
	{
		data=D->at(i);
		MisFit+=ObjFunc(x, grad, &data);
	}
	return MisFit/double(D->size());
}

/*--------------------------INVERSE PROBLEM SOLVER----------------------------------------*/
void Global_solver(OF_data Data, double UT, \
	struct Optim_option &option, vec &x, std::string &mesout)
{
	double current_time=omp_get_wtime();
	char message[1000];
	vec lb(NUMBER_OF_PARAMETERS), ub(NUMBER_OF_PARAMETERS);
	lb[0]=10.0;
	lb[1]=200.0;
    lb[2]=30.0;
    lb[3]=30.0;
	lb[4]=0.0;
	ub[0]=200.0;
	ub[1]=450.0;
    ub[2]=240.0;
    ub[3]=180.0;
	ub[4]=200.0;
	// nlopt::opt optg(nlopt::G_MLSL_LDS, NUMBER_OF_PARAMETERS);//GN_CRS2_LM
	nlopt::opt optg(nlopt::GN_ISRES, NUMBER_OF_PARAMETERS);//GN_CRS2_LM
	// nlopt::opt optg(nlopt::GN_CRS2_LM, NUMBER_OF_PARAMETERS);
	nlopt::opt optl(nlopt::LN_COBYLA, NUMBER_OF_PARAMETERS);
	optg.nlopt::opt::set_local_optimizer(optl);
	optg.nlopt::opt::set_xtol_rel(option.xtol_rel);
	optg.nlopt::opt::set_ftol_rel(option.ftol_rel);
	optg.nlopt::opt::set_lower_bounds(lb);
	optg.nlopt::opt::set_upper_bounds(ub);
	optg.nlopt::opt::set_maxeval(option.maxeval);
	// optl.nlopt::opt::set_xtol_abs(1.0);
	// optg.nlopt::opt::set_stopval(5.0);
	// printf("Data size is %d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", Data.size());
	// OF_data data=Data;
	// printf("Data size is %d:%d:%d for process #%d!\n", data.height.size(), data.isempty.size(), data.data.size(), omp_get_thread_num());

	// printf("Height:Data:PulseLength=%d:%d:%d\n", \
		data.height.size(), data.data.size(), data.pulselength);
	optg.set_min_objective(ObjFunc, &Data);
	try
	{
		nlopt::result result=optg.optimize(x, option.minf);
		// nlopt::algorithm algorithm=optg.get_algorithm();
		// optg.nlopt::opt::get_algorithm_name()
		// const char *nlopt_result_to_string(nlopt_result result);
		// option.neval=optg.nlopt::opt::get_numevals();
		option.time=omp_get_wtime()-current_time;
		option.opt_failed=false;
		// option.stopval=optg.nlopt::opt::get_stopval();
		// option.ftol_rel=optg.nlopt::opt::get_ftol_rel();
		// option.ftol_abs=optg.nlopt::opt::get_ftol_abs();
		// option.xtol_rel=optg.nlopt::opt::get_xtol_rel();

		sprintf(message, "%2.2f\t#%2.2d | %d %d | UT=%2.2d:%2.2d || %3.0f | %3.0f | %3.0f | %3.0f | %3.0f || %1.4e | %3.1f | %4.0d | %.1e| %.1e | %d || GN\n",
			UT, omp_get_thread_num(), Data.frequency, Data.pulselength,\
			int(floor(UT)), int(round((UT-floor(UT))*60.0)),\
			floor(x[0]), floor(x[1]), floor(x[2]), floor(x[3]), floor(x[4]),\
			option.minf, option.time, optg.nlopt::opt::get_numevals(),\
			optg.nlopt::opt::get_ftol_rel(), optg.nlopt::opt::get_xtol_rel(), \
			result);
		mesout=mesout+message;
	}
	catch(std::exception &e)
	{
		fprintf(stdout, "Thread #%2.2d (UT=%2.2d:%2.2d), Local optimization failed:%s, %d\n",
			omp_get_thread_num(), int(floor(UT)), int(round((UT-floor(UT))*60.0)), e.what(), optg.nlopt::opt::get_numevals());
		option.opt_failed=true;
		for(int i=0; i<NUMBER_OF_PARAMETERS; i++)
			x[i]=0.0/0.0;
	}
}

void Local_solver(OF_data Data, double UT, \
	struct Optim_option &option, vec &x, std::string &mesout)
{
	double current_time=omp_get_wtime();
	char message[1000];
	vec lb(NUMBER_OF_PARAMETERS), ub(NUMBER_OF_PARAMETERS);
	lb[0]=10.0;
	lb[1]=200.0;
    lb[2]=30.0;
    lb[3]=30.0;
	lb[4]=0.0;
	ub[0]=200.0;
	ub[1]=450.0;
    ub[2]=240.0;
    ub[3]=180.0;
	ub[4]=180.0;
	nlopt::opt optl(nlopt::LN_COBYLA, NUMBER_OF_PARAMETERS);//LN_BOBYQA
	optl.nlopt::opt::set_lower_bounds(lb);
	optl.nlopt::opt::set_upper_bounds(ub);
	optl.nlopt::opt::set_xtol_rel(option.xtol_rel);
	optl.nlopt::opt::set_ftol_rel(option.ftol_rel);
	optl.nlopt::opt::set_maxeval(option.maxeval);
	optl.nlopt::opt::set_min_objective(ObjFunc, &Data);
	try
	{
		nlopt::result result=optl.optimize(x, option.minf);
		// option.neval=optl.nlopt::opt::get_numevals();
		option.time=omp_get_wtime()-current_time;
		option.opt_failed=false;

		sprintf(message, "%2.2f\t#%2.2d | %d %d | UT=%2.2d:%2.2d || %3.0f | %3.0f | %3.0f | %3.0f | %3.0f || %1.4e | %3.1f | %4.0d | %.1e | %.1e| %d || LN\n",
			UT, omp_get_thread_num(), Data.frequency, Data.pulselength,\
			int(floor(UT)), int(round((UT-floor(UT))*60.0)),\
			floor(x[0]), floor(x[1]), floor(x[2]), floor(x[3]), floor(x[4]),\
			option.minf, option.time, optl.nlopt::opt::get_numevals(),\
			optl.nlopt::opt::get_ftol_rel(), optl.nlopt::opt::get_xtol_rel(), \
			result);
		mesout=mesout+message;
	}
	catch(std::exception &e)
	{
		fprintf(stdout, "Thread #%2.2d (UT=%2.2d:%2.2d), Local optimization failed:%s, %d\n",
			omp_get_thread_num(), int(floor(UT)), int(round((UT-floor(UT))*60.0)), e.what(), optl.nlopt::opt::get_numevals());
		option.opt_failed=true;
		for(int i=0; i<NUMBER_OF_PARAMETERS; i++)
			x[i]=0.0/0.0;
	}
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Temperature functions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

std::map<Plasma_pars, vec> Scatter_section()
{
	std::ifstream inp;
	std::map<Plasma_pars, vec> M;
	Plasma_pars p;
	std::string head;
	vec R;
	int count=0;
	double Te, Tr, lambda=3.0e5/155500.0;
	inp.open("Scatter_section.dat");
	while(!inp.eof())
	{
		getline(inp, head);
		vector<string> buf=LineSplit(head, "\t");
		Tr=atof(buf[0].c_str());
		Te=atof(buf[1].c_str());
		R.push_back(atof(buf[3].c_str()));
		count++;
		if(count==112)
		{
			count=0;
			p={lambda, Te, Tr};
			M[p]=R;
			R.clear();
		}
	}
	inp.close();
	return M;
}

vec T_profile(vec coefs, vec knots, vec height)
{
    vec Te, b=coefs;
	const size_t nbreak=knots.size();
    const size_t ncoeffs=coefs.size();
    gsl_bspline_workspace *bw;
    gsl_vector *B, *K;
    /* allocate a cubic bspline workspace (k=4) */
    bw=gsl_bspline_alloc(4, nbreak);
    B=gsl_vector_alloc(ncoeffs);
    K=gsl_vector_alloc(nbreak);
    // printf("Ncoeffs=%d; Nbreaks=%d\n", ncoeffs, nbreak);
    for(size_t i=0; i<nbreak; i++)
    {
        gsl_vector_set(K, i, knots[i]);
        // std::cout<<gsl_vector_get(K, i)<<"\t";
    }
    // std::cout<<"\n";
    // for(size_t k=0; k<ncoeffs; k++)
    // {
    //     std::cout<<b[k]<<"\t";
    // }
    // std::cout<<"\n";

    gsl_bspline_knots(K, bw);
    // gsl_bspline_knots_uniform(height, height[height.size()-1], bw);
    for(size_t r=0; r<height.size(); r++)
    {
        double te=0.0;
        if(height[r]>knots[0] and height[r]<knots[nbreak-1])
        {
            gsl_bspline_eval(height[r], B, bw);
            for(size_t k=0; k<ncoeffs; k++)
                te+=b[k]*gsl_vector_get(B, k);
        }
        else
            te=sqrt(-1.0);

        Te.push_back(te);
    }
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    return Te;
}

std::vector<vec> Radeq_ACF(std::vector<vec> B0, vec Far, vec height, size_t Length)
{
	// printf("[%d:%d], %d, %d, %d\n", B0.size(), B0[0].size(), Far.size(), height.size(), Length);
	std::vector<vec> B;
	B.resize(height.size(), vec (Length));
	size_t n=B.size();
    vec norm(n);
	for(size_t tau=0; tau<Length; tau++)
	{
		for(size_t r=0; r<n; r++)
		{
			for(size_t z=0; z<height.size(); z++)
			{
                /*B[r][tau]+=Far[z]\
                    *Pulse(height[r]-height[z], double(Length-tau)*0.15*DELTAT)\
                    *Pulse(height[r]-height[z]+double(tau)*0.15*DELTAT, \
					double(Length-tau)*0.15*DELTAT)*B0[z][tau];*/
                B[r][tau]+=Far[z]\
                    *Pulse(height[r]-height[z], double(Length)*0.15*DELTAT)\
                    *Pulse(height[r]-height[z]+double(tau)*0.15*DELTAT, \
					double(Length)*0.15*DELTAT)*B0[z][tau];
			}
            norm[r]=B[r][0];
		}
	}
    for(size_t tau=0; tau<Length; tau++)
	{
		for(size_t r=0; r<n; r++)
            B[r][tau]/=norm[r];
    }
	// printf("B size is %d:%d; CM size is %d:%d\n", B.size(), B[0].size(), B0.size(), B0[0].size());
    return B;
}

void Get_Tprofiles(vec x, Data_T D, vec &Te_p, vec &Tr_p,\
     std::vector<vec> &B, vec &A)
{
    double lambda=3.0e5/155500.0, Te, Tr;
    std::vector<vecc> CorMatrix;
    vec height, far, knots;
	std::vector<bool> isempty;
	CorMatrix=D.CorMatrix;
	height=D.height;
	isempty=D.isempty;
    far=D.far;
    knots=D.knots;
    unsigned nE=D.nE;
    unsigned nR=D.nR;
    A.resize(height.size());

    vec be, br;
    for(size_t i=0; i<nE; i++)
        be.push_back(x[i]);
    br.push_back(1.0);
    br.push_back(1.0);
    for(size_t i=0; i<nR; i++)
        br.push_back(x[i+nE]);
    br.push_back(1.0);
    br.push_back(1.0);
    Te_p=T_profile(be, knots, height);
    Tr_p=T_profile(br, knots, height);
    std::vector<vec> B0;
    for(size_t t=0; t<height.size(); t++)
    {
        Te=100.0*round(Te_p[t]/100.0);
        Tr=floor(Tr_p[t])+round(10.0*(Tr_p[t]-floor(Tr_p[t])))/10.0;
        vec R=D.Ssec.find({lambda, Te, Tr})->second;
        if(!R.empty())
            B0.push_back(R);
        else
            printf("These temperature values don't exist!\n");
    }
    B=Radeq_ACF(B0, far, height, D.Pulse_Length);
    for(size_t t=0; t<height.size(); t++)
    {
        double num=0.0, denum=0.0, a;
        if(!isempty[t])
        {
            for(size_t r=D.Clutter_nlags; r<D.Pulse_Length; r++)
            {
                num+=B[t][r]*real(CorMatrix[t][r]);
                denum+=B[t][r]*B[t][r];
            }
            a=num/denum;
        }
        else
            a=sqrt(-1.0);
        A[t]=a;
    }
}

double Ti_constraint(const vec &x, vec &grad, void *data)
{
    constraint_data *d=reinterpret_cast<constraint_data*>(data);
    size_t i=d->i;
    size_t j=d->j;
    return x[i]-x[j];
}

double Tr_constraint(const vec &x, vec &grad, void *data)
{
    constraint_data *d=reinterpret_cast<constraint_data*>(data);
    size_t i=d->i;
    size_t nE=d->nE;
    size_t nR=d->nR;
    double res;
    vec bx;
    for(size_t j=0; j<nE; j++)
        bx.push_back(x[j]);
    bx.push_back(1.0);
    bx.push_back(1.0);
    for(size_t j=0; j<nR; j++)
        bx.push_back(x[j+nE]);
    bx.push_back(1.0);
    bx.push_back(1.0);
    if(i==0 or i==nE-2)
    {
        // printf("x[%zu]=%f <= x[%zu]=%f\n", i, x[i], i+1, x[i+1]);
        return x[i]-x[i+1];
    }
    else if(i==1)
    {
        // printf("x[%zu]*x[%zu]=%f <= x[%zu]=%f\n",\
        //  i, nE+i-1, x[i]*x[nE+i-1], i+1, x[i+1]);
        return x[i]*x[nE+i-1]-x[i+1];
    }
    else if(i==nE-3)
    {
        // printf("x[%zu]=%f <= x[%zu]*x[%zu]=%f\n",\
        //  i, x[i], i+1, nE+i-2, x[i+1]*x[nE+i-2]);
        return x[i]-x[i+1]*x[nE+i-2];
    }
    else
    {
        // printf("x[%zu]*x[%zu]=%f <= x[%zu]*x[%zu]=%f\n",\
        //  i, nE+i-1, x[i]*x[nE+i-1], i+1, nE+i-2, x[i+1]*x[nE+i-2]);
        return x[i]*x[nE+i-1]-x[i+1]*x[nE+i-2];
    }
}

double Te_constraint(const vec &x, vec &grad, void *data)
{
    constraint_data *d=reinterpret_cast<constraint_data*>(data);
    size_t i=d->i;
    return x[i-1]-x[i];
}

double ObjFuncT(const vec &x, vec &grad, void *f_data)
{
    std::ofstream out;
    double chisq=0.0;
    double lambda=3.0e5/155500.0, Te, Tr;
    std::vector<vecc> CorMatrix;
    vec height, far, knots;
	std::vector<bool> isempty;
	Data_T* D;
	D=(Data_T*)f_data;
	CorMatrix=D->CorMatrix;
	height=D->height;
	isempty=D->isempty;
    far=D->far;
    knots=D->knots;
    unsigned nE=D->nE;
    unsigned nR=D->nR;
    // printf("%d, %d\n", x.size(), knots.size());

    vec be, br;
    for(size_t i=0; i<nE; i++)
        be.push_back(x[i]);
    br.push_back(1.0);
    br.push_back(1.2);
    for(size_t i=0; i<nR; i++)
        br.push_back(x[i+nE]);
    br.push_back(1.0);
    br.push_back(1.0);
    vec Te_p=T_profile(be, knots, height);
    vec Tr_p=T_profile(br, knots, height);
	// printf("far size is %d; height size is %d\n", far.size(), height.size());

    std::vector<vec> B, B0;
    for(size_t t=0; t<height.size(); t++)
    {
        Te=100.0*round(Te_p[t]/100.0);
        Tr=floor(Tr_p[t])+round(10.0*(Tr_p[t]-floor(Tr_p[t])))/10.0;
        // std::cout<<lambda<<"\t"<<Te<<"\t"<<Tr<<"\n";
        vec R=D->Ssec.find({lambda, Te, Tr})->second;
        if(!R.empty())
            B0.push_back(R);
    }

    B=Radeq_ACF(B0, far, height, D->Pulse_Length);
	// printf("B size is %d:%d; CM size is %d:%d\n", B.size(), B[0].size(), CorMatrix.size(), CorMatrix[0].size());
    for(size_t t=0; t<height.size(); t++)
    {
        double num=0.0, denum=0.0, a;
        if(!isempty[t])
        {
            for(size_t r=D->Clutter_nlags; r<D->Pulse_Length; r++)
            {
                num+=B[t][r]*real(CorMatrix[t][r]);
                denum+=B[t][r]*B[t][r];
            }
            a=num/denum;

            for(size_t r=D->Clutter_nlags; r<D->Pulse_Length; r++)
                chisq+=pow(real(CorMatrix[t][r])-B[t][r], 2.0);
        }
        else
            a=sqrt(-1.0);
    }

    std::string message;
    for(size_t i=0; i<nE; i++)
        message=message+std::to_string(x[i])+" | ";
    for(size_t i=0; i<nR; i++)
        message=message+std::to_string(x[i+nE])+" | ";
    message=message+"\n";

	if(D->Global_solver_on)
	{
		fprintf(stderr, "%s", message.c_str());
	}

    return sqrt(chisq)/double(D->Pulse_Length-D->Clutter_nlags);
}

void Global_solver_T(Data_T data, double UT, \
	Optim_option &option, vec &x, std::string &mesout)
{
	double current_time=omp_get_wtime();
    char message[1000];
    size_t npar=x.size();
	vec lb(npar), ub(npar);
    unsigned nE=data.nE;
    unsigned nR=data.nR;
	unsigned maxeval=option.maxeval;

    for(size_t i=0; i<nE; i++)
    {
        lb[i]=400.0;
        ub[i]=4000.0;
    }
    for(size_t i=0; i<nR; i++)
    {
        lb[i+nE]=1.0;
        ub[i+nE]=3.0;
    }
	nlopt::opt optg(nlopt::GN_ISRES, npar);
    nlopt::opt optl(nlopt::LN_COBYLA, npar);
    optg.nlopt::opt::set_local_optimizer(optl);
	optg.nlopt::opt::set_lower_bounds(lb);
	optg.nlopt::opt::set_upper_bounds(ub);
	optg.nlopt::opt::set_maxeval(maxeval);
	optl.nlopt::opt::set_xtol_rel(1.0e-3);
    // optg.nlopt::opt::set_maxtime(3600);
    optg.set_min_objective(ObjFuncT, &data);

    std::vector<constraint_data> de, dr, di;
    // for(size_t i=1; i<nE; i++)
    //     de.push_back({i, 0, nE, nR});
	// for(size_t i=0; i<de.size(); i++)
    //     optg.nlopt::opt::add_inequality_constraint(Te_constraint, &de[i], 10.0);
	for(size_t i=0; i<nE-1; i++)
    	di.push_back({i, i+1, nE, nR});
	for(size_t i=0; i<di.size(); i++)
        optg.nlopt::opt::add_inequality_constraint(Ti_constraint, &di[i], 1.0e-3);
    for(size_t i=0; i<nE-1; i++)
        dr.push_back({i, 0, nE, nR});
	for(size_t i=0; i<dr.size(); i++)
		optg.nlopt::opt::add_inequality_constraint(Tr_constraint, &dr[i], 1.0e-3);

	try
	{
		nlopt::result result=optg.optimize(x, option.minf);
		// option.neval=optg.nlopt::opt::get_numevals();
		option.neval=1;
		sprintf(message, "#%2.2d | UT=%2.2d:%2.2d || %2.3f | %3.0f | %4.0d GN!\n",\
		 omp_get_thread_num(), int(floor(UT)), int(round((UT-floor(UT))*60.0)),\
		 option.minf, omp_get_wtime()-current_time, option.neval);

		option.time=omp_get_wtime()-current_time;
		option.opt_failed=false;
	}
	catch(std::exception &e)
	{
		sprintf(message, "#%2.2d (UT=%2.2d:%2.2d), Global optimization failed:%s\n",
		 omp_get_thread_num(), int(floor(UT)), \
         int(round((UT-floor(UT))*60.0)), e.what());

		option.opt_failed=true;
		for(size_t i=0; i<npar; i++)
			x[i]=0.0/0.0;
	}
    mesout=mesout+message;
}

void Local_solver_T(Data_T data, double UT, \
	Optim_option &option, vec &x, std::string &mesout)
{
	double current_time=omp_get_wtime();
    char message[1000];
    size_t npar=x.size();
	vec lb(npar), ub(npar);
    unsigned nE=data.nE;
    unsigned nR=data.nR;
	printf("npar=%d\n", npar);
    for(size_t i=0; i<nE; i++)
    {
        lb[i]=400.0;
        ub[i]=4000.0;
    }
    for(size_t i=0; i<nR; i++)
    {
        lb[i+nE]=1.0;
        ub[i+nE]=3.0;
    }
    nlopt::opt optl(nlopt::LN_COBYLA, npar);
	optl.nlopt::opt::set_lower_bounds(lb);
	optl.nlopt::opt::set_upper_bounds(ub);
    optl.nlopt::opt::set_maxeval(option.maxeval);
    optl.nlopt::opt::set_xtol_rel(1.0e-4);
    optl.set_min_objective(ObjFuncT, &data);

    std::vector<constraint_data> de, dr, di;
    // for(size_t i=1; i<nE; i++)
    //     de.push_back({i, 0, nE, nR});
	// for(size_t i=0; i<de.size(); i++)
    //     optl.nlopt::opt::add_inequality_constraint(Te_constraint, &de[i], 1.0);
	for(size_t i=0; i<nE-1; i++)
    	di.push_back({i, i+1, nE, nR});
	for(size_t i=0; i<di.size(); i++)
        optl.nlopt::opt::add_inequality_constraint(Ti_constraint, &di[i], 1.0e-3);
    for(size_t i=0; i<nE-1; i++)
        dr.push_back({i, 0, nE, nR});

	for(size_t i=0; i<dr.size(); i++)
		optl.nlopt::opt::add_inequality_constraint(Tr_constraint, &dr[i], 0.0);

	try
	{
		nlopt::result result=optl.optimize(x, option.minf);
		// option.neval=optg.nlopt::opt::get_numevals();
		option.neval=1;
		sprintf(message, "#%2.2d | UT=%2.2d:%2.2d || %2.3f | %3.0f | %4.0d GN!\n",\
		 omp_get_thread_num(), int(floor(UT)), int(round((UT-floor(UT))*60.0)),\
		 option.minf, omp_get_wtime()-current_time, option.neval);
		option.time=omp_get_wtime()-current_time;
		option.opt_failed=false;
	}
	catch(std::exception &e)
	{
		printf("#%2.2d (UT=%2.2d:%2.2d), Local optimization failed:%s\n",
		 omp_get_thread_num(), int(floor(UT)), \
         int(round((UT-floor(UT))*60.0)), e.what());

		option.opt_failed=true;
		for(size_t i=0; i<npar; i++)
			x[i]=0.0/0.0;
	}
    mesout=mesout+message;
}

double polynom(double x, void *p)
{
	struct polynom_params *params=(struct polynom_params *)p;
	int n=(params->n);
	int trunc=(params->trunc);
	double a[n], norm;
	for(int i=0; i<n; i++)
	   a[i]=(params->A[i]);
	norm=gsl_poly_eval(a, n, 0.0);
	return gsl_poly_eval(a, n, x)/norm;
}

double polynom_deriv(double x, void *p)
{
	struct polynom_params *params=(struct polynom_params *)p;
	int n=(params->n);
	int T=(params->T);
	int trunc=(params->trunc);
	double a[n], b[n-1], norm;
	for(int i=0; i<n; i++)
	   a[i]=(params->A[i]);
	b[0]=2.0*a[2];
	b[1]=3.0*a[3];
	b[2]=4.0*a[4];
	b[3]=5.0*a[5];
	b[4]=6.0*a[6];
	double B=gsl_poly_eval(a, n, x);
	double dB=gsl_poly_eval(b, n-1, x);
	double D=(1.0-double(x)/double(T));
	double dD=-1.0/double(T);
	// double res=B*dD-dB*D;
	double res=dB;
	return res;
}

void polynom_fdf(double x, void *p, double *y, double *dy)
{
	struct polynom_params *params=(struct polynom_params *)p;
	int n=(params->n);
	int trunc=(params->trunc);
	double a[n], b[n-1], c[n-2];
	for(int i=0; i<n; i++)
	   a[i]=(params->A[i]);
	b[0]=2.0*a[2];
	b[1]=3.0*a[3];
	b[2]=4.0*a[4];
	b[3]=5.0*a[5];
	b[4]=6.0*a[6];
	c[0]=3.0*a[3];
	c[1]=8.0*a[4];
	c[2]=15.0*a[5];
	c[3]=24.0*a[6];

	*y=gsl_poly_eval(b, n-1, x);
	*dy=gsl_poly_eval(c, n-2, x);
}

double polyval(double x, vec A)
{
	A.insert(A.begin()+1, 0.0);
	double a[A.size()];
	std::copy(A.begin(), A.end(), a);
	// double res=A[0]+A[1]*pow(x, 2.0)+A[2]*pow(x, 3.0)+A[3]*pow(x, 4.0)+A[4]*pow(x, 5.0)+A[5]*pow(x, 6.0);
	return gsl_poly_eval(a, A.size(), x);
}

void ACF_pardet(double &t_0, double &t_min, double &a_min, double N, int trunc, vec p)
{
	double a[POLY_DEGREE+1];
	double x_lo=double(trunc),  r=0.0;
	p.insert(p.begin()+1, 0.0);
	std::copy(p.begin(), p.end(), a);
	struct polynom_params params={a, int(N), POLY_DEGREE+1, trunc};

	double t=double(trunc);
	double B=polynom(t, &params);
	while(t<N and B>0.0)
	{
		t+=0.1;
		B=polynom(t, &params)/(1.0-t/N);
	}
	t_0=t;
	
	t+=1.0;
	B=polynom_deriv(t, &params)/(1.0-t/N);
	while(t<N and B<0.0)
	{
		t+=0.1;
		B=polynom_deriv(t, &params)/(1.0-t/N);
	}
	t_min=t;

	a_min=-1.0*polynom(t, &params)/(1.0-t/N);

	// printf("t_0 = %f, t_min = %f, a_min = %f\n", t_0, t_min, a_min);
}

void ACFEST_pardet(double &t_0, double &t_min, double &a_min, double N, int trunc, vec p)
{
	double a[POLY_DEGREE+1];
	double x_lo=double(trunc),  r=0.0;
	p.insert(p.begin()+1, 0.0);
	std::copy(p.begin(), p.end(), a);
	struct polynom_params params={a, int(N), POLY_DEGREE+1, trunc};

	double t=double(trunc);
	double B=polynom(t, &params);
	while(t<N and B>0.0)
	{
		t+=0.1;
		B=polynom(t, &params);
	}
	t_0=t;
	
	t+=1.0;
	B=polynom_deriv(t, &params);
	while(t<N and B<0.0)
	{
		t+=0.1;
		B=polynom_deriv(t, &params);
	}
	t_min=t;

	a_min=-1.0*polynom(t, &params);

	// printf("t_0 = %f, t_min = %f, a_min = %f\n", t_0, t_min, a_min);
}

/* void ACF_pardet(double &t_0, double &t_min, double &a_min, double T, int trunc, vec p)
{
	double a[POLY_DEGREE+1];
	double b[POLY_DEGREE-1];
	double z[2*POLY_DEGREE];
	double y[2*POLY_DEGREE];
	double dz[2*POLY_DEGREE-4];
	bool first=true;
	gsl_poly_complex_workspace * w=gsl_poly_complex_workspace_alloc (POLY_DEGREE+1);
	gsl_poly_complex_workspace * dw=gsl_poly_complex_workspace_alloc (POLY_DEGREE-1);
	p.insert(p.begin()+1, 0.0);
	std::copy(p.begin(), p.end(), a);
	struct polynom_params params={a, int(T), POLY_DEGREE+1, trunc};
	// a={p[0], p[1], p[2], p[3], p[4], p[5], p[6]};
	gsl_poly_complex_solve(a, POLY_DEGREE+1, w, z);
	// if(!GSL_SUCCESS)
	// 	printf("Polynom roots failed!\n");

	for(size_t i=0; i<POLY_DEGREE; i++)
	{
		if(first and z[2*i]>0.0 and z[2*i+1]==0.0 and z[2*i]<T)
		{
			t_0=z[2*i];
			first=false;
		}
		// t_0=z[2];
		// first=false;
	}
	if(!first)
	{
		first=true;
		// a[0]=p[0]/T;
		// a[1]=2.0*p[2];
		// a[2]=(3.0*p[3]-p[2]/T);
		// a[3]=(4.0*p[4]-2.0*p[3]/T);
		// a[4]=(5.0*p[5]-4.0*p[4]/T);
		// a[5]=(6.0*p[6]-4.0*p[5]/T);
		// a[6]=-5.0*p[6]/T;
		// gsl_poly_complex_solve (a, POLY_DEGREE+1, w, y);

		b[0]=2.0*p[2];
		b[1]=3.0*p[3];
		b[2]=4.0*p[4];
		b[3]=5.0*p[5];
		b[4]=6.0*p[6];
		gsl_poly_complex_solve (b, POLY_DEGREE-1, dw, dz);

		for(size_t i=0; i<POLY_DEGREE; i++)
		{
			if(first and dz[2*i]>t_0 and dz[2*i+1]==0.0 and dz[2*i]<T)
			{
				t_min=dz[2*i];
				first=false;
			}
			// printf("x(%d)=%f+i*%f;\t", i, dz[2*i]*DELTAT, dz[2*i+1]*DELTAT);
		}
		// printf("\n");

		// size_t t=round(t_0);
		// double f_prev=polynom(double(t-1), &params);//(1.0-double(t-1)/double(T));
		// double f_cur=polynom(double(t), &params);//(1.0-double(t)/double(T));
		// while(t<T and f_cur<f_prev)
		// {
		// 	t++;
		// 	f_prev=polynom(double(t-1), &params);//(1.0-double(t-1)/double(T));
		// 	f_cur=polynom(double(t), &params);//(1.0-double(t)/double(T));
		// }
		
		// t_min=double(t);
		first=false;
		if(t_0<0.0 or t_0>=T)
			t_0=sqrt(-1.0);
		if(t_min<0.0 or t_min>=T)
			t_min=sqrt(-1.0);
		t_0*=double(DELTAT);
		a_min=polynom(t_min, &params);//(1.0-double(t_min)/double(T));
		t_min*=double(DELTAT);
	}
	else
	{
		// printf("ACF t_0 not defined!\n");
		t_0=sqrt(-1.0);
		t_min=sqrt(-1.0);
		a_min=sqrt(-1.0);
	}

	// if(!first)
	// {
		// std::copy(p.begin(), p.end(), a);
		// a_min=gsl_poly_eval(a, POLY_DEGREE+1, t_min);
		// a_min/=a[0];
		// a_min=polynom(double(t_min), &params);//(1.0-double(t_min)/double(T));
		
	// }
	// else
	// {
		// printf("ACF t_min not defined!\n");
		// t_min=sqrt(-1.0);
		// a_min=sqrt(-1.0);
	// }
	if(abs(a_min)>1.0)
		a_min=sqrt(-1.0);	

	gsl_poly_complex_workspace_free (w);
	gsl_poly_complex_workspace_free (dw);
} */

int repoly_f(const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t n=((struct data *)data)->n;
	size_t trunc=((struct data *)data)->trunc;
	double *y=((struct data *)data)->y;
	double A[POLY_DEGREE];
	for(size_t i=0; i<POLY_DEGREE; i++)
	{
		A[i]=gsl_vector_get(x, i);
	}
	for(double t=0; t<n; t+=1.0)
	{
		double Yi=A[0]+A[1]*pow(t, 2.0)+A[2]*pow(t, 3.0)+A[3]*pow(t, 4.0)+A[4]*pow(t, 5.0)+A[5]*pow(t, 6.0);
		gsl_vector_set(f, int(t), Yi-y[int(t)]);
	}

	return GSL_SUCCESS;
}

int repoly_df(const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t n=((struct data *)data)->n;
	size_t trunc=((struct data *)data)->trunc;
	double *y=((struct data *)data)->y;
	for(size_t i=0; i<n; i++)
	{
		double t=double(i);
		gsl_matrix_set(J, i, 0, 1.0);
		gsl_matrix_set(J, i, 1, pow(t, 2.0));
		gsl_matrix_set(J, i, 2, pow(t, 3.0));
		gsl_matrix_set(J, i, 3, pow(t, 4.0));
		gsl_matrix_set(J, i, 4, pow(t, 5.0));
		gsl_matrix_set(J, i, 5, pow(t, 6.0));
	}
	return GSL_SUCCESS;
}

/* void ReACF_interpolator(double *y, size_t n, int trunc, vec &Ar, vec &Er, double &chi)
{
	// printf("!\n");
	const gsl_multifit_fdfsolver_type *T=gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *s;
	int status, info;
	const size_t p=POLY_DEGREE;
	gsl_vector *res_f;
	double chi0;
	s=gsl_multifit_fdfsolver_alloc(T, n, p);

	const double xtol=1e-8;
	const double gtol=1e-8;
	const double ftol=0.0;

	gsl_matrix *J=gsl_matrix_alloc(n, p);
	gsl_matrix *covar=gsl_matrix_alloc(p, p);
	struct data d={n, trunc, y};
	gsl_multifit_function_fdf f;
	double x_init[POLY_DEGREE]={1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	gsl_vector *w=gsl_vector_alloc(n);
	for(double i=0; i<n; i++)
	{
		gsl_vector_set(w, i, pow(1.0-i/double(n), 2.0));
		if(i>n-trunc)
			gsl_vector_set(w, i, 0.5);
	}
	gsl_vector_view x=gsl_vector_view_array(x_init, p);

	f.f=&repoly_f;
	f.df=&repoly_df;
	f.n=n;
	f.p=p;
	f.params=&d;

	gsl_multifit_fdfsolver_wset(s, &f, &x.vector, w);
	// gsl_multifit_fdfsolver_set(s, &f, &x.vector);
	res_f=gsl_multifit_fdfsolver_residual(s);
	chi0=gsl_blas_dnrm2(res_f);
	status=gsl_multifit_fdfsolver_driver(s, 200, xtol, gtol, ftol, &info);
	gsl_multifit_fdfsolver_jac(s, J);
	gsl_multifit_covar(J, 0.0, covar);
	chi=gsl_blas_dnrm2(res_f);
	// gsl_blas_ddot(res_f, res_f, &chi);

	// #define FIT(i) gsl_vector_get(s->x, i)
	// #define ERR(i) sqrt(chi*gsl_matrix_get(covar, i, i))

	// Ar.resize(POLY_DEGREE+1);
	// Er.resize(POLY_DEGREE+1);
	for(size_t i=0; i<p; i++)
	{
		Ar[i]=gsl_vector_get(s->x, i);
		Er[i]=sqrt(chi*gsl_matrix_get(covar,i,i));
	}

	gsl_matrix_free(covar);
	gsl_matrix_free(J);
	gsl_multifit_fdfsolver_free(s);
	// printf("Real ACF interpolation is correct. %g.\n", y[10]);
} */

void ReACF_interpolator(const vec R, const unsigned trunc, vec &Ar, vec &Er, double &chi)
{
	const size_t p=POLY_DEGREE;
	const size_t n=R.size();

	gsl_matrix *X=gsl_matrix_alloc(n, p);
	gsl_vector *y=gsl_vector_alloc(n);
	gsl_vector *w=gsl_vector_alloc(n);

	for(size_t i=trunc; i<n; i++)
	{
		double t=double(i);
		gsl_matrix_set(X, t, 0, 1.0);
		gsl_matrix_set(X, t, 1, pow(double(t), 2.0));
		gsl_matrix_set(X, t, 2, pow(double(t), 3.0));
		gsl_matrix_set(X, t, 3, pow(double(t), 4.0));
		gsl_matrix_set(X, t, 4, pow(double(t), 5.0));
		gsl_matrix_set(X, t, 5, pow(double(t), 6.0));
		if(R[t]>1.1*R[trunc] or R[t]<-0.7*R[trunc])
			gsl_vector_set(y, t, 0.0);
		else
			gsl_vector_set(y, t, R[t]);

		// gsl_vector_set(w, t, 1.0);
		gsl_vector_set(w, t, pow(1.0-double(t)/double(n), 1.0));
		if(1.0-double(t)/double(n)<=0.0)
		{
			gsl_vector_set(y, t, 0.0);
			gsl_vector_set(w, t, 1.0);	
		}
		// gsl_vector_set(w, t, 1.0-double(t)/double(n));
	}

	gsl_multifit_linear_workspace *work=gsl_multifit_linear_alloc(n, p);
	gsl_vector *c=gsl_vector_alloc(p);
	gsl_vector *c_lcurve=gsl_vector_alloc(p);
	gsl_vector *reg_param=gsl_vector_alloc(n);
	gsl_vector *rho=gsl_vector_alloc(n);
	gsl_vector *eta=gsl_vector_alloc(n);
	gsl_vector *G=gsl_vector_alloc(n);
	gsl_matrix *cov=gsl_matrix_alloc(p, p);
	double lambda_l;
	double lambda_gcv;
	double G_gcv;
	size_t reg_idx, rank;
	double rcond;
	double rnorm, snorm;
	chi=0.0;
	gsl_vector_set(c, 0, 1.0);
	for(size_t i=1; i<p; i++)
		gsl_vector_set(c, i, 0.0);

	gsl_multifit_wlinear(X, w, y, c, cov, &chi, work);
	// gsl_multifit_wlinear_tsvd(X, w, y, 0.001, c, cov, &chi, &rank, work);
	
	gsl_vector *L=gsl_vector_alloc(p);
	gsl_matrix *Xs=gsl_matrix_alloc(n, p);
	gsl_vector *ys=gsl_vector_alloc(n);
	gsl_vector *cs=gsl_vector_alloc(p);
	// for(size_t i=0; i<p; i++)
	// 	gsl_vector_set(L, i, 1.0);
	// gsl_multifit_linear_wstdform1(L, X, w, y, Xs, ys, work);
	// gsl_multifit_linear_svd(Xs, work);
	// gsl_multifit_linear_solve(0.0, Xs, ys, cs, &rnorm, &snorm, work);
	// chi=pow(rnorm, 2.0);
	// gsl_multifit_linear_genform1(L, cs, c, work);

	for(size_t i=0; i<p; i++)
	{
		Ar[i]=gsl_vector_get(c, i);
		Er[i]=sqrt(chi*gsl_matrix_get(cov,i,i));
	}
	if(!isfinite(Ar[0]))
	{
		// fprintf(stderr, "ACF interpolation failed by thread %d!\t", omp_get_thread_num());
		for(size_t i=0; i<p; i++)
		{
			Ar[i]=0.0;
			Er[i]=0.0;
		}
	}

	// gsl_multifit_linear_lcurve(ys, reg_param, rho, eta, work);
	// gsl_multifit_linear_lcorner(rho, eta, &reg_idx);
	// lambda_l=gsl_vector_get(reg_param, reg_idx);
	// gsl_multifit_linear_solve(lambda_l, Xs, ys, cs, &rnorm, &snorm, work);
	// gsl_multifit_linear_genform1(L, cs, c, work);
	// for(size_t i=0; i<p; i++)
	// 	Ar[i]=gsl_vector_get(c, i);
	// chi=pow(rnorm, 2.0)+pow(lambda_l*snorm, 2.0);

	// gsl_multifit_linear_gcv(ys, reg_param, G, &lambda_l, &G_gcv, work);
	// gsl_multifit_linear_solve(lambda_l, Xs, ys, cs, &rnorm, &snorm, work);
	// gsl_multifit_linear_genform1(L, cs, c, work);
	// chi=pow(rnorm, 2.0)+pow(lambda_l*snorm, 2.0);
	// for(size_t i=0; i<p; i++)
	// {
	// 	Ar[i]=gsl_vector_get(c, i);
	// 	Er[i]=sqrt(chi*gsl_matrix_get(cov,i,i));
	// }

	gsl_multifit_linear_free(work);
	gsl_vector_free(c);
	gsl_vector_free(cs);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(ys);
	gsl_vector_free(c_lcurve);
	gsl_vector_free(reg_param);
	gsl_vector_free(rho);
	gsl_vector_free(eta);
	gsl_vector_free(G);
	gsl_vector_free(L);

	gsl_matrix_free(X);
	gsl_matrix_free(Xs);
	gsl_matrix_free(cov);
}


/* void ReACF_interpolator(vec y, size_t n, int trunc, vec &Ar, vec &Er, double &chi)
{
	// printf("!\n");
	const gsl_multifit_nlinear_type *T=gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *s;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters();
	fdf_params.trs=gsl_multifit_nlinear_trs_lm;

	int status, info;
	const size_t p=POLY_DEGREE;
	gsl_vector *res_f;

	const double xtol=1e-8;
	const double gtol=1e-8;
	const double ftol=0.0;

	gsl_matrix *J;
	gsl_matrix *covar=gsl_matrix_alloc(p, p);
	struct data d={n, trunc, y};
	// printf("%f %f\n", Ar[0], Ar[6]);
	// double x_init[POLY_DEGREE]={Ar[0], Ar[2], Ar[3], Ar[4], Ar[5], Ar[6]};
	double x_init[POLY_DEGREE]={1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double weights[n];
	for(size_t i=0; i<n; i++)
		weights[i]=double(n)/double(i+1);

	gsl_vector_view x=gsl_vector_view_array(x_init, p);
	gsl_vector_view w=gsl_vector_view_array(weights, n);

	fdf.f=&repoly_f;
	fdf.df=&repoly_df;
	fdf.fvv=NULL;
	fdf.n=n;
	fdf.p=p;
	fdf.params=&d;

	s=gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);
	gsl_multifit_nlinear_winit (&x.vector, &w.vector, &fdf, s);
	res_f=gsl_multifit_nlinear_residual(s);

	status=gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, NULL, NULL, &info, s);
	J=gsl_multifit_nlinear_jac(s);
	gsl_multifit_nlinear_covar(J, 0.0, covar);
	gsl_blas_ddot(res_f, res_f, &chi);

	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))

	// Ar.resize(POLY_DEGREE+1);
	// Er.resize(POLY_DEGREE+1);
	double dof=n-p;
	double c=GSL_MAX_DBL(1, sqrt(chi/dof));
	Ar={FIT(0), FIT(1), FIT(2), FIT(3), FIT(4), FIT(5)};
	Er={c*ERR(0), c*ERR(1), c*ERR(2), c*ERR(3), c*ERR(4), c*ERR(5)};

	// fprintf(stderr, "summary from method '%s/%s'\n",
	// 	gsl_multifit_nlinear_name(s),
	// 	gsl_multifit_nlinear_trs_name(s));
	// fprintf(stderr, "number of iterations: %zu\n",
	// 	gsl_multifit_nlinear_niter(s));
	// fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
	// fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
	// fprintf(stderr, "reason for stopping: %s\n",
	// 	(info == 1) ? "small step size" : "small gradient");

	gsl_matrix_free(covar);
	gsl_multifit_nlinear_free(s);
	// printf("Real ACF interpolation is correct. %f %f.\n", sqrt(chi/dof), ERR(0));
} */

void Get_Temperatures(double &Te, double &Tr, double t_0, double t_min, double A_min)
{
	double beta_e[5], beta_r[3];
	std::ifstream inp;
	inp.open("beta_e.dat");
	for(int i=0; i<5; i++)
		inp>>beta_e[i];
	inp.close();
	inp.open("beta_r.dat");
	for(int i=0; i<3; i++)
		inp>>beta_r[i];
	inp.close();
	Te=beta_e[0]/(t_0*t_0)+beta_e[1]/t_0+beta_e[2]/(t_min*t_min)+beta_e[3]/t_min+beta_e[4];
	Tr=A_min*A_min*beta_r[0]+A_min*beta_r[1]+beta_r[2];
	// Ti=Te/Tr;
}

std::vector<string> LineSplit(std::string line, std::string sep)
{
	std::vector<string> buf;
	std::size_t found=0;
	std::size_t it=0;
	while(it!=std::string::npos)
	{
		found=line.find_first_of(sep, found+1);
		if(it<1)
			buf.push_back(line.substr(it, found-it).c_str());
		else
			buf.push_back(line.substr(it+1, found-it-1).c_str());
		it=found;
	}
	if(buf.size()==0)
		fprintf(stderr, "Empty file data!\n");
	// cout<<"buf size "<<buf.size()<<"\n";
	// for(size_t i=0; i<buf.size(); i++)
	// 	printf("%s\t", buf[i].c_str());
	// printf("\n");
	return buf;
}

double F_spline_fit(vector<double> y, double x)
{
	int n=y.size();
	double Y[n], X[n], val=0.0;
	for(int i=0; i<n; i++)
	{
		Y[i]=y[i];
		X[i]=double(i);
	}
	gsl_interp_accel *acc=gsl_interp_accel_alloc();
    gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(spline, X, Y, n);
	val=gsl_spline_eval(spline, x, acc);
	gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
	return val;
}

void B_spline_fit(vector<double> &Cv, vector<double> &Cerr, vector<double> Xv,
	vector<double> Yv, vector<double> Wv)
{
	vector<double> beta;
	const size_t n=Xv.size();
	// printf("vector size %zu:%zu!\n", n, Yv.size());
	// if(Yv.size()!=n or Wv.size()!=n)
		// fprintf(stderr, "Vector size is wrong! %d\t%d\t%d\n", Xv.size(), Yv.size(), Wv.size());
	const size_t ncoeffs=NCOEFFS;
	const size_t nbreak=NBREAK;
	beta.resize(ncoeffs);
	Cv.resize(Xv.size());
	Cerr.resize(Xv.size());
	size_t i, j;
	if(n!=0)
	{
		gsl_bspline_workspace *bw;
		gsl_vector *B;
		gsl_vector *c, *w;
		gsl_vector *x, *y;
		gsl_matrix *X, *cov;
		gsl_multifit_linear_workspace *mw=gsl_multifit_linear_alloc(n, ncoeffs);
		double chisq, Rsq, dof, tss;
		/* allocate a cubic bspline workspace (k=4) */
		bw=gsl_bspline_alloc(4, nbreak);
		B=gsl_vector_alloc(ncoeffs);
		x=gsl_vector_alloc(n);
		y=gsl_vector_alloc(n);
		X=gsl_matrix_alloc(n, ncoeffs);
		c=gsl_vector_alloc(ncoeffs);
		w=gsl_vector_alloc(n);
		cov=gsl_matrix_alloc(ncoeffs, ncoeffs);
		/* this is the data to be fitted */
		for (i=0; i<n; ++i)
		{
			gsl_vector_set(x, i, Xv[i]);
			gsl_vector_set(y, i, Yv[i]);
			gsl_vector_set(w, i, Wv[i]);
		}
		/* use uniform breakpoints on [0, 15] */
		gsl_bspline_knots_uniform(Xv[0], Xv[n-1], bw);
		/* construct the fit matrix X */
		for (i=0; i<n; ++i)
		{
			double xi=gsl_vector_get(x, i);
			/* compute B_j(xi) for all j */
			gsl_bspline_eval(xi, B, bw);
			/* fill in row i of X */
			for (j=0; j<ncoeffs; ++j)
			{
				double Bj=gsl_vector_get(B, j);
				gsl_matrix_set(X, i, j, Bj);
			}
		}
		/* do the fit */
		gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);
		dof=n-ncoeffs;
		tss=gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
		Rsq=1.0-chisq/tss;
		// fprintf(stderr, "chisq/dof=%e, Rsq=%f\n", chisq/dof, Rsq);
		/* output the smoothed curve */
		{
			double xi, yi, yerr;
			for (i=0; i<Xv.size(); i++)
			{
				xi=Xv[i];
				gsl_bspline_eval(xi, B, bw);
				gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
				Cv[i]=yi;
				Cerr[i]=yi*yerr;
			}
		}

		for(j=0; j<ncoeffs; j++)
			beta[j]=gsl_vector_get(c, j);
		gsl_bspline_free(bw);
		gsl_vector_free(B);
		gsl_vector_free(x);
		gsl_vector_free(y);
		gsl_matrix_free(X);
		gsl_vector_free(c);
		gsl_vector_free(w);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
	}
	else
	{
		for(i=0; i<Xv.size(); i++)
		{
			Cv[i]=0.0/0.0;
			Cerr[i]=0.0/0.0;
		}
	}
}

void SkyNoise(double Freq, int DoY, map<int, double> &N_sky)
{
	ifstream inp;
	char buf[100];
	string dir="/home/tashlykov/data-out/sky_noise/", head;
	sprintf(buf, "%ssky_noise_20000101_%3.1fMHz.txt", dir.c_str(), Freq/1000.0);
	string FileName=buf;
	vector<string> vbuf;
	printf("+++++%s+++++\n", FileName.c_str());
	inp.open(FileName);
	if(!inp.is_open())
		fprintf(stderr, "File not found!\n");
	getline(inp, head);
	int hour, minute, it=0;
	double N, ut, ut0;
	while(!inp.eof())
	{
		getline(inp, head);
		vbuf=LineSplit(head, "\t");
		if(vbuf.size()>1)
		{
			// hour=atoi(vbuf[1].substr(0,3).c_str());
			// minute=atoi(vbuf[1].substr(4,2).c_str());
			N=atof(vbuf[3].c_str());
			ut0=atof(vbuf[2].c_str());
			ut=ut0+24.0*double(DoY)/365.2425;
			if(ut>=24.0)
				ut-=24.0;
			it=floor(ut*12.0);
			N_sky[it]=N;
			// printf("%s\t\t%d:%d %f\n", head.c_str(), hour, minute, N);
		}
		vbuf.clear();
	}
	inp.close();
}
