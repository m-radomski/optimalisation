#include <string.h>

struct Consts
{
	double DA;
	double DB;
	double FR;
	double TA;
	double TF;
};
typedef struct Consts Consts;

double cm2_to_m2(double cm2)
{
	return cm2 / 1e4;
}

double k_constant(double hole_diamater_m2)
{
    double a = 0.98;
    double b = 0.63;
    double g = 9.81;
    double pp = 1;

    return -a*b*hole_diamater_m2*sqrt((2*g)/pp);
}

double dVdt(double V, double t, double k)
{
	double result = 0.0;
	if(V > 0.0) 
	{
		result = k * sqrt(V);
	}

	return result;
}

double fVA(double VA, double VB, double TB, double t, Consts *c) 
{
	double result = dVdt(VA, t, k_constant(c->DA));
    return result;
}

double Fin(double VA, double VB, double TB, double t, Consts *c) 
{
	double result = fabs(fVA(VA, VB, TB, t, c)) + c->FR;
    return result;
}

double Tin(double VA, double VB, double TB, double t, Consts *c) 
{
	double result = ((c->TA*fabs(fVA(VA, VB, TB, t, c))) + (c->TF*c->FR)) / Fin(VA, VB, TB, t, c);
    return result;
}

double fVB(double VA, double VB, double TB, double t, Consts *c) 
{
	double result = dVdt(VB, t, k_constant(c->DB)) + Fin(VA, VB, TB, t, c);
    return result;
}

double fTB(double VA, double VB, double TB, double t, Consts *c) 
{
	double result = (Fin(VA, VB, TB, t, c)/VB) * (Tin(VA, VB, TB, t, c) - TB);
    return result;
}

double sym_get_max_temp(double diamaterAcm2)
{
    double h = 1.0;
    int tfinal = 1000;

    double *VA = (double *)malloc(tfinal * sizeof(double));
    double *VB = (double *)malloc(tfinal * sizeof(double));
    double *TB = (double *)malloc(tfinal * sizeof(double));
    double *t  = (double *)malloc(tfinal * sizeof(double));
	
	memset(VA, 0, tfinal * sizeof(double));
	memset(VB, 0, tfinal * sizeof(double));
	memset(TB, 0, tfinal * sizeof(double));
	memset(t,  0, tfinal * sizeof(double));

	Consts c_ = { 0 };
	Consts *c = &c_;
	c->DA = cm2_to_m2(diamaterAcm2);
	c->DB = cm2_to_m2(36.5);
	c->FR = 0.01;
	c->TA = 90.0;
	c->TF = 10.0;

    VA[0] = 5.0;
    VB[0] = 1.0;
    TB[0] = 10.0;
    t[0] = 0.0;

	double Tmax = 0.0;

	for(int i = 0; i < tfinal - 1; i++)
	{
        t[i+1] = t[i] + h;

        double k1VA = fVA(VA[i], VB[i],TB[i],t[i], c);
        double k1VB = fVB(VA[i], VB[i],TB[i],t[i], c);
        double k1TB = fTB(VA[i], VB[i],TB[i],t[i], c);

        double k2VA = fVA(VA[i]+(h/2)*k1VA, VB[i]+(h/2)*k1VB, TB[i]+(h/2)*k1TB, t[i]+h/2, c);
        double k2VB = fVB(VA[i]+(h/2)*k1VA, VB[i]+(h/2)*k1VB, TB[i]+(h/2)*k1TB, t[i]+h/2, c);
        double k2TB = fTB(VA[i]+(h/2)*k1VA, VB[i]+(h/2)*k1VB, TB[i]+(h/2)*k1TB, t[i]+h/2, c);

        double k3VA = fVA(VA[i]+(h/2)*k2VA, VB[i]+(h/2)*k2VB, TB[i]+(h/2)*k2TB, t[i]+h/2, c);
        double k3VB = fVB(VA[i]+(h/2)*k2VA, VB[i]+(h/2)*k2VB, TB[i]+(h/2)*k2TB, t[i]+h/2, c);
        double k3TB = fTB(VA[i]+(h/2)*k2VA, VB[i]+(h/2)*k2VB, TB[i]+(h/2)*k2TB, t[i]+h/2, c);

        double k4VA = fVA(VA[i]+(h  )*k3VA, VB[i]+(h  )*k3VB, TB[i]+(h  )*k3TB, t[i]+h, c);
        double k4VB = fVB(VA[i]+(h  )*k3VA, VB[i]+(h  )*k3VB, TB[i]+(h  )*k3TB, t[i]+h, c);
        double k4TB = fTB(VA[i]+(h  )*k3VA, VB[i]+(h  )*k3VB, TB[i]+(h  )*k3TB, t[i]+h, c);

        VA[i+1] = VA[i] + (h/6)*(k1VA + 2*k2VA + 2*k3VA + k4VA);
        VB[i+1] = VB[i] + (h/6)*(k1VB + 2*k2VB + 2*k3VB + k4VB);
        TB[i+1] = TB[i] + (h/6)*(k1TB + 2*k2TB + 2*k3TB + k4TB);

		Tmax = Tmax < TB[i] ? TB[i] : Tmax;
	}

	return Tmax;
}
