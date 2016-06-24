#include "class"
#include "math.h"

Class::Class()
{
    toXYZ(40*M_PI/180, -74*M_PI/180, 0);
    x1 = x; y1 = y; z1 = z;
    cout <<x1<<" "<< y1<<" "<< z1;

    toXYZ(55 *M_PI/180, 37*M_PI/180, 0);
    x2 = x; y2 = y; z2 = z;
    calculation_ncos_V(20);
}

void Class::calculation_ncos_V (double beta){
//Расчет направляющих косинусов вектора скорости в момент старта
//beta - угол наклона вектора скорости к плоскости горизонта, рад (подобран в ручную)
//Vxs, Vys, Vzs - направляющие косинусы вектора скорости

    double s1, s2, s0; //Вспомогательные переменные
    double N1, N12, N2; //Нормы, скалярное произведение

    beta = 20*M_PI/180;
    double dx, dy, dz;
    N1 = pow(x1,2) + pow(y1,2) + pow(z1,2);
    N12 = (x1 * x2 + y1 * y2 + z1 * z2) / N1;
    dx = x2 - x1 * N12;
    dy = y2 - y1 * N12;
    dz = z2 - z1 * N12;
    N2 = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
    N1 = sqrt(N1);

    s0 = sqrt(pow(tan(beta),2) + pow((1 - tan(beta)),2));
    s1 = tan(beta) / s0 / N1;
    s2 = (1 - tan(beta)) / s0 / N2;

    //Направляющие косинусы вектора скорости
    Vxs = s1 * x1 + s2 * dx;
    Vys = s1 * y1 + s2 * dy;
    Vzs = s1 * z1 + s2 * dz;
}

double Class::atan2(double s, double c){
    double u;
    if (c == 0){
        if (s == 1)
            return M_PI/2;
        else
            return 3*M_PI/2;
    }
    u = atan(s/c);
    if ((c>0) && (s>=0))
        return u;
    if ((c<0) && (s<=0))
        return (u+M_PI);
    if ((c>0) && (s<=0))
        return (2*M_PI + u);
    else
        return (u+M_PI);
}

void Class::toHAzL(double xx, double yy, double zz, double xp, double yp, double zp){ //Перевод Гринвичских x,y,z в высоту-азимут-дальность

    double c, f, rp; //вспомогательные
    double M[3][3];
    double xnez, ynez, znez;
    double dx, dy, dz;
    double Rxy;
    c = 1 / (1 - a);
    f = sqrt(pow(xp,2) + pow(yp,2));
    rp = sqrt(pow(xp,2) + pow(yp,2) + pow((pow(c,2)*zp), 2));

    M[0][0] = -pow(c,2) * xp * zp / f*rp;
    M[0][1] = -pow(c,2) * yp * zp / f*rp;
    M[0][2] = f / rp;
    M[1][0] = -yp / f;
    M[1][1] = xp / f;
    M[1][2] = 0;
    M[2][0] = xp / rp;
    M[2][1] = yp / rp;
    M[2][2] = pow(c,2) * zp / rp;

    dx = xx - xp;
    dy = yy - yp;
    dz = zz - zp;

    xnez = M[0][0] * dx + M[0][1] * dy + M[0][2] * dz;
    ynez = M[1][0] * dx + M[1][1] * dy + M[1][2] * dz;
    znez = M[2][0] * dx + M[2][1] * dy + M[2][2] * dz;

    L = sqrt(pow(xnez,2) + pow(ynez,2) + pow(znez,2));
    Rxy = sqrt(pow(xnez,2) + pow(ynez,2));
    if (Rxy!=0)
    {
        hz = atan(znez/Rxy);
        Az = atan2(ynez/Rxy, xnez/Rxy);
    }
    else{
        hz = 0;
        Az = 0;
    }
}

void Class::toXYZ(double fi, double la, double h){
    double Rxy, CF, SF, rr, aa;
    CF = cos(fi);
    aa = pow((1-a), 2);
    SF = sin(fi);
    rr = sqrt(pow(CF,2) + aa * pow(SF,2));
    Rxy = Re * CF/rr + h * CF;
    x = Rxy * cos(la);
    y = Rxy * sin(la);
    z = (Re * aa / rr + h) * SF;
}

void Class::toLatLongH(double xx, double yy, double zz){ //Перевод Гринвичских x,y,z в широту-долготу-высоту
    double Rxy, R, Fi0;
    double R2;
    Rxy = sqrt(pow(xx,2) + pow(yy,2));
    La = atan2(yy/Rxy, xx/Rxy);
    R2 = pow(xx,2) + pow(yy,2) + pow(zz,2);
    R = sqrt(R2);
    h = R - Re * (1 - a*pow(zz,2)/R2);
    Fi0 = atan(zz/Rxy);
    Fi = Fi0 + 2*a * (1 - h/R) * zz * Rxy/R2;
}
