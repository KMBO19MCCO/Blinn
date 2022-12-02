#include <iostream>
#include <vector>
#include <cmath> 
#include <excerpt/excerpt.h>


template <typename T>
bool is_equal(T a, T b)
{
    return fabs(a - b) <= ( ( fabs(a) > fabs(b) ? fabs(b) : fabs(a) ) * std::numeric_limits<T>::epsilon());
}

template <typename T>
std::vector<std::vector<T>> blinn(T A, T B, T C, T D)
{
    std::vector<std::vector<T>> roots;
    B /= A;
    C /= A;
    D /= A;
    A  = 1;
    
    B /= 3;
    C /= 3;
    
    T d1 = pr_product_difference<T>(A, C, B, B);
    T d2 = pr_product_difference<T>(A, D, B, C);
    T d3 = pr_product_difference<T>(B, D, C, C);
    T det = pr_product_difference<T>(4.0 * d1, d3, d2, d2);
    
    if (det <= 0)
    {
        T At, Cbar, Dbar;
        if (pow(B, 3) * D >= pow(C, 3) * A)
        {
            At = A;
            Cbar = d1;
            Dbar = pr_product_difference<T>(A, d2, 2 * B, d1);
        }
        else
        {
            At = D;
            Cbar = d3;
            Dbar = pr_product_difference<T>(2 * C, d3, D, d2);
        }
        
        //copysign(a, b) = sign(b) * a (from <cmath>)
        //cbrt(x) like pow(x, 1.0/3.0), but more accurate
        T T0 = -copysign(fabs(At) * sqrt(-det), Dbar);
        T T1 = -Dbar + T0;
        T p = cbrt(T1 / 2.0);
        T q = is_equal<T>(T0, T1) ? (-p) : (-Cbar / p);
        T xt1 = (Cbar <= 0) ? (p + q) : (-Dbar / (fma(p,p, fma(q,q, Cbar))));
        
        if (pow(B, 3) * D >= pow(C, 3) * A)
        {
            roots.push_back({xt1 - B, A});
        }
        else
        {
            roots.push_back({-D, xt1 + C});
        }
                                            
    }
    else
    {
        //---L---
        T CbarA = d1;
        T DbarA = pr_product_difference<T>(A, d2, 2 * B, d1);
        T thetaA = fabs(atan2(A *sqrt(det), -DbarA)) / 3.0;
        T xt1A = 2 * sqrt(-CbarA) * cos(thetaA);
        //xt3A = 2 * sqrt(-CbarA) * (cos(thetaA) / -2. - (sqrt(3.) / 2.)*sin(thetaA));
        T xt3A =  pr_product_difference<T>(-sqrt(-CbarA * 3.0), sin(thetaA), sqrt(-CbarA), cos(thetaA));
        T xtL = (xt1A + xt3A > 2 * B)? xt1A : xt3A;
        
        roots.push_back({xtL - B, A});
        
        //---S---
        T CbarD = d3;
        T DbarD = pr_product_difference<T>(2 * C, d3, D, d2);
        T thetaD = fabs(atan2(D * sqrt(det), -DbarD)) / 3.0;
        T xt1D = 2 * sqrt(-CbarD) * cos(thetaD);
        //T xt3D = 2 * sqrt(-CbarD) * (cos(thetaD) / -2.0 - (sqrt(3.0) / 2.0) *sin(thetaD));
        T xt3D = pr_product_difference<T>(-sqrt(-CbarD), cos(thetaD), sqrt(-CbarD * 3.0), sin(thetaD)); 
                     
        T xtS = (xt1D + xt3D < 2 * C)? xt1D : xt3D;
        roots.push_back({ -D, xtS + C });
                     
        //---M---
        T E = roots[0][1] * roots[1][1];
        T F = pr_product_difference<T>(-roots[0][0], roots[1][1], roots[0][1], roots[1][0]);
        T G = roots[0][0] * roots[1][0];

        roots.push_back({pr_product_difference<T>(C, F, B, G),  pr_product_difference<T>(C, E, B, F)});

    }
    return roots;
}

int main()
{
    #ifdef FP_FAST_FMA
      std::cout << "FMA enabled\n";
    #else
      std::cout << "FMA disabled\n";
    #endif
    
    
    unsigned int p = 3;
    std::vector<double> roots(p);
    std::vector<double> coefs(p + 1);
    
    int code = generate_polynomial<double>(p, 0, 2, 0, 1e-5, 0.0, 1.0, roots, coefs);
    if (code < 0)
    {
        std::cout << "Error in generate_polynomial\n";
    }
    
    std::vector<std::vector<double>> my_roots = blinn(coefs[3],coefs[2],coefs[1], coefs[0]);
    std::vector<double> norm_roots;
    for (auto v : my_roots)
    {
        norm_roots.push_back(v[0] / v[1]);
        
    }

    double max_abs = 0.0;
    double max_rel = 0.0;
    int rv = compare_roots<double>(3, 3, norm_roots, roots, max_abs, max_rel);
    
    //---Output---
    std::cout << "max_abs: " << max_abs << '\n';
    std::cout << "max_rel: " << max_rel << '\n';
    std::cout << "rv: " << rv << '\n';

    std::cout << "Roots: ";
    for(auto elem : roots)
    {
        std::cout << elem <<  " ";
    }
    std::cout << "\n";
    
    std::cout << "My_roots: ";
    for(auto elem : norm_roots)
    {
        std::cout << elem <<  " ";   
    }
    std::cout << "\n";
    
}
