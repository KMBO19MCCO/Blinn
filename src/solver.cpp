#include <cmath>
#include <BlinnSolver/solver.hpp>

#include <Eigen/Dense>
#include <excerpt/excerpt.h>


//---Aliases---
template<typename T>
using Vector2T = Eigen::Matrix<T, 2, 1>;

template<typename T>
using Vector4T = Eigen::Matrix<T, 4, 1>;

template<typename T>
using Matrix2T = Eigen::Matrix<T, 2, 2>;

template<typename T>
using Matrix4T = Eigen::Matrix<T, 4, 4>;
//-------------

//compare floating-point numbers
template<typename T>
bool Solver<T>::is_equal(T a, T b)
{
    return std::fabs(a - b) <= ( ( std::fabs(a) > std::fabs(b) ? std::fabs(b) : std::fabs(a) ) * std::numeric_limits<T>::epsilon());
}

//---CONSTRUCTOR--------
//I PART: normalize coefs! B and C divide by 3, all coefs (include B and C) divide by A.
template<typename T>
void Solver<T>::normalize_coefs()
{
    m_B /= static_cast<T>(3.0);
    m_C /= static_cast<T>(3.0);

    if(is_equal(m_A, static_cast<T>(0.0)))
    {
        m_B /= m_A;
        m_C /= m_A;
        m_D /= m_A;
        m_A = static_cast<T>(1.0);
    }

    m_coefs << m_A, m_B, m_C, m_D;
}

//II PART of constructor
template<typename T>
void Solver<T>::calc_params()
{
        m_t = static_cast<T>(1.0);
        m_t2 = std::pow(m_t, 2);

        m_u = static_cast<T>(0.0);
        m_u2 = std::pow(m_u, 2);
        
        auto tmp = [&](T X, T Y, T Z) 
        {return (pr_product_difference<T>(m_t2, X, -m_u2, Z) + (2 * m_t * m_u * Y));};
        
        m_s = -tmp(m_B, m_C, m_D);
        m_v = tmp(m_A, m_B, m_C);
        
        m_s2 = std::pow(m_s, 2);
        m_v2 = std::pow(m_v, 2);
}

//III PART of constructor
template<typename T>
void Solver<T>::first_covariant()
{
    m_T_big(0, 0) = std::pow(m_t, 3);
    m_T_big(0, 1) = static_cast<T>(3.0) * m_u * m_t2;
    m_T_big(0, 2) = static_cast<T>(3.0) * m_t * m_u2;
    m_T_big(0, 3) = std::pow(m_u, 3);
     
    m_T_big(1, 0) = m_t2 * m_s;
    m_T_big(1, 1) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_u * m_s, -m_t2, m_v);
    m_T_big(1, 2) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_u * m_v, -m_u2, m_s);
    m_T_big(1, 3) = m_u2 * m_v;
    
    m_T_big(2, 0) = m_t2 * m_s;
    m_T_big(2, 1) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_s * m_v, -m_s2, m_u);
    m_T_big(2, 2) = pr_product_difference<T>(static_cast<T>(2.0) * m_u, m_s * m_v, -m_v2, m_t);
    m_T_big(2, 3) = m_u2 * m_v;

    m_T_big(3, 0) = std::pow(m_s, 3);
    m_T_big(3, 1) = static_cast<T>(3.0) * m_v * m_s2;
    m_T_big(3, 2) = static_cast<T>(3.0) * m_s * m_v2;
    m_T_big(3, 3) = std::pow(m_v, 3);
    
    m_coefs_tilda = m_T_big * m_coefs;
}

//IV part of constructor
template<typename T>
void Solver<T>::second_covariant()
{
    //T coef = std::pow(pr_product_difference<T>(m_t, m_v, m_s, m_u), 2);
    //Calc delta1, delta2 and delta3 
    m_delta1 = pr_product_difference<T>(m_A, m_C, m_B, m_B);
    m_delta2 = pr_product_difference<T>(m_A, m_D, m_B, m_C);
    m_delta3 = pr_product_difference<T>(m_B, m_D, m_C, m_C);

    //---Init matrix H - Hessian---
    Matrix2T<T> H;
    H << static_cast<T>(2.0) * m_delta1, m_delta2, m_delta2, static_cast<T>(2.0) * m_delta3;
    m_detH   = pr_product_difference<T>(4 * m_delta1, m_delta3, m_delta2, m_delta2);

    m_T_lit << m_t, m_u, m_s, m_v;

    //Calc delta tilda
    m_H_tilda = m_T_lit * H * m_T_lit.transpose(); // maybe add * coef
    
    /*
    m_delta1_tilda = mat_delta_tilda(0,0) / static_cast<T>(2.0);
    m_delta2_tilda = mat_delta_tilda(0,1);
    m_delta3_tilda = mat_delta_tilda(1,1) / static_cast<T>(2.0);
    m_detH_tilda = coef * m_detH;
    */
}

 
//V part of constructor
template<typename T>
void Solver<T>::third_covariant()
{
    m_Aj = pr_product_difference<T>(std::pow(m_A, 2), m_D, static_cast<T>(3.0)* m_A, m_B * m_C) 
               + static_cast<T>(2.0) * std::pow(m_B, 3);   
    m_Bj = pr_product_difference<T>(std::pow(m_B, 2), m_C, static_cast<T>(2.0), m_A * std::pow(m_C, 2)) + (m_A * m_B * m_D);
    m_Cj = pr_product_difference<T>(std::pow(m_B, 2), static_cast<T>(2.0) * m_D, m_B, std::pow(m_C, 2)) - (m_A * m_C * m_D);
    m_Dj = pr_product_difference<T>(static_cast<T>(3.0) * m_B, m_C * m_D, static_cast<T>(2.0), std::pow(m_C, 3)) - (m_A * std::pow(m_D, 2));

    Vector4T<T> J;
    J << m_Aj, m_Bj, m_Cj, m_Dj;
    m_TJ = m_T_big * J;
}

//Finally, constructor
template<typename T>
Solver<T>::Solver(const std::vector<T>& coefs)
{
        m_A = coefs[3];
        m_B = coefs[2];
        m_C = coefs[1];
        m_D = coefs[0];
        
        normalize_coefs();
        calc_params();

        first_covariant();
        second_covariant();
        third_covariant();
}
//-----------------------

//---CASES OF SOLVING---
template<typename T>
int Solver<T>::solve(std::vector <T>& roots)
{
    T A_line = m_coefs_tilda[0];
    //T B_line = static_cast<T>(0.0);
    //T C_line = (m_t2 * m_H_tilda(0,0) + m_t * m_u * m_H_tilda(0,1) + m_u2 * m_H_tilda(1,1) ) / 2.0;
    T C_line = ( pr_product_difference<T>(m_t2, m_H_tilda(0,0), -m_u2, m_H_tilda(1,1)) + (m_t * m_u * m_H_tilda(0,1)) ) / static_cast<T>(2.0);
    T D_line = m_TJ[0];
    
    if(m_detH < static_cast<T>(0.0))
    {
        T r = det_less_zero(A_line, C_line, D_line);
        roots.push_back(r);
        return 1;
    }
    else
    {
        auto L = L_root();
        auto S = S_root();
        auto M = M_root(L.first, L.second, S.first, S.second);
        
        roots.push_back(L.first / L.second);
        roots.push_back(S.first / S.second);
        roots.push_back(M.first / M.second);
        return 3;
    }
}


template<typename T>
T Solver<T>::det_less_zero(T At, T Cbar, T Dbar)
{
    //T At = (std::pow(m_B, 3) * m_D >= std::pow(m_C, 3) * m_A)? m_A : m_D;

    T T0 = -std::copysign(std::fabs(At) * std::sqrt(-m_detH), Dbar);
    T T1 = -Dbar + T0;
    T p = std::cbrt(T1 / static_cast<T>(2.0));
    T q = is_equal(T0, T1) ? (-p) : (-Cbar / p);
    
    T xt1 = (Cbar <= static_cast<T>(0.0)) ? (p + q) : (-Dbar / (fma(p,p, fma(q,q, Cbar))));
    
    Vector2T<T> solution_tmp;
    solution_tmp << xt1, static_cast<T>(1.0);
    auto solution = solution_tmp.transpose() * m_T_lit;
    
    return solution(0,0) / solution(0,1);
}

template<typename T>
std::pair<T,T> Solver<T>::L_root()
{
    //CbarA: the minus sign that occurs in formulas for the root L in a variable
    T CbarA = -m_delta1;
    T DbarA = pr_product_difference<T>(m_A, m_delta2, static_cast<T>(2.0) * m_B, m_delta1);
    T thetaA = std::fabs(std::atan2(m_A * std::sqrt(m_detH), -DbarA)) / static_cast<T>(3.0);
    T xt1A = 2 * std::sqrt(CbarA) * cos(thetaA);
    //xt3A = 2 * sqrt(-CbarA) * (cos(thetaA) / -2. - (sqrt(3.) / 2.)*sin(thetaA));
    T xt3A =  pr_product_difference<T>(-std::sqrt(CbarA * static_cast<T>(3.0)), sin(thetaA), std::sqrt(CbarA), cos(thetaA));
    T xtL = (xt1A + xt3A > static_cast<T>(2.0) * m_B)? xt1A : xt3A;

    return std::make_pair(xtL - m_B, m_A);
}

template<typename T>
std::pair<T,T> Solver<T>::S_root()
{
    T CbarD = m_delta3;
    T DbarD = pr_product_difference<T>(static_cast<T>(2.0) * m_C, m_delta3, m_D, m_delta2);
    T thetaD = std::fabs(std::atan2(m_D * std::sqrt(m_detH), -DbarD)) / static_cast<T>(3.0);
    T xt1D = static_cast<T>(2.0) * std::sqrt(-CbarD) * cos(thetaD); 
    //T xt3D = 2 * sqrt(-CbarD) * (cos(thetaD) / -2.0 - (sqrt(3.0) / 2.0) *sin(thetaD));
    T xt3D = pr_product_difference<T>(-std::sqrt(-CbarD), cos(thetaD), std::sqrt(-CbarD * 3.0), sin(thetaD)); 
    T xtS = (xt1D + xt3D < static_cast<T>(2.0) * m_C)? xt1D : xt3D;

    return std::make_pair(-m_D, xtS + m_C);

}

template<typename T>
std::pair<T,T> Solver<T>::M_root(T root00, T root01, T root10, T root11)
{
        T E = root01 * root11;
        T F = pr_product_difference<T>(-root00, root11, root01, root10);
        T G = root00 * root10;
        
        return std::make_pair(pr_product_difference<T>(m_C, F, m_B, G), 
                              pr_product_difference<T>(m_C, E, m_B, F));
}

//----------------------------------------------


//Technical Chocolate
template class Solver<long double>;
template class Solver<double>;
template class Solver<float>;
