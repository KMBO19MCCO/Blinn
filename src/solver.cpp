#include <cmath>
#include <complex>
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
    //We consider the initial equation of the form: Ax^3 + 3Bx^2 + 3Cx + D = 0
    //For the transition, divide B and C into 3.
    m_B /= static_cast<T>(3.0);
    m_C /= static_cast<T>(3.0);

    //normalization
    if(is_equal(m_A, static_cast<T>(0.0)))
    {
        m_B /= m_A;
        m_C /= m_A;
        m_D /= m_A;
        m_A = static_cast<T>(1.0);
    }
    
    //square of coefs
    m_A2 = m_A * m_A;
    m_B2 = m_B * m_B;
    m_C2 = m_C * m_C;
    m_D2 = m_D * m_D;
    
    
    m_coefs << m_A, m_B, m_C, m_D;
}

//II PART of constructor
template<typename T>
void Solver<T>::calc_params()
{
        //In this implementation, we consider cases when t = 1, u = 0
        m_t = static_cast<T>(1.0);
        m_t2 = m_t * m_t;

        m_u = static_cast<T>(0.0);
        m_u2 = m_u * m_u;
        
        auto tmp = [&](T X, T Y, T Z) 
        {return (pr_product_difference<T>(m_t2, X, -m_u2, Z) + (2 * m_t * m_u * Y));};

        //Calculate s and v 
        m_s = -tmp(m_B, m_C, m_D);
        m_v = tmp(m_A, m_B, m_C);
        
        m_s2 = m_s * m_s;
        m_v2 = m_v * m_v;
}

//III PART of constructor
template<typename T>
void Solver<T>::first_covariant()
{
    //m_T_big is transition matrix for coefs A,B,C,D and for J ( J in third_covariant() )
    m_T_big(0, 0) = m_t2 * m_t;
    m_T_big(0, 1) = static_cast<T>(3.0) * m_u * m_t2;
    m_T_big(0, 2) = static_cast<T>(3.0) * m_t * m_u2;
    m_T_big(0, 3) = m_u2 * m_u;
     
    m_T_big(1, 0) = m_t2 * m_s;
    m_T_big(1, 1) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_u * m_s, -m_t2, m_v);
    m_T_big(1, 2) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_u * m_v, -m_u2, m_s);
    m_T_big(1, 3) = m_u2 * m_v;
    
    m_T_big(2, 0) = m_t2 * m_s;
    m_T_big(2, 1) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_s * m_v, -m_s2, m_u);
    m_T_big(2, 2) = pr_product_difference<T>(static_cast<T>(2.0) * m_u, m_s * m_v, -m_v2, m_t);
    m_T_big(2, 3) = m_u2 * m_v;

    m_T_big(3, 0) = m_s2 * m_s;
    m_T_big(3, 1) = static_cast<T>(3.0) * m_v * m_s2;
    m_T_big(3, 2) = static_cast<T>(3.0) * m_s * m_v2;
    m_T_big(3, 3) = m_v2 * m_v;
    
    m_coefs_tilda = m_T_big * m_coefs;
}

//IV part of constructor
template<typename T>
void Solver<T>::second_covariant()
{
    //T coef = std::pow(pr_product_difference<T>(m_t, m_v, m_s, m_u), 2); //we dont use in code
    
    //deltas are elements of H Matrix - second covariant, also we use det(H) in solve (look solve()) to recogbize the type of roots; 

    //m_delta1 = std::fma(m_A, m_C, m_B2); //!!!change to FMA breaks the program, idk why
    m_delta1 = pr_product_difference<T>(m_A, m_C, m_B, m_B); 
    m_delta2 = pr_product_difference<T>(m_A, m_D, m_B, m_C);
    //m_delta3 = std::fma(m_B, m_D, m_C2); //!!!change to FMA breaks the program, idk why
    m_delta3 = pr_product_difference<T>(m_B, m_D, m_C, m_C); 

    //init matrix H and calc Det(H)
    Matrix2T<T> H;
    H << static_cast<T>(2.0) * m_delta1, m_delta2, m_delta2, static_cast<T>(2.0) * m_delta3;
    m_detH   = pr_product_difference<T>(4 * m_delta1, m_delta3, m_delta2, m_delta2);

    //transition matrix for TU coords: use for Hessian and roots [x~, w~]
    m_T_lit << m_t, m_u, m_s, m_v;

    //Calc H~ - matrix H in TU coords
    m_H_tilda = m_T_lit * H * m_T_lit.transpose(); // you may add * coef to get C_line
}
 
//V part of constructor
template<typename T>
void Solver<T>::third_covariant()
{
    /*
    J is 3 covariant
    We need calculate J coefs
    In common matrix form:
        J = [t u] * H * [s v].transpose()
    But we will write as in the article Blinn (5 part) with small transformations, since these are identical expressions
    */
    
    //AC - B^2 = m_delta1
    //AD - BC  = m_delta2
    //BD - C^2 = m_delta3
    
    //Aj = (A^2) * D - 3ABC + 2 * (B^3) = A^2 * D - ABC + 2(B^3) - 2ABC = A(AD - BC) - 2B(AC - B^2) =
    //= A * delta2 - 2B * delta1; 
    m_Aj = pr_product_difference<T>(m_A, m_delta2, static_cast<T>(2.0) * m_B, m_delta1);
    
    //Bj = (B^2) * C - 2A(C^2) + ABD = (B^2) * C - A(C^2) + ABD - A(C^2) = -C(AC - B^2) + A(BD - C^2) = 
    //= A * m_delta3 - C * m_delta1
    m_Bj = pr_product_difference<T>(m_A, m_delta3, m_C, m_delta1);
    
    //Cj = 2D(B^2) - B(C^2) - ACD = D(B^2) - B(C^2) +  D(B^2) - ACD = B(BD - C^2) - D(AC - B^2) =
    //= B * m_delta3 - D * m_delta1
    m_Cj = pr_product_difference<T>(m_B, m_delta3, m_D, m_delta1);
    
    //Dj = 3BCD - 2(C^3) - A(D^2) = 2BCD - 2(C^3) + BCD - A(D^2) = 2C(BD - C^2) - D(AD - BC) =
    // = 2C * delta3 - D * delta2
    m_Dj = pr_product_difference<T>(static_cast<T>(2.0) * m_C, m_delta3, m_D, m_delta2);

    Vector4T<T> J;
    J << m_Aj, m_Bj, m_Cj, m_Dj;
    
    //TJ is vector J in TU coords;
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
    //Coefs_line are coefs after depression and scaling  
    T A_line = m_coefs_tilda[0];
    //T B_line = static_cast<T>(0.0); //we dont use this var in code
    
    //C_line = [(t^2) * H(0,0) + tu * H(0,1) + (u^2) * H(1,1)] * 0.5 = 
    //=[t^2 * 2 * delta1 + tu * delta2 + (u^2) * 2 * delta3] * 0.5 = 
    //= (t^2 * delta1) + (tu * 0.5 * delta2) + ((u^2) * delta3) = 
    //= (t^2 * delta1) + (tu * 0.25 * delta2) + (tu * 0.25 * delta2) + ((u^2) * delta3) =
    //= t(t*delta1 + 0.25u * delta2) + u(u *delta3 + 0.25t * delta2)
    T C_line = m_t * pr_product_difference<T>(m_t, m_delta1, static_cast<T>(-0.25) * m_u, m_delta2) + m_u * pr_product_difference<T>(m_u, m_delta3, static_cast<T>(-0.25) * m_t, m_delta2);
    T D_line = m_TJ[0];
    
    //if det(H) < 0 then we solve case 11_ (one real root and one сomplex conjugate root)
    if(m_detH < static_cast<T>(0.0))
    {
        T r = det_less_zero(A_line, C_line, D_line);
        roots.push_back(r);
        return 1;
    }
    //if det(H) > 0 then we solve case 111(3 real roots)
    //if det(H) == 0 then we solve case 11^2 (1 real root with multiplicity 1 and 1 real root multiplicity 2) 
    //if (det(H) == 0) and (delta1 == delta2 == delta3 == 0) then we solve case 1 (1 real root with multiplicity 3) 
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
    //example of call this function in solve(): det_less_zero(A_line, C_line, D_line)
    //At is A_line
    //Cbar is C_line
    //Dbar is D_line
    
    //T0, T1 are temp vars for calculate p and q.
    T T0 = -std::copysign(std::fabs(At) * std::sqrt(-m_detH), Dbar);
    T T1 = -Dbar + T0;
    
    //p and q are temp vars for calc xt1 - root in TU coord
    T p = std::cbrt(T1 / static_cast<T>(2.0));
    T q = is_equal(T0, T1) ? (-p) : (-Cbar / p);
    
    //xt1 - root in TU coord
    T xt1 = (Cbar <= static_cast<T>(0.0)) ? (p + q) : (-Dbar / (fma(p,p, fma(q,q, Cbar))));
    
    //create vector for multiple with transition matrix T_lit - transition matrix for TU coords;
    Vector2T<T> solution_tmp;
    solution_tmp << xt1, static_cast<T>(1.0);
    auto solution = solution_tmp.transpose() * m_T_lit;
    
    //Return on x/w - because our task find only X-roots. We can do it, because cubic is homogeneous.
    return solution(0,0) / solution(0,1);
}

//L-root - algorithm for t = 1, u = 0
template<typename T>
std::pair<T,T> Solver<T>::L_root()
{
    //Cbar is C_line (for t = 1, u = 0)
    //CbarA: the minus sign that occurs in formulas for the root L in a variable
    T CbarA = -m_delta1;
    
    //DbarA = A * delta_2 -2B * delta1 = mA_j (for t = 1, u = 0)
    T DbarA = m_Aj;
    
    //thetaA is angle for cos and sin.
    T thetaA = std::fabs(std::atan2(m_A * std::sqrt(m_detH), -DbarA)) / static_cast<T>(3.0);
    
    //temp vars for calc root xtL
    T xt1A = 2 * std::sqrt(CbarA) * cos(thetaA);
    T xt3A =  pr_product_difference<T>(-std::sqrt(CbarA * static_cast<T>(3.0)), sin(thetaA), std::sqrt(CbarA), cos(thetaA));
    
    //root in TU coords
    T xtL = (xt1A + xt3A > static_cast<T>(2.0) * m_B)? xt1A : xt3A;

    //return undepressed roots x~ and ~w. ~w need for M-root.
    return std::make_pair(xtL - m_B, m_A);
}

//S-root - algorithm for t = 0, u = 1
template<typename T>
std::pair<T,T> Solver<T>::S_root()
{
    //CbarD is C_line (for t = 0, u = 1)
    T CbarD = m_delta3;
    //DbarD = 2C * delta3 - D * delta2 = Dj (for t = 0, u = 1)
    T DbarD = m_Dj;
    
    //thetaD is angle for cos and sin.
    T thetaD = std::fabs(std::atan2(m_D * std::sqrt(m_detH), -DbarD)) / static_cast<T>(3.0);
    
    //temp vars for calc root XtS
    T xt1D = static_cast<T>(2.0) * std::sqrt(-CbarD) * std::cos(thetaD); 
    
    //xt3D = 2 * sqrt(-CbarD) * [cos(thetaD) / -2.0 - (sqrt(3.0) / 2.0) *sin(thetaD)] =
    //= sqrt(-CbarD) * [-cos(thetaD) - sqrt(3) * sin(thetaD)] =  -sqrt(-CbarD) * cos(thetaA) - sqrt(-Cbar * 3) * sin(thetaA)
    T xt3D = pr_product_difference<T>(-std::sqrt(-CbarD), std::cos(thetaD), std::sqrt(-CbarD * static_cast<T>(3.0)), std::sin(thetaD)); 
    
    //root in TU coords
    T xtS = (xt1D + xt3D < static_cast<T>(2.0) * m_C)? xt1D : xt3D;

    //return undepressed roots x~ and ~w. ~w need for M-root.
    return std::make_pair(-m_D, xtS + m_C);

}

//M-root - fusion of S-root and L-root. There is property about the sum of the roots in depressed-scaled form.
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
