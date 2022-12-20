#include <vector>
#include <Eigen/Dense>

//"Using" keyword in header file is very bad
//I wrote this comment only for your understanding code in solver.cpp
/*---Aliases---
template<typename T>
using Vector2T = Eigen::Matrix<T, 2, 1>;

template<typename T>
using Vector4T = Eigen::Matrix<T, 4, 1>;

template<typename T>
using Matrix2T = Eigen::Matrix<T, 2, 2>;

template<typename T>
using Matrix4T = Eigen::Matrix<T, 4, 4>;
---------------*/

template <typename T>
class Solver
{
public:

    //Constructor
    Solver(const std::vector<T>& coefs);
    
    //Solver
    int solve(std::vector <T>& roots);
    
    //Compare floating-point numbers
    bool is_equal(T a, T b);
private:

    //---PARTS OF CONSTRUCTOR AND VARs---
    //I PART: normalize coefs! B and C divide by 3, all coefs (include B and C) divide by A.
    void normalize_coefs();
    //coefs of polynom
    Eigen::Matrix<T, 4, 1> m_coefs;
    T m_A;
    T m_B;
    T m_C;
    T m_D;

    //II PART: calc params and sqare of params. By default, t = 1 and u = 0;
    void calc_params();
    //params
    T m_t;
    T m_u;
    T m_s;
    T m_v;

    //square of params
    T m_t2;
    T m_u2;
    T m_s2;
    T m_v2;

    //III PART: First Covariant: Calc big Transition matrix
    //New coefs (with tilda) = m_T_big * coefs;
    void first_covariant();
    Eigen::Matrix<T, 4, 4> m_T_big;
    Eigen::Matrix<T, 4, 1> m_coefs_tilda;
    
    //IV PART: Second Covariant: Calc 
    void second_covariant();
    T m_delta1;
    T m_delta2;
    T m_delta3;
    T m_detH;
    
    Eigen::Matrix<T, 2, 2> m_T_lit;
    Eigen::Matrix<T, 2, 2> m_H_tilda;
    
    //V PART: Third covariant
    void third_covariant();
    //Result of third covariant
    T m_Aj;
    T m_Bj;
    T m_Cj;
    T m_Dj;
    Eigen::Matrix<T, 4, 1> m_TJ;
    //-----------------------------
    
    //---Cases OF SOLVING---
    T det_less_zero(T At, T Cbar, T Dbar);
    
    //Roots for det > 0
    std::pair<T,T> L_root();
    std::pair<T,T> S_root();
    std::pair<T,T> M_root(T root00, T root01, T root10, T root11);
    //----------------------
};
