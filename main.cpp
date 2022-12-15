#include <iostream>
#include <vector>
#include <cmath> 
#include <excerpt/excerpt.h>
#include <Eigen/Dense>

//---Aliases---
template<typename T>
using Vector4T = Eigen::Matrix<T, 4, 1>;

template<typename T>
using Matrix2T = Eigen::Matrix<T, 2, 2>;

template<typename T>
using Matrix4T = Eigen::Matrix<T, 4, 4>;
//-------------

template <typename T>
class Blinn
{
public:
    Blinn(T A, T B, T C, T D, T t, T u, T k): m_A(A), m_B(B), m_C(C), m_D(D), m_t(t), m_u(u), m_k(k)
    {
        //calc square of params
        m_t2 = std::pow(m_t, 2);
        m_u2 = std::pow(m_u, 2);
        
        //calc s and v
        auto tmp = [&](T X, T Y, T Z) 
                        {return m_k * (pr_product_difference<T>(m_t2, X, -m_u2, Z) + (2 * m_t * m_u * Y));};
        
        m_s = -tmp(B, C, D);
        m_v = tmp(A, B, C);
        
        //calc square of params
        m_s2 = std::pow(m_s, 2);
        m_v2 = std::pow(m_v, 2);
        
        
        //---first covariant---
        first_covariant();
        
        //---second covariant---
        second_covariant();
        
        //---Third covariant---
        third_covariant();
        
    } 
    
    ~Blinn(){}
    
    int solve(std::vector<T> &roots)
    {
        //soon...
        return 0;
    }

private:
    //coefs of polynom
    T m_A;
    T m_B;
    T m_C;
    T m_D;
    
    //params
    T m_t;
    T m_u;
    T m_k;
    T m_s;
    T m_v;
    
    //frequently used degrees
    T m_t2;
    T m_u2;
    T m_s2;
    T m_v2;
    
    //Result of first covariant
    Vector4T<T> m_coefs;
    
    //Result of second covariant
    T m_delta1;
    T m_delta2;
    T m_delta3;
    T m_detH;
    
    T m_delta1_tilda;
    T m_delta2_tilda;
    T m_delta3_tilda;
    T m_detH_tilda;
    
    //Result of third covariant
    T m_Aj;
    T m_Bj;
    T m_Cj;
    T m_Dj;
    
    
    
    //---Functions for Constructor------
    void first_covariant()
    {
        Vector4T<T> coefs;
        coefs << m_A, m_B, m_C, m_D;
        
        Matrix4T<T> T_big;
        
        T_big(0, 0) = std::pow(m_t, 3);
        T_big(0, 1) = static_cast<T>(3.0) * m_u * m_t2;
        T_big(0, 2) = static_cast<T>(3.0) * m_t * m_u2;
        T_big(0, 3) = std::pow(m_u, 3);
        
        T_big(1, 0) = m_t2 * m_s;
        T_big(1, 1) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_u * m_s, -m_t2, m_v);
        T_big(1, 2) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_u * m_v, -m_u2, m_s);
        T_big(1, 3) = m_u2 * m_v;
        
        T_big(2, 0) = m_t2 * m_s;
        T_big(2, 1) = pr_product_difference<T>(static_cast<T>(2.0) * m_t, m_s * m_v, -m_s2, m_u);
        T_big(2, 2) = pr_product_difference<T>(static_cast<T>(2.0) * m_u, m_s * m_v, -m_v2, m_t);
        T_big(2, 3) = m_u2 * m_v;
            
        T_big(3, 0) = std::pow(m_s, 3);
        T_big(3, 1) = static_cast<T>(3.0) * m_v * m_s2;
        T_big(3, 2) = static_cast<T>(3.0) * m_s * m_v2;
        T_big(3, 3) = std::pow(m_v, 3);
        
        m_coefs = T_big * coefs;
    }
    
    void second_covariant()
    {
        T coef = std::pow(pr_product_difference<T>(m_t, m_v, m_s, m_u), 2);
        
        Matrix2T<T> mat_params;
        mat_params << m_t, m_u, m_s, m_v;
        
        //---Init mat_delta---
        m_delta1 = pr_product_difference<T>(m_A, m_C, m_B, m_B);
        m_delta2 = pr_product_difference<T>(m_A, m_D, m_B, m_C);
        m_delta3 = pr_product_difference<T>(m_B, m_D, m_C, m_C);
        m_detH   = pr_product_difference<T>(4 * m_delta1, m_delta3, m_delta2, m_delta2);
        
        Matrix2T<T> mat_delta;
        mat_delta << static_cast<T>(2.0) * m_delta1, m_delta2, m_delta2, static_cast<T>(2.0) * m_delta3;
        
        //Calc delta tilda
        Matrix2T<T> mat_delta_tilda;
        mat_delta_tilda = coef * mat_params * mat_delta * mat_params.transpose();
        
        m_delta1_tilda = mat_delta_tilda(0,0) / static_cast<T>(2.0);
        m_delta2_tilda = mat_delta_tilda(0,1);
        m_delta3_tilda = mat_delta_tilda(1,1) / static_cast<T>(2.0);
        m_detH_tilda = coef * m_detH; 
    }
    
    void third_covariant()
    {
        m_Aj = pr_product_difference<T>(std::pow(m_A, 2), m_D, static_cast<T>(3.0)* m_A, m_B * m_C) 
               + static_cast<T>(2.0) * std::pow(m_B, 3);   
        m_Bj = pr_product_difference<T>(std::pow(m_B, 2), m_C, static_cast<T>(2.0), m_A * std::pow(m_C, 2)) + (m_A * m_B * m_D);
        m_Cj = pr_product_difference<T>(std::pow(m_B, 2), static_cast<T>(2.0) * m_D, m_B, std::pow(m_C, 2)) - (m_A * m_C * m_D);
        m_Dj = pr_product_difference<T>(static_cast<T>(3.0) * m_B, m_C * m_D, static_cast<T>(2.0), std::pow(m_C, 3)) 
            - (m_A * std::pow(m_D, 2));
            
    }
};
////////////////////


template<typename T>
std::pair<T, T> testPolynomial(unsigned int roots_count)
{
    T t = 1;
    T u = 1;
    T max_absolute_error, max_relative_error; 
    
    std::vector<T> roots_computed(roots_count); 
    std::vector<T> roots(roots_count);
    std::vector<T> coefficients(roots_count + 1); 
    generate_polynomial<T>(roots_count, 0, roots_count, 0, 1e-10, -1.0, 1.0, roots, coefficients);
 
    Blinn<T> solver(coefficients[3], coefficients[2], coefficients[1], coefficients[0], t, u, 1);
    
    int cnt_real_roots = solver.solve(roots_computed); // находим число вещ.корней и сами корни
    if (cnt_real_roots != 0 && cnt_real_roots != -1)
    { 
        //если корни найдены и мы не попали в какой-то исключ.случай, то сравниваем найденные корни с истинными
        compare_roots<T>(roots_computed.size(), roots.size(), roots_computed, roots, max_absolute_error, max_relative_error);
         
    } 
    else
    {
        max_absolute_error = 0;
        max_relative_error = 0;
    }
    
    return std::pair<T, T>(max_absolute_error, max_relative_error);
}

int main()
{
    #ifdef FP_FAST_FMA
      std::cout << "FMA enabled\n";
    #else
      std::cout << "FMA disabled\n";
    #endif

    float max_absolut_deviation = 0;
    float max_relative_deviation = 0;
    
    for (auto i = 0; i < 10'000'000; ++i) 
    {
        auto deviation = testPolynomial<float>(3);
        
        if (deviation.first > max_absolut_deviation) 
        {
            max_absolut_deviation = deviation.first;
        }
        if (deviation.second > max_relative_deviation) 
        {
            max_relative_deviation = deviation.second;
        }
    }
    
    std::cout << std::endl << "MAX_ABSOLUT_deviation = " << max_absolut_deviation << std::endl;
    std::cout << std::endl << "MAX_RELATIVE_deviation = " << max_relative_deviation << std::endl;
    
}
