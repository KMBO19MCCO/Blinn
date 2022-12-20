#include <iostream>
#include <vector>
#include <utility>
	
#include <BlinnSolver/solver.hpp>
#include <excerpt/excerpt.h>

template<typename T>
std::pair<T, T> testPolynomial(unsigned int roots_count)
{

    T max_absolute_error, max_relative_error; 
    
    std::vector<T> roots_computed(roots_count); 
    std::vector<T> roots(roots_count);
    std::vector<T> coefficients(roots_count + 1); 
    
    generate_polynomial<T>(roots_count, 0, roots_count, 0, 1e-10, -1.0, 1.0, roots, coefficients);
    
    Solver<T> solver(coefficients);
    
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

    long double max_absolut_deviation = 0;
    long double max_relative_deviation = 0;
    
    for (auto i = 0; i < 1'000'000; ++i) 
    {
         
        auto deviation = testPolynomial<long double>(3);
        
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
